### Image chips/patches for DL  

library(terra)
library(sf) 
library(dplyr)
library(tidyr)
library(stringr)
library(tidyterra)
library(future)
library(future.apply)


set.seed(11)

########################################################################################

args <- c(
    64,
    "Data/Training_Data/02_Done_Reviewed_NWI_Data/"
)

args = commandArgs(trailingOnly = TRUE) # arguments are passed from terminal to here

cat("these are the arguments: \n", 
    "1) Cluster number for HUC groups:", args[1], "\n", 
    "2) path the reviewed training data :", args[2]
)


setGDALconfig("GDAL_PAM_ENABLED", "FALSE") # does not create aux.xml files
########################################################################################
l_dem <- list.files("Data/TerrainProcessed/HUC_DEMs/", pattern = ".tif", full.names = TRUE) 
l_dem_cluster <- l_dem[str_detect(l_dem, paste0("cluster_", args[1])) & !str_detect(l_dem, "wbt")]
l_dem_cluster_nums <- str_extract(l_dem, "(?<=cluster_)\\d+(?=_)") |> unique()

l_chm <- list.files("Data/CHMs/HUC_CHMs/", pattern = ".tif", full.names = TRUE) 
l_chm_cluster <- l_chm[str_detect(l_chm, paste0("cluster_", args[1]))]

l_naip <- list.files("Data/NAIP/HUC_NAIP_Processed/", pattern = ".tif", full.names = TRUE) 
l_naip_cluster <- l_naip[str_detect(l_naip, paste0("cluster_", args[1]))]

l_terr <- list.files("Data/TerrainProcessed/HUC_TerrainMetrics/", 
                             pattern = paste0("cluster_", args[1], "_huc"), 
                             full.names = TRUE)
l_terr_cluster <- l_terr[str_detect(l_terr, "local") & !str_detect(l_terr, "10m|1000m")]
l_sat_cluster <- list.files("Data/Satellite/HUC_Processed_NY_Sentinel_Indices/", 
                            pattern = paste0("cluster_", args[1]),
                            full.names = TRUE)

l_wet <- list.files(args[2], pattern = ".gpkg$", full.names = TRUE) 
l_wet_cluster <- l_wet[str_detect(l_wet, paste0("cluster_", args[1]))]

print(l_wet_cluster)
length(l_naip_cluster) == length(l_dem_cluster)
########################################################################################
set.seed(420)
filter_min_distance <- function(points, min_dist = 128) {
    n <- nrow(points)
    keep <- rep(TRUE, n)
    
    for (i in 1:(n - 1)) {
        if (keep[i]) {
            too_close <- st_is_within_distance(points[i, ], points[(i + 1):n, ], dist = min_dist)[[1]]
            if (length(too_close) > 0) {
                keep[(i + too_close)] <- FALSE
            }
        }
    }
    
    points[keep, ]
}
chip_patch_create <- function(wetland_file){
    setGDALconfig("GDAL_PAM_ENABLED", "FALSE")
    huc_num <- str_extract(wetland_file, "(?<=huc_)\\d+")
    huc_poly <- sf::st_read("Data/NY_HUCS/NY_Cluster_Zones_250_NAomit_6347.gpkg", quiet = TRUE,
                                  query = paste0("SELECT * FROM NY_Cluster_Zones_250_NAomit_6347 WHERE huc12 = '", huc_num, "'"))
    
    target_wetlands <- st_read(wetland_file, quiet = TRUE) # target wetlands
    tw_centroid <- st_centroid(target_wetlands) |> st_geometry() |> st_cast(to = "MULTIPOINT") #centroid cast to multipoint
    
    tw_boundary <- st_boundary(target_wetlands) |> st_cast("LINESTRING") # cast to linestring
    perimeter <- st_length(tw_boundary)
    n_points <- pmax(0, as.integer(perimeter / 500)) # Divide the perimeter by 500m, round to nearest integer
    tw_b_line <- st_line_sample(tw_boundary, density = 1/500) # put n_points on the border
    
    tw_bl_point <- st_cast(tw_b_line, "POINT")
    tw_bl_point <- tw_bl_point[!st_is_empty(tw_bl_point)]
    tw_c_point <- st_cast(tw_centroid, "POINT")

    ### upland points
    rand_pts <- st_sample(huc_poly, 100)
    target_wetlands_buffer <- st_buffer(target_wetlands, dist = 250)
    rand_pts_intersect <- st_intersects(rand_pts, target_wetlands_buffer, sparse = FALSE)
    pts_outside_target <- rowSums(rand_pts_intersect) == 0
    upl_pts <- rand_pts[pts_outside_target, ]

    ###combine points
    tw_bl_c_cmb <- rbind(
        st_sf(geometry = tw_bl_point),
        st_sf(geometry = tw_c_point),
        st_sf(geometry = upl_pts)
    )
    tw_bl_c_cmb_f <- filter_min_distance(tw_bl_c_cmb, 128) # filter out points that are too close
    tw_bl_c_cmbbuff <- st_buffer(tw_bl_c_cmb_f, dist = 64, endCapStyle = "SQUARE") |>  # buffer the points
        dplyr::mutate(huc_num = huc_num) |>
        dplyr::select(huc_num, geometry) |>
        vect()
    # writeVector(tw_bl_c_cmbbuff,
    #             paste0("Data/R_Patches_Vector/", "cluster_", args[1], "_huc_", huc_num, ".gpkg" ))

    dem_rast <- l_dem_cluster[grepl(huc_num, l_dem_cluster)] |> rast()
    set.names(dem_rast, "DEM")
    chm_rast <- l_chm_cluster[grepl(huc_num, l_chm_cluster)] |> rast()
    sat_rast <- l_sat_cluster[grepl(huc_num, l_sat_cluster)] |> rast()
    terr_rast <- l_terr_cluster[grepl(huc_num, l_terr_cluster)] |> rast()
    tw_rast <- target_wetlands |> vect() |>
        terra::rasterize(y = dem_rast, field = "MOD_CLASS")

    stack <- c(dem_rast, terr_rast, chm_rast, sat_rast , tw_rast)
    for(i in seq_len(nrow(tw_bl_c_cmbbuff)))
    # for(i in seq_len(3))
        {
        fn <- paste0("Data/R_Patches/", "cluster_", args[1], "_huc_", huc_num, "_patch_", i, ".tif" )
        if(!file.exists(fn)){
            crop(stack, tw_bl_c_cmbbuff[i], mask = TRUE,
                 filename = fn,
                 overwrite = TRUE)
        } else {
            next
        }

    }

    return(NULL)

}

### Parallel

if(future::availableCores() > 16){
    corenum <-  4
} else {
    corenum <-  (future::availableCores())
}
print(corenum)
options(future.globals.maxSize= 32.0 * 1e9)
# plan(multisession, workers = corenum)
plan(future.callr::callr)

future_lapply(l_wet_cluster, chip_patch_create, future.seed = TRUE, 
              future.packages = c("terra", "sf", "dplyr", "tidyr", "stringr"),
              future.globals = TRUE)

### Non-parallel
# system.time({lapply(l_wet_cluster, chip_patch_create)})


### Checks
list_patches <- list.files("Data/R_Patches/", full.names = T)
lp <- lapply(list_patches, FUN = \(x) {rast(x) |> nlyr()}) |> unlist()
lapply(list_patches, FUN = \(x) {rast(x) |> nlyr()}) |> unlist() |> table()

le <- lapply(list_patches, FUN = \(x) {rast(x, lyrs = "MOD_CLASS") |> values() |> is.na() |> all()}) |> unlist() 

list_patches[le == TRUE]
