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
########################################################################################
set.seed(420)
chip_patch_create <- function(wetland_file){
    setGDALconfig("GDAL_PAM_ENABLED", "FALSE")
    huc_num <- str_extract(wetland_file, "(?<=huc_)\\d+")

    target_wetlands <- st_read(wetland_file, quiet = TRUE) # target wetlands
    tw_centroid <- st_centroid(target_wetlands) |> st_geometry() |> st_cast(to = "MULTIPOINT") #centroid cast to multipoint

    tw_exterior <- st_exterior_ring(target_wetlands) # line border of wetland
    tw_boundary <- st_boundary(tw_exterior) |> st_cast("LINESTRING") # cast to linestring
    perimeter <- st_length(tw_boundary)
    n_points <- pmax(0, as.integer(perimeter / 500)) # Divide the perimeter by 500m, round to nearest integer
    tw_b_line <- st_line_sample(tw_boundary, n = n_points) # put n_points on the border

    tw_bl_point <- st_cast(tw_b_line, "POINT")
    tw_bl_point <- tw_bl_point[!st_is_empty(tw_bl_point)]
    tw_c_point <- st_cast(tw_centroid, "POINT")
    tw_bl_c_cmb <- rbind(
        st_sf(geometry = tw_bl_point),
        st_sf(geometry = tw_c_point)
    )
    tw_bl_c_cmbbuff <- st_buffer(tw_bl_c_cmb, dist = 64, endCapStyle = "SQUARE") |>  # buffer the points
        dplyr::mutate(huc_num = huc_num) |> 
        dplyr::select(huc_num, geometry) |> 
        vect()
    
    dem_rast <- l_dem_cluster[grepl(huc_num, l_dem_cluster)] |> rast()
    set.names(dem_rast, "DEM")
    chm_rast <- l_chm_cluster[grepl(huc_num, l_chm_cluster)] |> rast()
    sat_rast <- l_sat_cluster[grepl(huc_num, l_sat_cluster)] |> rast()
    terr_rast <- l_terr_cluster[grepl(huc_num, l_terr_cluster)] |> rast()
    tw_rast <- target_wetlands |> vect() |> terra::rasterize(y = dem_rast, field = "MOD_CLASS") 
    
    stack <- c(dem_rast, terr_rast, chm_rast, sat_rast , tw_rast)
    for(i in seq_len(nrow(tw_bl_c_cmbbuff)))
    # for(i in seq_len(3))    
        {
        crop(stack, tw_bl_c_cmbbuff[i], mask = TRUE, 
             filename = paste0("Data/R_Patches/", "cluster_", args[1], "_huc_", huc_num, "_patch_", i, ".tif" ), 
             overwrite = TRUE)
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

