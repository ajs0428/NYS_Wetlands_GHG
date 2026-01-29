library(terra)
library(sf)

clusters <- st_read("Data/NY_HUCS/NY_Cluster_Zones_250_NAomit_6347.gpkg")

pond_files <- list.files("Data/Ponds/data/shape_files/", full.names = TRUE, pattern = ".shp")
pond_polys <- lapply(pond_files, st_read, quiet = TRUE) |> 
    dplyr::bind_rows() |> 
    dplyr::select(Shape_Leng, Shape_Area, geometry) |> 
    st_transform("EPSG:6347")

clusters_with_ponds <- st_intersection(pond_polys[1,], clusters) 

hucs_with_ponds <- clusters |> 
    dplyr::filter(huc12 %in% clusters_with_ponds$huc12)
plot(hucs_with_ponds)

dem_list <- list.files("Data/TerrainProcessed/HUC_DEMs/", full.names = TRUE) |> 
    str_subset("wbt", negate = TRUE)
dem_huc_ponds_list <- dem_list[grepl(paste(hucs_with_ponds$huc12, collapse = "|"), dem_list)]

chm_list <- list.files("Data/CHMs/HUC_CHMs/", full.names = TRUE)
chm_huc_ponds_list <- chm_list[grepl(paste(hucs_with_ponds$huc12, collapse = "|"), chm_list)]


buff_extract_fun <- function(pond_file){
    pond <- st_read(pond_file, quiet = TRUE) |> 
        st_transform("EPSG:6347")
    pond_name <- basename(pond_file) |> str_remove(pattern = ".shp")
    pond_huc <-  st_intersection(pond, clusters)[["huc12"]]
    
    pond_10 <- st_buffer(pond, dist = 10) |> st_difference(y = pond) |> vect()
    pond_25 <- st_buffer(pond, dist = 25) |> st_difference(y = pond) |> vect()
    pond_50 <- st_buffer(pond, dist = 50) |> st_difference(y = pond) |> vect()
    
    pond_dem_file <- dem_list[grepl(paste(pond_huc, collapse = "|"), dem_list)]
    pond_chm_file <- chm_list[grepl(paste(pond_huc, collapse = "|"), chm_list)]
    
    dem_chm <- c(rast(pond_dem_file, win = ext(terra::buffer(pond_50, 5))), 
                 rast(pond_chm_file, win = ext(terra::buffer(pond_50, 5))))
    set.names(dem_chm, c("DEM", "CHM"))
    
    pond_10_mask <- terra::mask(x = dem_chm, mask = pond_10, names = c("DEM_10m", "CHM_10m"))
    pond_25_mask <- terra::mask(x = dem_chm, mask = pond_25, names = c("DEM_25m", "CHM_25m"))
    pond_50_mask <- terra::mask(x = dem_chm, mask = pond_50, names = c("DEM_50m", "CHM_50m"))
    
    fn <- paste0("Data/Ponds/PondBuffers/", "Pond_", pond_name, "_10m25m50m_DEM_CHM.tif")
    writeRaster(c(pond_10_mask, pond_25_mask, pond_50_mask),
                filename = fn, overwrite = TRUE)
    return(NULL)
}

system.time({lapply(pond_files, buff_extract_fun)})


### Metrics to extract
    # mean canopy height: mean
    # max canopy height: max
    # percent over 3m: perc_gt_3m
    # mean elevation: mean (also)
    # percent no canopy: perc_eq_0m
dem_chm_buff_files <- list.files("Data/Ponds/PondBuffers/", pattern = ".tif", full.names = TRUE)

prop_over_3m <- function(x, na.rm = TRUE) {
    if (na.rm) {
        x <- x[!is.na(x)]
    }
    sum(x > 3) / length(x)
}
prop_over_0m <- function(x, na.rm = TRUE) {
    if (na.rm) {
        x <- x[!is.na(x)]
    }
    sum(x == 0) / length(x)
}

dem_chm_metrics_extract_fun <- function(demchmFile){
    demchm <- rast(demchmFile)
    pond_name <- sub("^[^_]+_(.+?)_Polygon.*", "\\1", basename(demchmFile))
    mami <- terra::global(demchm, fun = c("mean", "max"), na.rm = TRUE) |> 
        tibble::rownames_to_column() |> 
        dplyr::rename_with(~ifelse(.x == "rowname", "buffer", .x)) 
    over3 <- terra::global(demchm |> select(starts_with("CHM")), fun = prop_over_3m) |> 
        tibble::rownames_to_column() |> 
        dplyr::rename_with(~ifelse(.x == "rowname", "buffer", .x)) |> 
        dplyr::rename_with(~ifelse(.x == "global", "perc_gt_3m", .x))
    eq_0 <- terra::global(demchm |> select(starts_with("CHM")), fun = prop_over_0m)  |> 
        tibble::rownames_to_column() |> 
        dplyr::rename_with(~ifelse(.x == "rowname", "buffer", .x)) |> 
        dplyr::rename_with(~ifelse(.x == "global", "perc_eq_0m", .x))
    
    full_table <- left_join(mami, over3, by = join_by("buffer")) |> 
        left_join(eq_0, by = join_by("buffer")) |> 
        mutate(pond_name = pond_name) |> 
        select(pond_name, everything())
    return(full_table)
}

system.time({t2 <- lapply(dem_chm_buff_files, dem_chm_metrics_extract_fun) |> dplyr::bind_rows()})
t2

readr::write_csv(t2, "Data/Ponds/PondBuffers/DEM_CHM_Pond_Buffer_Metrics.csv")

for(i in dem_chm_buff_files){
    plot(rast(i))
}
