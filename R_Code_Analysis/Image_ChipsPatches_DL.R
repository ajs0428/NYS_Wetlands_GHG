### Image chips/patches for DL  

library(terra)
library(sf) 
library(tidyverse)
library(tidyterra)
library(future)
library(future.apply)


set.seed(11)

########################################################################################

args <- c(
    208,
    "Data/NY_HUCS/NY_Cluster_Zones_250_NAomit_6347.gpkg",
    "Data/NWI/NY_NWI_6347.gpkg"
)


args = commandArgs(trailingOnly = TRUE) # arguments are passed from terminal to here

cat("these are the arguments: \n", 
    "1) Cluster number for HUC groups:", args[1], "\n", 
    "2) path to the overall zones or study areas :", args[2], "\n",
    "3) Path to the NWI wetlands:", args[3], "\n"
)

########################################################################################
l_dem <- list.files("Data/TerrainProcessed/HUC_DEMs/", pattern = ".tif", full.names = TRUE) 
l_dem_cluster <- l_dem[str_detect(l_dem, paste0("cluster_", args[1])) & !str_detect(l_dem, "wbt")]
l_dem_cluster_nums <- str_extract(l_dem, "(?<=cluster_)\\d+(?=_)") |> unique()

l_chm <- list.files("Data/CHMs/HUC_CHMs/", pattern = ".tif", full.names = TRUE) 
l_chm_cluster <- l_chm[str_detect(l_chm, paste0("cluster_", args[1]))]

l_wet <- list.files("Data/Training_Data/Targeted_Wetlands_For_Field_Validation/", pattern = ".gpkg$", full.names = TRUE) 
l_wet_cluster <- l_wet[str_detect(l_wet, paste0("cluster_", args[1]))]


# Combine the points from the centroid and outer perimeter 
    # buffer and create bboxes 
    # remove overlapping bboxes
    # save as a single vector 
    # crop raster stacks
    # save raster stack chips/patches separately with unique IDs
tw <- st_read(l_wet_cluster[1], quiet = TRUE) # target wetlands
tw_c <- st_centroid(tw[1,]) |> st_geometry() |> st_cast(to = "MULTIPOINT") #centroid cast to multipoint
tw_cb <- st_buffer(tw_c, dist = 128) # buffer the centroid (not needed)
tw_cbb <- st_bbox(tw_cb) |> st_as_sfc() # bounding box of the centroid (not needed)
plot(c(st_geometry(tw[1,]), st_geometry(tw_c), st_geometry(tw_cb)))

tw_e <- st_exterior_ring(tw[1,]) # line border of wetland
tw_b <- st_boundary(tw_e) |> st_cast("LINESTRING") # cast to linestring
tw_bl <- st_line_sample(tw_b, 2) # put 2 points on the border
tw_bl_c <- st_combine(c(tw_bl, tw_c)) # combine with centroid
plot(c(st_geometry(tw_bl_c), st_geometry(tw_e), st_geometry(tw[1,])), lwd = c(1, 1, 3), lty = c(1,1,3), border = c("red", "blue", "purple"))

