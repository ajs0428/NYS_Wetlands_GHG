### Image chips/patches for DL  

library(terra)
library(sf) 
library(dplyr)
library(tidyr)
library(stringr)
library(tidyterra)
library(readr)
library(future)
library(future.apply)


set.seed(11)

########################################################################################

args <- c(
    64,
    "Data/R_Patches_Vector_Reviewed/", #Path to GIS reviewed wetland vector patches
    128
)

args = commandArgs(trailingOnly = TRUE) # arguments are passed from terminal to here

cat("these are the arguments: \n", 
    "1) Cluster number for HUC groups:", args[1], "\n", 
    "2) path the reviewed training data :", args[2], "\n",
    "3) patch size :", args[3], "\n"
)


setGDALconfig("GDAL_PAM_ENABLED", "FALSE") # does not create aux.xml files but maybe needed
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

logpath <- "Data/R_Patches_Vector/Vector_Patch_Checklist.csv"
########################################################################################
fct_df <- data.frame(ID = 0:4, MOD_CLASS = c("EMW", "FSW", "OWW", "SSW", "UPL"))
patchsize = as.numeric(args[3])
########################################################################################
set.seed(420)

rast_chip_patch_create <- function(wetland_file){
    ## Setup vars
    if(grepl("NWI", basename(wetland_file))){
        sourceWetlands <- "NWI"
    } else if(grepl("NHP", basename(wetland_file))){
        sourceWetlands <- "NHP"
    } else if(grepl("Laba", basename(wetland_file))){
        sourceWetlands <- "Laba"
    } else {
        sourceWetlands <- "Other"
    }
    patchsize = as.numeric(args[3])
    huc_num <- str_extract(wetland_file, "(?<=huc_)\\d+")
    # huc_poly <- sf::st_read("Data/NY_HUCS/NY_Cluster_Zones_250_NAomit_6347.gpkg", quiet = TRUE,
    #                               query = paste0("SELECT * FROM NY_Cluster_Zones_250_NAomit_6347 WHERE huc12 = '", huc_num, "'"))

    dem_rast <- l_dem_cluster[grepl(huc_num, l_dem_cluster)] |> rast()
    set.names(dem_rast, "DEM")
    chm_rast <- l_chm_cluster[grepl(huc_num, l_chm_cluster)] |> rast()
    sat_rast <- l_sat_cluster[grepl(huc_num, l_sat_cluster)] |> rast()
    terr_rast <- l_terr_cluster[grepl(huc_num, l_terr_cluster)] |> rast()
    naip_rast <- l_naip_cluster[grepl(huc_num, l_naip_cluster)] |> rast()
    set.names(naip_rast, c("r", "g", "b", "nir", "n_ndvi", "n_ndwi"))
    
    stack <- c(dem_rast, terr_rast, chm_rast, sat_rast, naip_rast)
    stack_fn <- paste0("Data/HUC_Raster_Stacks/HUC_DL_Stacks/", "cluster_", args[1], "_huc_", huc_num, "_stack.tif")
    if(!file.exists(stack_fn)){
        writeRaster(stack, filename = stack_fn, overwrite = TRUE)
    }
    
    ### Union all the polygons then rejoin and separate as groups 
        ### so that each patch of touching polygons is a separate 
            ### object that can be used to crop the rasters
    tw <- st_read(l_wet_cluster[[1]])
    tw_union <- tw |>
        st_union() |>
        st_cast("POLYGON") |>
        st_as_sf() |>
        mutate(group_id = row_number())
    st_geometry(tw_union) <- "geom"
    tw_union_area <- tw_union |> 
        mutate(area = as.numeric(st_area(geom))) |> 
        filter(area >= ((patchsize*2)**2)-0.5) #remove patches that are smaller than the 256*256 dimensions
    tw_grouped_list <- tw |> st_join(tw_union_area, left = FALSE) |>
        group_split(group_id)
    
    #### Each patch should be a separate file that is patchsize*2 x patchsize*2
    for(i in seq_along(tw_grouped_list)){
        
        tw_vect <- vect(tw_grouped_list[[i]]) 
        dem_crop <- crop(dem_rast, tw_vect, touches = TRUE)
        tw_rast <- tw_vect  |>
            terra::rasterize(y = dem_crop, field = "MOD_CLASS", touches = TRUE)
        tw_rast_lc <- levels(tw_rast)[[1]][[2]] #character vector of levels present
        tw_rast_ln <- levels(tw_rast)[[1]][[1]] #numbers/integers of levels present
        fct_n <- fct_df[fct_df$MOD_CLASS %in% tw_rast_lc, ][,1] # subset the levels present from the full factor dataframe
        tw_rast_sub <- subst(tw_rast, from = tw_rast_ln, to = fct_n, raw = TRUE)
        #tw_rast_sub_int <- terra::as.int(tw_rast_sub)
        levels(tw_rast_sub) <- fct_df
        
        fn <- paste0("Data/R_Patches/", sourceWetlands,"_cluster_", args[1], "_huc_", huc_num, "_patch_", i, "_", patchsize*2, "m.tif" )
        fn_labels <- paste0("Data/R_Patches_Labels/", "labels_only_", sourceWetlands, "_cluster_", args[1], "_huc_", huc_num, "_patch_", i, "_", patchsize*2, "m.tif" )

        # Regular Patches with all predictors
        if(!file.exists(fn)){
            cropped_stack <- crop(stack, tw_vect, mask = TRUE)
            cropped_stack_labeled <- c(cropped_stack, tw_rast_sub)
            writeRaster(cropped_stack_labeled, filename = fn, overwrite = TRUE)
        } else {
            message("Already file ", fn)
        }
        #Labels only patches NO predictors
        if(!file.exists(fn_labels)){
            writeRaster(tw_rast_sub, filename = fn_labels, overwrite = TRUE)
        } else {
            message("Already file ", fn_labels)
        }
    }

    return(NULL)

}



### Parallel

# if(future::availableCores() > 16){
#     corenum <-  4
# } else {
#     corenum <-  (future::availableCores())
# }
# print(corenum)
# options(future.globals.maxSize= 32.0 * 1e9)
# # plan(multisession, workers = corenum)
# plan(future.callr::callr)
# 
# future_lapply(l_wet_cluster, rast_chip_patch_create, 
#               future.seed = TRUE,
#               future.packages = c("terra", "sf", "dplyr", "tidyr", "stringr", "purrr"),
#               future.globals = TRUE)

### Non-parallel
system.time({lapply(l_wet_cluster, rast_chip_patch_create)})


# l_patches <- list.files("Data/R_Patches_Vector")
# 
# check_df <- data.frame(patch_file_name = l_patches,
#                        reviewer = rep("NAME", length(l_patches)),
#                        boundaries_altered = rep("TBD", length(l_patches)),
#                        confidence = rep("TBD", length(l_patches)))
# 
# readr::write_csv(check_df, "Data/R_Patches_Vector/Vector_Patch_Checklist.csv")
# ### Checks
# list_patches <- list.files("Data/R_Patches_Labels/", full.names = T)
# lapply(list_patches, \(x) rast(x))
# lp <- lapply(list_patches, FUN = \(x) {rast(x) |> nlyr()}) |> unlist()
# # lapply(list_patches, FUN = \(x) {rast(x) |> nlyr()}) |> unlist() |> table()
# 
# le <- lapply(list_patches, FUN = \(x) {rast(x, lyrs = "MOD_CLASS") |> values() |> unique() |> nrow()}) |> unlist()
# 
# list_patches[le == 1]
# list_patches[lp < 27]



