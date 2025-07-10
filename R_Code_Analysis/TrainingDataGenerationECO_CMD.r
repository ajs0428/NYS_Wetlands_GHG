#!/usr/bin/env Rscript
cat("args should be:", "\n",
    "1) path to full NWI dataset", "\n", 
    "2) path to areas of interest (Ecoregions)", "\n",
    "3) the number of areas to subset", "\n")

args = commandArgs(trailingOnly = TRUE) # arguments are passed from terminal to here
cat("these are the arguments: ", (args), "\n")

# test if there is at least one argument: if not, return an error
if (length(args)<3) {
    stop("At least three arguments must be supplied (input file).n", call.=FALSE)
 }

library(terra)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyterra))
library(future)
library(future.apply)

corenum <-  future::availableCores() -2
stopifnot("Too many cores specified" = corenum <= future::availableCores()[[1]])
plan(multicore, workers = corenum) # number of cores is decided by availableCores() and should work with Slurm scheduler
# this should probably be an argument for bash


set.seed(11)


ny_nwi <- vect(args[1])
ny_areas <- vect(args[2]) |> terra::project(crs(ny_nwi))
num_of_areas <- as.numeric(args[3])

total_num_areas <- nrow(unique(values(ny_areas[, "US_L3CODE"])))[[1]]


if (num_of_areas > total_num_areas) {
    stop(cat("Too many areas specified the total is: ", total_num_areas, "\n"), .call = FALSE)
}

area_ids <- strsplit(unique(values(ny_areas[, "US_L3CODE"]))[1:num_of_areas,], " ")

# The function should take a arguments to subset ecoregions and produce NWI sample points within them
training_pts_func <- function(ids, areas = ny_areas, nwi = ny_nwi) {
    
    # The target area is the single ecoregion from within NY State
    target <- tidyterra::filter(areas, US_L3CODE == as.numeric(ids[[1]]))
    
    # target_name is for saving the files
    target_name <- str_replace_all(unique(target$US_L3NAME), pattern = "[ -]", "_")
    # print(target_name)
    
    # target_area is the area dissolved/aggregated to a single polygon
    target_area <- target |> 
        terra::aggregate()
    target_area$US_L3NAME <- target_name
    #print(names(target_area))
    
    # Reduce/crop the NWI wetlands to within target_area 
    nwi_area_crop <- nwi |> terra::crop(y = target_area)
    # print(nwi_area_crop)
    
    # Remove all wetland types that are estuarine, lacustrine, and other
    nwi_area_filter <- nwi_area_crop |>
        dplyr::filter(WETLAND_TY == "Freshwater Forested/Shrub Wetland" | WETLAND_TY == "Freshwater Emergent Wetland")

    # The number of emergent wetland polygons in NWI
    numEMW <- length(nwi_area_filter[nwi_area_filter$WETLAND_TY == "Freshwater Emergent Wetland"])
    numFSW <- length(nwi_area_filter[nwi_area_filter$WETLAND_TY == "Freshwater Forested/Shrub Wetland"])
    
    # The total number of points generated in an area
    pts <- terra::spatSample(target_area, method = "random", size = 5E5) 
    
    # The NWI wetland points are created by sampling 80% the NWI polygons
    nwi_pts_wet <- nwi_area_filter |>
        terra::buffer(-10) |> # a negative buffer should remove points on the lines
        terra::as.points() |>
        sample(size = 0.80*(numEMW+numFSW))  |>
        dplyr::mutate(MOD_CLASS = case_when(WETLAND_TY == "Freshwater Forested/Shrub Wetland" ~ "FSW",
                                            WETLAND_TY == "Freshwater Emergent Wetland" ~ "EMW",
                                            .default = "Other"),
                      COARSE_CLASS = "WET") |>
        dplyr::select(MOD_CLASS, COARSE_CLASS)
    
    # Upland points are defined as outside the NWI polygons 
        # Might have some comission error/included wetlands, so there are many of these
    nwi_pts_upl <- terra::mask(pts, nwi_area_crop |> buffer(10), inverse = TRUE) |>
        dplyr::mutate(MOD_CLASS = "UPL",
                      COARSE_CLASS = "UPL") |>
        dplyr::select(MOD_CLASS, COARSE_CLASS)

    nwi_pts_all <- rbind(nwi_pts_upl, nwi_pts_wet)
    
    # The number of wetland points to balance the classes a bit
    numFSW_pts <- nrow(nwi_pts_all[nwi_pts_all$MOD_CLASS == "FSW", ])
    numEMW_pts <- nrow(nwi_pts_all[nwi_pts_all$MOD_CLASS == "EMW", ])
    
    # If there are fewer than half of emergent vs. forested/scrub/shrub then supplement the points by sampling
        # additional emergent polygons
    if(numEMW_pts < 0.5*numFSW_pts){
        suppPoints <- nwi_area_filter |>
            dplyr::filter(str_detect(WETLAND_TY, "Emergent")) |>
            terra::buffer(-10) |> # a negative buffer should remove points on the lines
            terra::as.points() |>
            sample(size = c(0.5*numFSW_pts, 1000)[which.max(c(0.5*numFSW_pts, 1000))])  |>
            dplyr::mutate(MOD_CLASS = "EMW",
                          COARSE_CLASS = "WET") |>
            dplyr::select(MOD_CLASS, COARSE_CLASS)
        
        nwi_pts_all_supp <- rbind(nwi_pts_all, suppPoints)
    } else { #don't change if > half of forested/scrub/shrub
        nwi_pts_all_supp <- nwi_pts_all
    }
    
    # Summary of point distribution
    print(nwi_pts_all_supp |> as.data.frame() |>  dplyr::group_by(MOD_CLASS) |> dplyr::summarise(count = n()))

    # Saving files for both the study area and the training data points 
    writeVector(target_area, paste0("Data/NY_Ecoregions/",
                                    target_name,
                                   #"_",
                                   #ids,
                                   "_ECO.gpkg"),
                overwrite = TRUE)
    writeVector(nwi_pts_all_supp, paste0("Data/Training_Data/",
                                         target_name,
                                         #"_",
                                         #ids,
                                         "_training_pts.gpkg"), overwrite = TRUE)
}

system.time({future_lapply(area_ids, training_pts_func, future.seed=TRUE)})

# test_nwi <- vect("Data/NWI/NY_Wetlands_6350.gpkg")
# test_areas <- vect("Data/ny_eco_l4/ny_eco_l4.shp") |> terra::project(crs(test_nwi))
# system.time(training_pts_func(ids = 59, nwi = test_nwi, areas = test_areas))
