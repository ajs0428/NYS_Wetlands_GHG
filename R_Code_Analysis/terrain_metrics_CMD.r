#!/usr/bin/env Rscript

cat("args should be:", "\n",
    "1) path to DEM folder", "\n", 
    "2) an odd integer for the window size and scale", "\n",
    "3) the number of cores for the computer to use", "\n")

args = commandArgs(trailingOnly = TRUE) # arguments are passed from terminal to here
cat("these are the arguments: ", args, "\n")


stopifnot("The number of arguments is less than 3" = length(args) < 3)

options(rgl.useNULL = TRUE) # prevents an rgl warning
library(terra)
library(MultiscaleDTM)
library(future)
library(future.apply)

corenum <- as.numeric(args[3])
stopifnot("Too many cores specified" = corenum <= future::availableCores()[[1]])
plan(multicore, workers = corenum) # number of cores is decided by availableCores() and should work with Slurm scheduler
                        # this should probably be an argument for bash

path_var = args[1] #the first argument should be a file path
win_dim = as.numeric(args[2]) # the second argument should be an odd integer for the window size and scale

stopifnot("The dimensions of the window are even" = win_dim%%2 != 0)

dems <- list.files(args[1], 
                   pattern = ".img$",
                   full.names = TRUE)


terrain_calc_func <- function(dem_path, window = c(win_dim,win_dim)){
    #print(window)
    name <- tools::file_path_sans_ext(basename(dem_path))
    
    dem_rast <- terra::rast(dem_path)
    #print(dem_rast)
    #print(paste0(path_var, "/",  name, "_DMV_15.tif"))
    
    s <- MultiscaleDTM::Pfit(dem_rast,
           w = window,
           metrics = "pslope",
           include_scale = TRUE,
           na.rm = TRUE
           )
    writeRaster(s, filename = paste0(path_var, "/", name, "_slope_", win_dim, ".tif"),
           overwrite = T)
    t <- MultiscaleDTM::DMV(dem_rast,
                        w = window,
                        shape = "rectangle",
                        stand = "range",
                        include_scale = TRUE,
                        na.rm = TRUE
                        )
    writeRaster(t, filename = paste0(path_var, "/", name, "_DMV_", win_dim, ".tif"),
                overwrite = T)
    q <- MultiscaleDTM::Qfit(dem_rast,
                       w = window,
                       na.rm = TRUE,
                       include_scale = TRUE,
                       metrics = c("meanc"))
    writeRaster(q, filename = paste0(path_var, "/", name, "_curve_",  win_dim, ".tif"),
                       overwrite = T)
}

system.time({future_lapply(dems, terrain_calc_func, future.seed = T)})
#system.time({lapply(dems, terrain_calc_func)})

