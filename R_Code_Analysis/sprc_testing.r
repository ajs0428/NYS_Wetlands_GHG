library(terra)

l <- list.files("Data/DEMs/imgs/test/", pattern = ".img$", full.names = TRUE)

z <- sprc(l)

h <- vect("Data/NY_HUCS/NY_HUCS_08_6350.gpkg", query = "SELECT * FROM NY_HUCS_08_6350 WHERE huc12 = '020200010102'")  |> 
              terra::project(crs(z[1]))
#hc <- crop(h, vect(ext(z), crs = crs(z[1])))
#hf <- tidyterra::filter(hc, huc12 == "042900050201")
zc <- crop(z, h)
zc

#plet(c(vect(ext(zc), crs = "EPSG:26918"), hf))


hs <- vect("Data/NY_HUCS/NY_HUCS_08_6350.gpkg")[1:5]
hsnames <- as.list(hs[["huc12"]][[1]])
f <- function(x){
    #print(x)
    print(hs["huc12" == x])  
}

lapply(hsnames, f)




writeFunc <- function(sprc) {
    for (i in 1:length(sprc)){
        if(terra::sources(sprc[i]) != ""){
            base::file.copy(terra::sources(sprc[i]), 
                            to = "/Users/Anthony/Data and Analysis Local/NYS_Wetlands_GHG/Data/DEMs/imgs/test_watershed/",
                            overwrite = TRUE)
        } else {
            name <- stringr::str_extract(tools::file_path_sans_ext(terra::sources(sprc[i])), "[^/]*$")
            writeRaster(sprc[i], filename = paste0("/Users/Anthony/Data and Analysis Local/NYS_Wetlands_GHG/Data/DEMs/imgs/test_watershed/",
                                                   name,
                                                   ".img"),
                        overwrite = TRUE)
        }
    }
}

writeFunc(zc)

ls <- list.files("Data/DEMs/imgs/test_watershed", pattern = "slope.*tif$", full.names = TRUE)
ls
zs <- sprc(ls)
zs

slpm <- terra::mosaic(zs)
writeRaster(slpm, "Data/DEMs/imgs/test_watershed/slp_test_watershed.tif")
