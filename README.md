# NYS GHG: Modeling for upscaling greenhouse gas fluxes from wetlands in New York State

The project focuses on using chamber measurement data from 48 weland sites across NYS to build models that can be applied across new geospatial maps of wetlands built by a deep learning model

## Workflow
For the most part, data used here are either specifically GHG data collected for the project or processed by another project repository repsonsible for executing geospatial data processing (NYS_Wetlands_Data). 

Geospatial data should be extracted using the 48 site locations and appended to dataframes in .csv format

Modeling is done in R scripts or Quarto notebooks (.qmd)

TBD: How to handle upscaling the model across the new production maps. Potentially, model upscaling scripts read geospatial map data from NYS_Wetlands_Prod and NYS_Wetlands_Data project folders above this project to use as the model predictors e.g.

$$GHG ~ DEM + CHM + TWI$$

Where DEM, CHM, and TWI are raster layer bands and model predictors.

Also TBD: Other data sources have come by extracting climate data from Google Earth Engine but haven't been super useful in models. 


## Current outline

`R_Code_Analysis/` holds modeling analysis scripts, and currently (11/17/2026) `Geospatial_Extracted_Modeling.R` is the main modeling script




