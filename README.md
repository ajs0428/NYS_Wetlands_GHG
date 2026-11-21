# NYS GHG: Modeling for upscaling greenhouse gas fluxes from wetlands in New York State

The project focuses on using chamber measurement data from 48 weland sites across NYS to build models that can be applied across new geospatial maps of wetlands built by a deep learning model

## Workflow
For the most part, data used here are either specifically GHG data collected for the project or processed by another project repository repsonsible for executing geospatial data processing (`NYS_Wetlands_Data/`). 

Geospatial data should be extracted using the 48 site locations and appended to dataframes in .csv format

Modeling is done in R scripts or Quarto notebooks (.qmd)

TBD: How to handle upscaling the model across the new production maps. Potentially, model upscaling scripts read geospatial map data from NYS_Wetlands_Prod and NYS_Wetlands_Data project folders above this project to use as the model predictors e.g.

$$GHG ~ DEM + CHM + TWI$$

Where DEM, CHM, and TWI are raster layer bands and model predictors.

Other data sources have come by extracting climate data from Google Earth Engine but haven't been super useful in models. 


## Current outline

`Data/` Holds all data 

- `Data/GHG/` holds GHG data observations. Right now it's a bit messy but the main data sources should be `NYDEC_GHG_Sites_6347.gpkg` which contains field-collected predictor data but **not** geospatial predictor data. The `Latitude` and `Longitude` locations in `NYDEC_GHG_Sites_6347.gpkg` provide location for these observations that can be used to extract geospatial predictors.
- Other iterations of ``NYDEC_GHG_Sites_6347*` have geospatial data extracted already, with the modeling dataframe output as `NYDEC_GHG_Model_Data.csv`. However these should be re-processed with new predictor extraction strategies (See TBD)

`R_Code_Analysis/` holds modeling analysis scripts, and currently (11/17/2026) 

- `GHG_Processing.R` – processes, explores, and creates dataframes for modeling
- `GHG_Modeling.R` - models GHG (CO2 and CH4) as a function of geospatial variables

## TBD

- Python notebook using Google Earth Engine to re-extract the climate and weather data from GHG point locations
- R script to extract re-extract raster predictors from `NYS_Wetlands_Data/` and potentially wetland probability from `NYS_Wetlands_Prod` inference maps made by model version `p027`
    - Need to decide how to extract data since 1m resolution pixel extraction does not capture representative conditions in the wetland sample location
    - Idea: extract using a buffer or moving window approach where points are the center pixel and raster predictor metrics are extracted with functions for the surrounding pixels. This would translate to model prediction using a moving window on each wetland pixel from the maps in `NYS_Wetlands_Prod`
- R modeling script needs to decide on model structure and how to capture auto-correlation
    - can explore autocorrelation in `GHG_Processing.R` 




