# NYS GHG: Modeling for upscaling greenhouse gas fluxes from wetlands in New York State

The project focuses on using chamber measurement data from 48 weland sites across NYS to build models that can be applied across new geospatial maps of wetlands built by a deep learning model

## Workflow
For the most part, data used here are either specifically GHG data collected for the project or processed by another project repository repsonsible for executing geospatial data processing (`NYS_Wetlands_Data/`). 

Geospatial data should be extracted using the 48 site locations and appended to dataframes in .csv format

Modeling is done in R scripts or Quarto notebooks (.qmd)

TBD: How to handle upscaling the model across the new production maps. Potentially, model upscaling scripts read geospatial map data from NYS_Wetlands_Prod and NYS_Wetlands_Data project folders above this project to use as the model predictors e.g.

$$GHG ~ DEM + CHM + TWI$$

Where DEM, CHM, and TWI are raster layer bands and model predictors.

Climate data (ERA5-Land, PRISM normals) are extracted from Google Earth Engine in `Python_Code_Analysis/GHG_Climate_Extraction.ipynb`.

Modeling findings so far are summarized in `CURRENT_RESULTS.md`.


## Current outline

`Data/` Holds all data 

- `Data/GHG/` holds GHG data observations. Right now it's a bit messy but the main data sources should be `NYDEC_GHG_Sites_6347.gpkg` which contains field-collected predictor data but **not** geospatial predictor data. The `Latitude` and `Longitude` locations in `NYDEC_GHG_Sites_6347.gpkg` provide location for these observations that can be used to extract geospatial predictors.
- `NYDEC_GHG_Sites_EcoEPA_6347.gpkg` (2026-11-25) is the site file used by the modeling workflow: GPS-corrected locations, flux observations (96 rows = 48 sites x 2 rounds), and `Ecoregion_EPA`, the EPA Level III ecoregion name. `GHG_Raster_Extraction.R` reads it and writes `NYDEC_GHG_Site_RasterPredictors.csv`, where `Ecoregion` is the EPA name, `US_L3CODE` its code, and `Ecoregion_NYDEC` the original sampling-stratum label (one-to-one with the four sampled EPA ecoregions).
- `Data/Ecoregions/NY_ECO_L3_6347.gpkg` holds the EPA Level III ecoregion polygons for NY (9 ecoregions; 4 sampled). These zone the statewide upscaling.
- Other iterations of ``NYDEC_GHG_Sites_6347*` have geospatial data extracted already, with the modeling dataframe output as `NYDEC_GHG_Model_Data.csv`. However these should be re-processed with new predictor extraction strategies (See TBD)

`R_Code_Analysis/` holds modeling analysis scripts, and currently (11/25/2026) 

- `GHG_Processing.R` – processes, explores, and creates dataframes for modeling
- `GHG_Raster_Extraction.R` – extracts 10/30/50 m window raster predictors at the sites
- `GHG_Raster_EDA.qmd` – predictor collinearity and univariate flux relationships
- `GHG_Dredge.qmd` – all-subsets model selection on site means, with permutation tests and nested CV
- `GHG_Model_SiteRanef.qmd` – visit-level models with a site random effect
- `GHG_Model_SiteMeans.qmd` – ecoregion-only vs. complex CH4 models for upscaling
- `Archive/GHG_Modeling.R` – earlier modeling script

## TBD

- Python notebook using Google Earth Engine to re-extract the climate and weather data from GHG point locations
- R script to extract re-extract raster predictors from `NYS_Wetlands_Data/` and potentially wetland probability from `NYS_Wetlands_Prod` inference maps made by model version `p027`
    - Need to decide how to extract data since 1m resolution pixel extraction does not capture representative conditions in the wetland sample location
    - Idea: extract using a buffer or moving window approach where points are the center pixel and raster predictor metrics are extracted with functions for the surrounding pixels. This would translate to model prediction using a moving window on each wetland pixel from the maps in `NYS_Wetlands_Prod`
- R modeling script needs to decide on model structure and how to capture auto-correlation
    - can explore autocorrelation in `GHG_Processing.R` 




