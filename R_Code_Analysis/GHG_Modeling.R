library(sf)
library(tidyverse)
library(here)

ghg_dat <- read.csv("Data/GHG/NYDEC_GHG_Model_Data.csv")

ch4_mod <- lm(
  CH4_avg ~
    # SSW_prob + EMW_prob + FSW_prob +
    Ecoregion *
    DEM +
    SSW_prob +
    temperature_2m_28d +
    total_precipitation_sum_28d +
    volumetric_soil_water_layer_1_28d,
  data = ghg_dat
)
summary(ch4_mod)

co2_mod <- lm(
  CO2_avg ~ twi + total_precipitation_sum_28d + FSW_prob + EMW_prob,
  data = ghg_dat
)
summary(co2_mod)
