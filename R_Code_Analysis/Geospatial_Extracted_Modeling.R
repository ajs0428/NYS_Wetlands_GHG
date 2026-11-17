library(sf)
library(tidyverse)
library(here)

dat_features <- st_read(
  "Data/GHG/NYDEC_GHG_Sites_6347_ExtractedPredictors.gpkg"
) |>
  as_tibble() |>
  select(-c(10:12, 14:80))
dat_predictions <- st_read(
  "Data/GHG/NYDEC_GHG_Sites_6347_ExtractedPredictions.gpkg"
) |>
  as_tibble() |>
  select(-c(1:6, 8:83))

dat_aug <- read.csv("Data/GHG/ERA5_Climate_Extraction_all_site_data.csv") |>
  mutate(Site_ID_round = paste0(Site_ID, "_R", Rnd_Nmb)) |>
  select(Site_ID_round, c(15:26))

names(dat_features)
names(dat_predictions)
names(dat_aug)

names(dat_features) %in% names(dat_predictions)

dat_featpred <- left_join(
  dat_features,
  dat_predictions,
  by = join_by(Site_ID_round)
) |>
  left_join(
    dat_aug,
    by = join_by(Site_ID_round)
  ) |>
  select(-ends_with(".x")) |>
  group_by(Site_ID) |>
  reframe(
    across(where(is.numeric), mean),
    across(where(~ !is.numeric(.x)), first),
    .groups = "drop"
  ) |>
  mutate(
    Geomorph_local = case_when(
      Geomorph_local == 1 ~ "Flat",
      Geomorph_local == 2 ~ "Peak",
      Geomorph_local == 3 ~ "Ridge",
      Geomorph_local == 4 ~ "Shoulder",
      Geomorph_local == 5 ~ "Spur",
      Geomorph_local == 6 ~ "Slope",
      Geomorph_local == 7 ~ "Hollow",
      Geomorph_local == 8 ~ "Footslope",
      Geomorph_local == 9 ~ "Valley",
      Geomorph_local == 10 ~ "Depression",
      .default = "None"
    ),
    DL_class = as.factor(DL_class)
  ) |>
  na.omit()


names(dat_featpred)

cordat <- dat_featpred |>
  select(
    where(is.numeric),
    -Round_Number,
    -r,
    -g,
    -b,
    -nir,
    -ends_with("_lo")
  ) |>
  cor()
corrplot::corrplot(cordat)


dat_featpred |>
  ggplot(aes(x = volumetric_soil_water_layer_1_28d, y = CH4_avg)) +
  #geom_boxplot() +
  geom_jitter()


ch4_mod <- lm(
  CH4_avg ~
    # SSW_prob + EMW_prob + FSW_prob +
    Ecoregion *
    DEM +
    SSW_prob +
    temperature_2m_28d +
    total_precipitation_sum_28d +
    volumetric_soil_water_layer_1_28d,
  data = dat_featpred
)
summary(ch4_mod)

co2_mod <- lm(
  CO2_avg ~ twi + total_precipitation_sum_28d + FSW_prob + EMW_prob,
  data = dat_featpred
)
summary(co2_mod)
