library(tidyverse)
library(sf)
library(terra)
library(corrplot)
library(lme4)
library(lmerTest)
library(MuMIn)
library(mgcv)
library(DHARMa)

# ---------------------------------------------------------------------------
# 0. Data prep
# ---------------------------------------------------------------------------
# NOTE: many of the field/soil columns come out of the gpkg as *character*
# ("Top10cm_perc_OM", water table, soil temp, ...). na.omit() therefore does
# not drop their missing values and they are silently unusable as predictors.
# Coerce them explicitly.
chr_to_num <- c(
  "Amb_CH4",
  "Amb_CO2",
  "pw_NH3_N_mg_L",
  "surf_NH3_N_mg_L",
  "pw_NO3_2_N_mg_L",
  "surf_NO3_2_N_mg_L",
  "pw_SRP_mg_L",
  "surf_SRP_mg_L",
  "pw_Acidified_SRP_mg_L",
  "Top6cm_perc_OM",
  "Top6cm_DryWeight_g_cc",
  "Top10cm_perc_OM",
  "Top10cm_DryWeight_g_cc",
  "Top20cm_perc_OM",
  "Top20cm_DryWeight_g_cc",
  "mean_pH",
  "mean_soil_temp",
  "mean_water_temp",
  "mean_DO_percent",
  "average_depth_to_water_table",
  "average_depth_to_water_table_over_20",
  "fraction_standing_water",
  "Water_Table_Category_score_average",
  "mean_woody_cover_abund",
  "mean_graminoids_abund",
  "mean_forbs_abund",
  "mean_non_vascular_abund",
  "mean_bare_ground_abund",
  "species_richness",
  "growth_habit_richness",
  "mean_veg_height"
)

d_raw <- st_read(
  "Data/FieldData/GHG/NYDEC_GHG_Sites_6347_ExtractedPredictors.gpkg",
  quiet = TRUE
) |>
  as_tibble() |>
  mutate(across(any_of(chr_to_num), ~ suppressWarnings(as.numeric(.x))))

# --- response construction -------------------------------------------------
# CH4_avg is strongly right skewed with 12 small negative values
# (min = -0.0177, max = 13.59). Two usable responses:
#
#   CH4_pos   = CH4_avg + SHIFT   strictly positive -> Gamma / lognormal GLMM
#   CH4_r5    = sign(y)*|y|^(1/5) the true (odd) 5th root; keeps the negatives
#
# The old code had `CH4_avg_plus = CH4_avg + -0.02`, which SUBTRACTS 0.02 and
# leaves 27 values <= 0, and `(CH4_avg)^(1/5)`, which returns NaN for every
# negative value (R will not take a fractional power of a negative number).
SHIFT <- 0.02 # > |min(CH4_avg)|; min(CH4_pos) becomes 0.00234

d <- d_raw |>
  mutate(
    CH4_pos = CH4_avg + SHIFT,
    CH4_log = log(CH4_pos),
    CH4_r5 = sign(CH4_avg) * abs(CH4_avg)^(1 / 5), # odd root, negatives kept
    CH4_slog = sign(CH4_avg) * log1p(abs(CH4_avg)),
    Site_ID = as.factor(Site_ID),
    ndvi = (nir - r) / (nir + r),
    ndwi = (g - nir) / (g + nir),
    ndvi_lo = (nir_lo - r_lo) / (nir_lo + r_lo),
    ndwi_lo = (g_lo - nir_lo) / (g_lo + nir_lo)
  )

# Geomorphon classes 3, 4, 7, 10 have 1-4 observations each. Left as-is they
# make dredge produce unstable / singular fits and blow up the SEs. Lump
# anything with < 5 obs into "other" and use flats (1) as the reference.
# Geomorphons categories: 1: Flat, 2: Peak, 3: Ridge, 4: Shoulder, 5: Spur, 6: Slope, 7: Hollow, 8: Footslope, 9: Valley, 10: Pit,
# Revised the grouping categories to include Peak+Ridge and Slope+Hollow based on geomorphology similarity
geo_n <- table(d$Geomorph_local)
geo_rare <- names(geo_n)[geo_n < 3]
geo_pkrd <- c("3", "4")
geo_slphol <- c("6", "7")
d <- d |>
  mutate(
    Geomorph_local = case_when(
      as.character(Geomorph_local) %in% geo_pkrd ~ "PeakRidge",
      as.character(Geomorph_local) %in% geo_slphol ~ "SlopeHollow",
      .default = as.character(Geomorph_local)
    ),
    Geomorph_local = relevel(factor(Geomorph_local), ref = "1")
  )
table(d$Geomorph_local)

# Coarse 2-level version. Class 9 (valley) is the only geomorphon level that is
# ever significant, and dropping to 2 levels is what makes Geomorph x continuous
# interactions estimable -- see section 1b.
d <- d |>
  mutate(
    Geo2 = factor(
      ifelse(Geomorph_local == "9", "valley", "other"),
      levels = c("other", "valley")
    )
  )

# Geomorphon classes 3, 4, 7, 10 have 1-4 observations each. Left as-is they
# make dredge produce unstable / singular fits and blow up the SEs. Lump
# anything with < 5 obs into "other" and use flats (1) as the reference.
# Geomorphons categories: 1: Flat, 2: Peak, 3: Ridge, 4: Shoulder, 5: Spur, 6: Slope, 7: Hollow, 8: Footslope, 9: Valley, 10: Pit,
# Revised the grouping categories to include Peak+Ridge and Slope+Hollow based on geomorphology similarity
geo_n <- table(d$Geomorph_local)
geo_rare <- names(geo_n)[geo_n < 3]
geo_pkrd <- c("3", "4")
geo_slphol <- c("6", "7")
d <- d |>
  mutate(
    Geomorph_local = case_when(
      as.character(Geomorph_local) %in% geo_pkrd ~ "PeakRidge",
      as.character(Geomorph_local) %in% geo_slphol ~ "SlopeHollow",
      .default = as.character(Geomorph_local)
    ),
    Geomorph_local = relevel(factor(Geomorph_local), ref = "1")
  )
table(d$Geomorph_local)

# Coarse 2-level version. Class 9 (valley) is the only geomorphon level that is
# ever significant, and dropping to 2 levels is what makes Geomorph x continuous
# interactions estimable -- see section 1b.
d <- d |>
  mutate(
    Geo2 = factor(
      ifelse(Geomorph_local == "9", "valley", "other"),
      levels = c("other", "valley")
    )
  )

# Raster-derived predictors (the ones available for wall-to-wall inference)
rast_vars <- c(
  "DEM",
  "slope_local",
  "TPI_local",
  "meanc_local",
  "dmv_local",
  "flowacc",
  "twi",
  "CHM",
  "pct_below_1m",
  "pct_1m_to_5m",
  "pct_above_5m",
  "ndvi",
  "ndwi",
  "ndvi_lo",
  "ndwi_lo",
  "r_lo",
  "g_lo",
  "b_lo",
  "nir_lo"
)
# Field/soil predictors (NOT mappable, but see section 4)
soil_vars <- c(
  "Top10cm_DryWeight_g_cc",
  "Top10cm_perc_OM",
  "Water_Table_Category_score_average",
  "fraction_standing_water"
)

# Subset columns BEFORE na.omit so unrelated missing columns don't cost rows
d_mod <- d |>
  select(
    Site_ID,
    Site_ID_round,
    Round_Number,
    CH4_avg,
    CH4_pos,
    CH4_log,
    CH4_r5,
    CH4_slog,
    CO2_avg,
    Geomorph_local,
    Geo2,
    all_of(rast_vars)
  ) |>
  na.omit() |>
  droplevels()

d_soil <- d |>
  select(
    Site_ID,
    Site_ID_round,
    Round_Number,
    CH4_avg,
    CH4_pos,
    CH4_log,
    CH4_r5,
    CH4_slog,
    Geomorph_local,
    Geo2,
    all_of(rast_vars),
    all_of(soil_vars)
  ) |>
  na.omit() |>
  droplevels()

cat("n raster-only:", nrow(d_mod), " sites:", nlevels(d_mod$Site_ID), "\n")
cat("n with soil:  ", nrow(d_soil), " sites:", nlevels(d_soil$Site_ID), "\n")


d_clim <- read.csv("Data/FieldData/GHG/NYDEC_GHG_Sites_ERA5.csv") |>
  mutate(
    Site_ID = as.factor(Site_ID),
    Site_ID_round = as.factor(paste0(Site_ID, "_R", Round_Numb))
  ) |>
  select(Site_ID_round, ends_with("7d"), bio01_C, bio12_mm)
glimpse(d_clim)
d_clim$Site_ID_round %in% d$Site_ID_round


d_full <- left_join(d_mod, d_clim, by = join_by(Site_ID_round))
glimpse(d_full)

#########
extra_stats <- function(m) {
  r2 <- tryCatch(
    suppressWarnings(MuMIn::r.squaredGLMM(m)),
    error = function(e) matrix(NA_real_, 1, 2)
  )
  # r.squaredGLMM returns >1 row for some GLMM families (delta / lognormal /
  # trigamma). Prefer the "delta" row when present.
  i <- if (!is.null(rownames(r2)) && "delta" %in% rownames(r2)) "delta" else 1L
  sing <- tryCatch(as.numeric(lme4::isSingular(m)), error = function(e) {
    NA_real_
  })
  msg <- tryCatch(
    length(m@optinfo$conv$lme4$messages),
    error = function(e) NA_real_
  )
  # Rank deficiency matters more than singularity for factor x continuous
  # terms: lme4 silently DROPS columns and fits a different model than the one
  # you wrote, with only a warning. Flag it explicitly.
  rankdef <- tryCatch(
    {
      X <- lme4::getME(m, "X")
      ncol(X) - qr(X)$rank
    },
    error = function(e) NA_real_
  )
  c(
    R2m = unname(r2[i, 1]),
    R2c = unname(r2[i, 2]),
    singular = sing,
    conv_warn = as.numeric(msg),
    rankdef = as.numeric(rankdef)
  )
}

########

mod_full <- lmer(
  CO2_avg ~
    Geomorph_local +
    (pct_below_1m + ndvi_lo + dmv_local + bio12_mm + bio01_C)^2 +
    (1 | Site_ID),
  data = d_full,
  REML = FALSE,
  na.action = "na.fail"
)
summary(mod_full)

dredge_mods <- MuMIn::dredge(
  mod_full,
  rank = "AICc",
  m.lim = c(1, 8),
  extra = extra_stats
)

head(dredge_mods, 10) # now carries R2m, R2c, singular, conv_warn

# Keep only clean fits, then re-rank. Use MuMIn's own subset method with
# recalc.delta/recalc.weights rather than assigning to the columns directly --
# `$<-` strips the "model.selection" class and breaks get.models() downstream.
dredge_ok <- subset(
  dredge_mods,
  extra_stats.singular == 0 &
    extra_stats.conv_warn == 0 &
    extra_stats.rankdef == 0,
  recalc.weights = TRUE,
  recalc.delta = TRUE
)

cat(
  "\ndredged:",
  nrow(dredge_mods),
  " singular/non-converged dropped:",
  nrow(dredge_mods) - nrow(dredge_ok),
  "\n"
)
head(dredge_ok, 10)

# Best clean model, and the best clean model that also clears an R2m threshold
best <- MuMIn::get.models(dredge_ok, subset = 1)[[1]]
summary(best)
MuMIn::r.squaredGLMM(best)

# Best clean model that also clears an R2m threshold -- empty for raster-only
subset(dredge_ok, extra_stats.R2m > 0.4) |> head()

formula(best)

mod_test <- lmer(
  CH4_log ~
    Geomorph_local + pct_below_1m * dmv_local * bio01_C + (1 | Site_ID),
  data = d_full,
  REML = FALSE,
  na.action = "na.fail"
)
summary(mod_test)
MuMIn::r.squaredGLMM(mod_test)

##### model Co2

mod_co2 <- lmer(
  CO2_avg ~
    Geomorph_local + ndvi * dmv_local * bio01_C + (1 | Site_ID),
  data = d_full,
  REML = FALSE,
  na.action = "na.fail"
)
summary(mod_co2)
MuMIn::r.squaredGLMM(mod_co2)
