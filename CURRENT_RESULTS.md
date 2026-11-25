# Current results

Status of the GHG flux modeling as of 2026-11-25. Project setup, data and
workflow are in `README.md`. This file covers only what the analyses found.
Numbers come from renders of the notebooks in `R_Code_Analysis/` on
2026-11-25. Re-render a notebook to update its numbers.

## Bottom line

**Upscale CH4 by EPA Level III ecoregion. Don't use the more complex
site-covariate models to upscale.**

* Ecoregion matches how the study was designed (12 sites in each of 4
  sampled ecoregions), and statewide polygons for it exist
  (`Data/Ecoregions/NY_ECO_L3_6347.gpkg`). It is the only model tested that
  predicts held-out sites at least as well as the overall mean.
* Ecoregion alone is a weak model, though: in-sample R² = 0.12, F-test
  p = 0.12. Its value is that it holds up out of sample, not that it predicts
  well.
* Models with raster and climate covariates fit the 48 sites much better
  (R² 0.25–0.46, lower AICc). But every one of them predicts held-out sites
  **worse than the mean** (negative cross-validated R²), or depends on a
  single site with leverage near 1. Their apparent skill doesn't carry over
  to new sites, so they would be a poor basis for statewide maps.

All CH4 results are for site-mean CH4 (mean of the 2 sampling rounds,
n = 48), modeled as `log(CH4 + c)` with c = 0.031 (4 of the 48 site means
are negative). Predictive R² uses the training-fold mean as its baseline, so
R² ≤ 0 means a model predicts no better than that mean.

## 1. Choosing the upscaling model (`GHG_Model_SiteMeans.qmd`)

### Model comparison (CH4)

| Model | k | dAICc | R² (in-sample) | cvR² LOO | cvR² spatial (25 km) | cvR² leave-ecoregion-out |
|---|---|---|---|---|---|---|
| null | 2 | 6.2 | 0 | 0 | 0 | 0 |
| `Ecoregion` (fixed) | 5 | 7.1 | 0.12 | 0.00 | −0.14 | 0.00 |
| `(1 \| Ecoregion)` | 3 | 8.1 | 0.08 (cond.) | −0.01 | −0.08 | 0.00 |
| `Ecoregion + WETprob × DMV_sd + GEO_modefrac` | 9 | **0.0** | **0.40** | −0.10 | −0.07 | −0.13 |
| same, ecoregion random | 7 | 2.4 | 0.24 marg. / 0.37 cond. | −0.28 | −0.19 | −0.12 |
| covariates only, no ecoregion | 6 | 1.8 | 0.26 | −0.77 | −0.79 | −0.14 |

AICc and in-sample R² favor the complex model. Every cross-validation scheme
says it predicts new sites worse than the ecoregion-only model, and worse
than the mean.

### Ecoregion effects

* Ecoregion means on the model scale range from −2.20 (Eastern Great Lakes
  Lowlands, lowest) to −0.78 (Northeastern Highlands). Back-transformed,
  that is median flux of about 0.08 to 0.43 in `CH4_avg` units.
* There is little between-ecoregion variance. The REML estimate is τ² = 0.20
  against a residual variance of 2.28 (ICC = 0.08). The profile 95% interval
  for τ includes 0, and the LRT against the null gives p = 0.27. τ is
  estimated from only 4 ecoregions, so it is poorly determined.
* Fixed vs. random ecoregion: the random effect shrinks each ecoregion mean
  about 49% of the way toward the grand mean. For wetlands in unsampled
  ecoregions, both versions predict −1.30 on the model scale (about 0.24
  median flux). The random-effect version also supplies a variance for those
  areas, which the fixed version can't.
* The 5 unsampled EPA ecoregions cover 10,891 km², 8.7% of NY land area
  (North Central Appalachians, Erie Drift Plain, Atlantic Coastal Pine
  Barrens, Ridge and Valley, Northern Piedmont).

### Why the complex model is a poor upscaling choice

* **It doesn't hold out of sample.** In-sample R² is 0.40 (adjusted 0.30),
  but cross-validated R² is negative under all three schemes. The fixed
  formula came out of the dredge search, so even these CV numbers are
  optimistic.
* **One site drives it.** ADKG24-AR3235 (WETprob = 0.19, DMV_sd = 0.57, both
  extreme) has leverage 0.98 and Cook's D 2.7. Without that site, the
  WETprob × DMV_sd interaction is no longer significant (p goes from 0.006
  to 0.16). Its LOO cvR² goes from −0.10 to 0.16, so the fitted surface is
  largely the shape of this one site.
* **The covariates partly stand in for ecoregion.** Ecoregion still adds to
  the covariates (F-test p = 0.03), and the covariates add to ecoregion
  (p = 0.003) in-sample. Without ecoregion, the covariates' LOO cvR² falls
  to −0.77.
* **Residual spatial autocorrelation is weak either way.** For the
  ecoregion model, Moran's I = 0.16 (p = 0.11, neighbors within 25 km). The
  complex models show slightly negative I. So spatial dependence doesn't
  decide between the models.

## 2. All-subsets selection (`GHG_Dredge.qmd`)

Dredge was run over a screened pool of 22 predictors (10 m window raster
statistics plus PRISM and ERA5 climate; VIF < 10; pairs with |ρ| ≥ 0.7 kept
apart). It fit 30,153 models of up to 5 terms per gas.

| Response | Best model | In-sample R² | Gain over null (AICc) | Permutation p | Nested LOO cvR² |
|---|---|---|---|---|---|
| CH4 | `DMV_sd + NDVI_lo_sd + era5_tmean_14d` | 0.25 | 6.6 | 0.30 | **−0.75** |
| CO2 | `CHM_mean + NDVI_mean + NDWI_lo_mean + prism_ppt_annual + era5_ppt_14d` | 0.56 | 26.9 | **0.01** | **−0.10** |

* **CH4: no signal beyond chance.** A search over shuffled responses finds a
  best model this good 30% of the time. Once the selection step is
  cross-validated, the chosen models predict much worse than the mean. The
  top set has 14 models within 2 AICc and no clear winner. DMV_sd (the SD of
  microtopographic depth) is in all 14, and in the best model of at least 45
  of the 48 CV folds, which makes it the most consistent CH4 predictor.
* **CO2: a real signal that doesn't predict.** The permutation test passes
  (p = 0.01), and CHM_mean plus precipitation terms recur across folds. But
  nested CV is still negative.
* **Hand-picked 3-way candidate** `DMV_sd * WETprob_mean * era5_tmean_14d`:
  in-sample R² = 0.46, cvR² of 0.02–0.09 across LOO, spatial and
  leave-ecoregion-out CV. The search-adjusted permutation test gives
  p = 0.075: the best of 1,540 screened 3-way models on shuffled data has a
  median R² of 0.39. The same high-leverage site (ADKG24-AR3235, leverage
  0.98) sits under the fitted surface.

## 3. Visit-level models with a site random effect (`GHG_Model_SiteRanef.qmd`)

This notebook keeps both visits (n = 96) with `(1 | Site_ID)`, and splits
weather into a site-mean part and a within-site deviation part. Its render
had not finished when this file was written (lmer dredge plus restricted
permutations). Results are still to be added.

## 4. Exploratory analysis (`GHG_Raster_EDA.qmd`)

* **Univariate CH4 signals are weak.** The strongest are DEM_sd_30m
  (ρ = −0.43), DEM_sd_10m (−0.38), DMV_sd_10m (−0.35) and WETprob_mean_10m
  (+0.33). None survives Benjamini–Hochberg correction (all p_BH ≥ 0.16).
* **CO2 tracks vegetation.** NDVI_mean is ρ ≈ +0.52 to +0.56 across windows
  and NDWI_mean ≈ −0.46 to −0.52, all with p_BH ≤ 0.011. CHM_mean is about
  +0.30.
* **Climate:** CO2 relates to 14-day ERA5 precipitation (ρ = −0.49) and its
  anomaly (−0.51), and to PRISM annual temperature (+0.36 to +0.40). These
  are significant after BH. No climate variable relates to CH4 after BH (best
  is PRISM month precipitation, ρ = 0.29, p_BH = 0.32). Climate has been of
  limited use for the CH4 model.
* **Collinearity:** within a window, about 15 of 22 raster predictors can be
  kept with VIF < 5. Strongly paired variables are WETprob mean/SD
  (ρ = −0.92), CHM mean/SD (0.89), NDVI/NDWI (|ρ| ≈ 0.85), DEM_sd/DMV_sd
  (0.82) and twi_mean/HAND_sd (−0.81). Temperature variables are nearly
  redundant within and across PRISM and ERA5.
* **Extraction quality:** windows are essentially complete (validfrac ≥
  0.999). Points sit at most about 6 m from their 10 m cell center.

## Open questions and next steps

* **Upscaling notebook:** back-transform ecoregion predictions to mean flux
  (smearing or lognormal correction, including τ² for unsampled
  ecoregions). Add a site-level bootstrap for uncertainty, and weight by
  wetland area per ecoregion (Σ WETprob × pixel area from `NYS_Wetlands_Prod`
  p027 maps).
* **Fixed vs. random ecoregion:** point predictions are nearly identical. The
  choice comes down to how uncertainty for the 5 unsampled ecoregions is
  handled.
* **CO2:** a predictive model has not been chosen. The vegetation and
  precipitation signals are real but don't cross-validate yet. An
  ecoregion-only comparison like the CH4 one hasn't been run.
* **Complex model as a sensitivity check:** if it is reported next to the
  upscaled estimate, compare the covariate ranges at the 48 sites with
  statewide wetland pixels first, given how much one site shapes the fit.
