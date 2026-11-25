### GHG_Raster_Extraction.R
###
### Extracts block-summary raster predictors at the 48 GHG chamber sites from
### the per-HUC12 1 m sources in NYS_Wetlands_Data (terrain, hydro, HAND, CHM,
### NAIP, leaf-off ortho) and the statewide wetland probability map in
### NYS_Wetlands_Prod.
###
### Window definition -- built to match a 10 m prediction map exactly:
###   * The 10 m grid is each HUC's 1 m DEM grid aggregated by 10 from the DEM
###     origin (every source, incl. the Prod probability maps, shares that grid).
###     A 10 m map made with terra::aggregate(stack, fact = 10) per HUC has
###     these same cells.
###   * A site is assigned the 10 m cell its point falls in. The point is
###     usually not at the cell centre; `offset_m` records the distance.
###   * A W m window is the W/10 x W/10 block of 10 m cells centred on that
###     cell (W must be an odd multiple of 10). Its mean and SD are those of
###     all valid 1 m pixels inside it. On the map that is
###     aggregate(sum, sum of squares, count) -> focal(n x n, sum), so a
###     prediction pipeline reproduces these values exactly.
###   * SD is the population SD of the 1 m pixels: sqrt(E[x^2] - E[x]^2).
###   * ndvi_lo / ndwi_lo are derived per 1 m pixel from the ortho with
###     huc_stack.R's ortho_lo_indices(), then summarized like other bands.
###   * Windows that reach past the HUC raster edge are padded with NA; the
###     per-window `validfrac_*` column reports the share of 1 m pixels used.
###   * Geomorphons (Geomorph_local, categorical classes 1-10) are summarized
###     as the share of the window's valid 1 m pixels in four groups that sum
###     to 1 (GEO_<group>_frac_<W>m), plus the majority group (GEO_mode_<W>m;
###     ties go to the group listed first in GEO_GROUPS) and its share
###     (GEO_modefrac_<W>m). The groups keep every class in a group the sites
###     cover, so a map built from them never meets an unsampled level. A
###     statewide tally (2026-11, wet cells p >= 0.5) found ~47% flat, ~27%
###     concave, ~21% footslope and ~5% convex majority cells, against 22 /
###     15 / 7 / 4 sites (10 m majority). On the map: per-class count via aggregate(sum of
###     indicator) -> focal(n x n, sum) / valid count, same as the means.
###
### Raster predictors are static, so output is one row per site (Site_ID), to
### be joined to the per-round GHG data in GHG_Processing.R.
###
### Ecoregion: `Ecoregion` is the EPA Level III name (US_L3NAME) from the
### sites file, checked against a spatial join to the EPA polygons;
### `US_L3CODE` is its code. The original NYDEC sampling-stratum label is kept
### as `Ecoregion_NYDEC` (one-to-one with the four sampled EPA ecoregions).
###
### Run from the project root with an R that has terra + sf (R 4.4.3 on BioHPC):
###   /programs/R-4.4.3/bin/Rscript R_Code_Analysis/GHG_Raster_Extraction.R

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
})

## --- Settings ---------------------------------------------------------------
DATA_ROOT <- "/ibstorage/anthony/NYS_Wetlands_Data"
PROD_PROBS <- "/ibstorage/anthony/NYS_Wetlands_Prod/Data/HUC_DL_Predictions_Prod/p027_focal_bin_p21/ensemble_cal"
## Sites with GPS corrections and EPA Level III ecoregion labels (Ecoregion_EPA)
SITES <- "Data/GHG/NYDEC_GHG_Sites_EcoEPA_6347.gpkg"
## EPA Level III ecoregions clipped to NY; used to check the site labels and add
## the US_L3CODE that the upscaling step joins on
ECOREGIONS <- "Data/Ecoregions/NY_ECO_L3_6347.gpkg"
HUC_ZONES <- file.path(
  DATA_ROOT,
  "Data/NY_HUCS/NY_Cluster_Zones_250_CROP_NAomit_6347.gpkg"
)
OUT_CSV <- "Data/GHG/NYDEC_GHG_Site_RasterPredictors.csv"

CELL <- 10 # map resolution (m) = aggregation factor on the 1 m grid
WINDOWS <- c(10, 30, 50) # window widths (m); odd multiples of CELL

## Curated predictors: band name in the huc_stack.R layers -> output prefix
BANDS <- c(
  DEM = "DEM",
  twi = "twi",
  hand = "HAND",
  CHM = "CHM",
  TPI_local = "TPI",
  dmv_local = "DMV",
  ndvi = "NDVI",
  ndwi = "NDWI",
  ndvi_lo = "NDVI_lo",
  ndwi_lo = "NDWI_lo",
  WET_prob = "WETprob"
)

## Geomorphon classes (rgeomorphon codes) -> groups. Hollow and depression are
## water-receiving like valley; the sites have none of either as a majority
## class, so they ride with valley. Order sets the GEO_mode tie-break.
GEO_GROUPS <- list(
  flat = 1, # flat
  concave = c(7, 9, 10), # hollow, valley, depression
  footslope = 8, # footslope
  convex = 2:6 # peak, ridge, shoulder, spur, slope
)
stopifnot(
  setequal(unlist(GEO_GROUPS), 1:10),
  !anyDuplicated(unlist(GEO_GROUPS))
)

stopifnot(all(WINDOWS %% CELL == 0), all((WINDOWS / CELL) %% 2 == 1))

## --- Reuse the band contract from NYS_Wetlands_Data -------------------------
source(file.path(DATA_ROOT, "R_Code_Analysis/huc_stack.R"))
abs_dirs <- lapply(huc_source_dirs(), \(d) file.path(DATA_ROOT, d))

## --- Sites -> HUC12 / cluster -----------------------------------------------
sites <- st_read(SITES, quiet = TRUE) |>
  group_by(Site_ID) |>
  slice(1) |> # both rounds share coordinates
  ungroup() |>
  select(Site_ID, Ecoregion = Ecoregion_EPA, Ecoregion_NYDEC = Ecoregion, Wetland_Type)

## Check the EPA labels against the polygons and add the L3 code
eco <- st_read(ECOREGIONS, quiet = TRUE) |>
  select(US_L3CODE, US_L3NAME) |>
  st_transform(st_crs(sites))
sites <- st_join(sites, eco, join = st_within)
if (anyDuplicated(sites$Site_ID) || anyNA(sites$US_L3NAME) ||
    any(sites$US_L3NAME != sites$Ecoregion)) {
  stop("Site Ecoregion_EPA labels don't match a spatial join to ", ECOREGIONS)
}
sites <- sites |> select(-US_L3NAME) |> relocate(US_L3CODE, .after = Ecoregion)

zones <- st_read(HUC_ZONES, quiet = TRUE) |> select(huc12, cluster)
sites <- st_join(sites, zones, join = st_within)
if (anyDuplicated(sites$Site_ID)) {
  stop("A site falls in more than one HUC zone")
}
if (anyNA(sites$huc12)) {
  stop(
    "Sites outside HUC zones: ",
    paste(sites$Site_ID[is.na(sites$huc12)], collapse = ", ")
  )
}
message(nrow(sites), " sites in ", n_distinct(sites$huc12), " HUC12s")

## --- One site -----------------------------------------------------------------
## Reads the 1 m pixels of the largest window (aligned to the 10 m blocks) for
## every band, then summarizes each nested window.
extract_site <- function(x, y, huc, cluster) {
  paths <- huc_source_paths(huc, cluster, dirs = abs_dirs, profile = "prod")
  if (!huc_sources_ready(paths, huc)) {
    stop("Missing sources for HUC ", huc)
  }
  # Open sources directly rather than via huc_layers(): it log-transforms
  # flowacc over the whole HUC (eager, minutes per HUC), and flowacc is unused.
  pick1 <- \(f) {
    if (length(f) > 1) {
      message("  multiple files, using ", basename(f[1]))
    }
    f[1]
  }
  lyrs <- lapply(paths, \(f) rast(pick1(f)))
  names(lyrs$dem) <- "DEM"
  names(lyrs$hand) <- "hand"
  names(lyrs$ortho) <- paste0(names(lyrs$ortho), "_lo")

  prob_f <- file.path(
    PROD_PROBS,
    sprintf("DLpred_binary_cluster_%s_huc_%s_probs.tif", cluster, huc)
  )
  if (!file.exists(prob_f)) {
    stop("No probability map: ", prob_f)
  }
  prob <- rast(prob_f)[[1]] # band 1 = WET probability (calibrated)
  names(prob) <- "WET_prob"

  # 10 m block containing the point, counted from the DEM origin
  ref <- lyrs$dem
  x0 <- xmin(ref)
  y0 <- ymax(ref)
  col <- floor((x - x0) / CELL)
  row <- floor((y0 - y) / CELL)
  k <- (max(WINDOWS) / CELL - 1) / 2
  e <- ext(
    x0 + (col - k) * CELL,
    x0 + (col + k + 1) * CELL,
    y0 - (row + k + 1) * CELL,
    y0 - (row - k) * CELL
  )

  grab <- function(r) extend(crop(r, e, snap = "near"), e) # NA-pad at raster edge
  ortho <- ortho_lo_indices(grab(lyrs$ortho))
  stk <- c(
    grab(lyrs$dem),
    grab(lyrs$terr),
    grab(lyrs$hydro),
    grab(lyrs$hand),
    grab(lyrs$chm),
    grab(lyrs$naip),
    ortho,
    grab(prob)
  )
  if (!compareGeom(stk, grab(ref), stopOnError = FALSE)) {
    stop("Grid mismatch in HUC ", huc)
  }
  stk <- stk[[names(BANDS)]]

  side <- (2 * k + 1) * CELL # pixels per side of the largest window
  stopifnot(nrow(stk) == side, ncol(stk) == side)
  arr <- as.array(stk) # rows x cols x bands, row 1 = north
  # Geomorphons kept apart from `arr` so validfrac_* is unchanged by them
  geo <- round(as.matrix(grab(lyrs$terr[["Geomorph_local"]]), wide = TRUE))
  stopifnot(nrow(geo) == side, ncol(geo) == side)

  out <- list(
    block_x = x0 + (col + 0.5) * CELL,
    block_y = y0 - (row + 0.5) * CELL
  )
  out$offset_m <- sqrt((x - out$block_x)^2 + (y - out$block_y)^2)

  for (w in WINDOWS) {
    half <- w / 2
    idx <- (side / 2 - half + 1):(side / 2 + half)
    sub <- arr[idx, idx, , drop = FALSE]
    for (b in seq_along(BANDS)) {
      v <- sub[,, b]
      v <- v[is.finite(v)]
      m <- if (length(v)) mean(v) else NA_real_
      out[[sprintf("%s_mean_%dm", BANDS[[b]], w)]] <- m
      out[[sprintf("%s_sd_%dm", BANDS[[b]], w)]] <-
        if (length(v)) sqrt(max(mean(v^2) - m^2, 0)) else NA_real_
    }
    out[[sprintf("validfrac_%dm", w)]] <- mean(apply(sub, c(1, 2), \(p) {
      all(is.finite(p))
    }))

    g <- geo[idx, idx]
    g <- g[is.finite(g) & g %in% 1:10]
    fr <- vapply(
      GEO_GROUPS,
      \(cl) if (length(g)) mean(g %in% cl) else NA_real_,
      0
    )
    for (grp in names(GEO_GROUPS)) {
      out[[sprintf("GEO_%s_frac_%dm", grp, w)]] <- fr[[grp]]
    }
    out[[sprintf("GEO_mode_%dm", w)]] <-
      if (length(g)) names(GEO_GROUPS)[which.max(fr)] else NA_character_
    out[[sprintf("GEO_modefrac_%dm", w)]] <- if (length(g)) {
      max(fr)
    } else {
      NA_real_
    }
  }
  as_tibble(out)
}

## --- Run ------------------------------------------------------------------------
xy <- st_coordinates(sites)
res <- lapply(seq_len(nrow(sites)), function(i) {
  message(sprintf(
    "[%2d/%d] %s (HUC %s)",
    i,
    nrow(sites),
    sites$Site_ID[i],
    sites$huc12[i]
  ))
  extract_site(xy[i, 1], xy[i, 2], sites$huc12[i], sites$cluster[i])
})

out <- sites |>
  st_drop_geometry() |>
  mutate(x = xy[, 1], y = xy[, 2]) |>
  bind_cols(bind_rows(res))

write.csv(out, OUT_CSV, row.names = FALSE)
message("Wrote ", OUT_CSV, ": ", nrow(out), " sites x ", ncol(out), " columns")
