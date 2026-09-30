# =============================================================================
# 03_analyses.R
#
# Ecoregion-scale analyses combining the ecoregionalization (script 01) and
# the Chl-a trends (script 02).
#
# Pouly et al. (2026). Multidecadal trends in phytoplankton biomass and
# phenology in Southern Ocean water masses. JGR: Oceans.
#
# Inputs : outputs of scripts 01 and 02
# Outputs: outputs/03_analyses/
#   cell_ecoregion.rds          ecoregion of every 4 km pixel
#   trends_with_ecoregion.rds   trends of scripts 02 + ecoregion
#   table_S1_by_ecoregion.csv   Table S1
#   summary_basin.csv           basin-wide statistics quoted in the text
#   fig4_biomass_by_ecoregion.rds, fig4_sen_by_ecoregion.rds   (Figure 4)
#   relative_change.rds         relative change in integrated Chl-a (Figure S8)
#   climatologies.rds           mean Chl-a and mean fPC1-fPC4 scores (Figures 5, S7, S10)
#   two_period_cycles.rds       two-period comparison (Figures 7, 8, S12, S13)
#   explained_variance.csv      ecoregions vs front-based zones
#   gradient_tests.csv          coincidence of trend gradients with boundaries
# =============================================================================


# ---- 0. Packages, paths and inputs --------------------------------------------

library(dplyr)
library(sf)
library(terra)
library(fda)
library(here)

dir_01  <- "G:/papier_phyto_change_2026/test code/ecoregion"
dir_02  <- "G:/papier_phyto_change_2026/test code/chla_fit"
out_dir <- "G:/papier_phyto_change_2026/test code/chla_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

eco_file  <- file.path(dir_01, "ecoregions", "ecoregions_k8.shp")
eco_field <- "ecoregion"          # "zone" if you use your original layer_k8.shp
sd_file   <- file.path(dir_01, "ecoregions", "cluster_sd.tif")

n_steps       <- 23
min_seasons   <- 10
early_window  <- c(2010, 2015)    # seasons 1998/99 to 2008/09
recent_window <- c(2015, 2019)    # seasons 2014/15 to 2023/24

coords  <- readRDS(file.path(dir_02, "data","grid_coords.rds"))
bio     <- readRDS(file.path(dir_02, "data","biomass_matrix.rds"))
scm     <- readRDS(file.path(dir_02, "data","scores_matrices.rds"))
pca     <- readRDS(file.path(dir_02, "data","fpca_phenology.rds"))
cells   <- bio$cells
seasons <- bio$seasons

save_out <- function(x, name) saveRDS(x, file.path(out_dir, paste0(name, ".rds")))

# ---- 1. Ecoregion of each pixel ------------------------------------------------

# Each 4 km pixel is assigned to the ecoregion containing its center
eco <- st_read(eco_file, quiet = TRUE) |> st_transform(4326)
pts <- vect(as.matrix(coords[cells, c("lon", "lat")]), crs = "EPSG:4326")
ex  <- terra::extract(vect(eco), pts)
ex  <- ex[!duplicated(ex$id.y) & !is.na(ex[[eco_field]]), ]

cell_eco <- data.frame(cell = cells, coords[cells, ], ecoregion = NA_integer_)
cell_eco$ecoregion[ex$id.y] <- as.integer(ex[[eco_field]])
save_out(cell_eco, "cell_ecoregion")

trends <- lapply(c(biomass = "biomass", fPC1 = "fPC1", fPC2 = "fPC2"), function(v) {
  readRDS(file.path(dir_02, "data",sprintf("trends_%s.rds", v))) |>
    left_join(cell_eco[, c("cell", "ecoregion")], by = "cell")
})
save_out(trends, "trends_with_ecoregion")


# ---- 2. Summary statistics (Table S1 and values quoted in the text) -------------

summarise_trends <- function(df, name) {
  df |>
    filter(!is.na(slope)) |>
    mutate(sig = signif %in% TRUE) |>
    summarise(variable     = name,
              n_total      = n(),
              n_signif     = sum(sig),
              median_slope = median(slope[sig]),
              mean_slope   = mean(slope[sig]),
              perc_signif  = 100 * n_signif / n_total,
              perc_pos     = 100 * sum(sig & slope > 0) / n_signif,
              perc_neg     = 100 * sum(sig & slope < 0) / n_signif,
              .groups = "drop")
}

table_S1 <- bind_rows(lapply(names(trends), function(v)
  trends[[v]] |> filter(!is.na(ecoregion)) |> group_by(ecoregion) |> summarise_trends(v)))

summary_basin <- bind_rows(lapply(names(trends), function(v) bind_rows(
  summarise_trends(trends[[v]], v) |> mutate(domain = "all pixels"),
  summarise_trends(filter(trends[[v]], !is.na(ecoregion)), v) |> mutate(domain = "within ecoregions"))))

write.csv(table_S1,      file.path(out_dir, "table_S1_by_ecoregion.csv"), row.names = FALSE)
write.csv(summary_basin, file.path(out_dir, "summary_basin.csv"),        row.names = FALSE)


# ---- 3. Figure 4: mean trajectory of the pixels with a significant trend -------

sig_bio <- trends$biomass |> filter(signif %in% TRUE, !is.na(ecoregion))

fig4 <- bind_rows(lapply(split(match(sig_bio$cell, cells), sig_bio$ecoregion), function(r) {
  M <- bio$biomass[r, , drop = FALSE]
  data.frame(season   = seasons,
             mean_bio = colMeans(M, na.rm = TRUE),
             q25      = apply(M, 2, quantile, 0.25, na.rm = TRUE),
             q75      = apply(M, 2, quantile, 0.75, na.rm = TRUE),
             n_pixels = colSums(!is.na(M)))
}), .id = "ecoregion") |>
  mutate(ecoregion = as.integer(ecoregion))

# Sen's slope of the ecoregion averages (descriptive only: the pixels were
# selected for having a significant trend, so no p-value is computed here)
sen4 <- fig4 |>
  filter(is.finite(mean_bio)) |>
  group_by(ecoregion) |>
  group_modify(~ {
    fit <- mblm::mblm(mean_bio ~ season, dataframe = as.data.frame(.x), repeated = FALSE)
    data.frame(intercept = coef(fit)[1], slope = coef(fit)[2])
  }) |>
  ungroup()

save_out(fig4, "fig4_biomass_by_ecoregion")
save_out(sen4, "fig4_sen_by_ecoregion")


# ---- 4. Relative change in integrated Chl-a (Figure S8) ------------------------

n_years  <- max(seasons) - min(seasons)                  # 25 years
n_valid  <- rowSums(!is.na(bio$biomass))
mean_bio <- ifelse(n_valid >= min_seasons, rowMeans(bio$biomass, na.rm = TRUE), NA_real_)

rel_change <- sig_bio |>
  mutate(mean_biomass = mean_bio[match(cell, cells)],
         perc_change  = 100 * slope * n_years / mean_biomass) |>
  filter(is.finite(perc_change))
save_out(rel_change, "relative_change")


# ---- 5. Climatologies (Figures 5c-d, S7 and S10) --------------------------------

# Mean productive-season Chl-a = mean of the fitted seasonal curve
# (integral divided by the length of the season, n_steps - 1)
clim <- data.frame(cell = cells, coords[cells, ],
                   n_seasons = n_valid,
                   chl_mean  = mean_bio / (n_steps - 1))

n_valid_pc <- rowSums(!is.na(scm$scores$PC1))
clim$n_seasons_pc <- n_valid_pc
for (m in names(scm$scores)) {
  v <- rowMeans(scm$scores[[m]], na.rm = TRUE)
  v[n_valid_pc < min_seasons] <- NA_real_
  clim[[paste0("mean_", m)]] <- v
}
save_out(clim, "climatologies")


# ---- 6. Two-period comparison of the seasonal cycle -----------------------------

# For each ecoregion, sign of trend and season: mean normalized curve of the
# pixels with a significant fPC1 (or fPC2) trend
B <- eval.basis(seq_len(n_steps), pca$basis)

sig_pc <- lapply(trends[c("fPC1", "fPC2")], function(tr)
  tr |>
    filter(signif %in% TRUE, slope != 0, !is.na(ecoregion)) |>
    transmute(cell, ecoregion, trend = ifelse(slope > 0, "positive", "negative")))

curves <- list()
for (sy in seasons) {
  lab <- sprintf("%d-%d", sy, sy + 1)
  X <- readRDS(file.path(dir_02, "data","chunks", sprintf("coefs_std_%s.rds", lab)))
  m <- readRDS(file.path(dir_02, "data","chunks", sprintf("meta_std_%s.rds", lab)))
  for (v in names(sig_pc)) {
    si <- sig_pc[[v]]
    k  <- match(m$cell, si$cell)
    ok <- which(!is.na(k))
    if (length(ok) == 0) next
    idx <- split(ok, paste(si$ecoregion[k[ok]], si$trend[k[ok]], sep = "|"))
    mc  <- vapply(idx, function(j) rowMeans(X[, j, drop = FALSE]), numeric(nrow(X)))
    cv  <- B %*% mc
    key <- do.call(rbind, strsplit(colnames(cv), "|", fixed = TRUE))
    curves[[length(curves) + 1]] <- data.frame(
      mode      = v,
      ecoregion = as.integer(rep(key[, 1], each = n_steps)),
      trend     = rep(key[, 2], each = n_steps),
      season    = sy,
      n_pixels  = rep(lengths(idx), each = n_steps),
      time      = rep(seq_len(n_steps), ncol(cv)),
      value     = as.vector(cv))
  }
}
curves <- bind_rows(curves)

smooth_col <- function(x, y, span = 0.5) {
  if (sum(is.finite(y)) < 4) return(y)
  as.numeric(predict(loess(y ~ x, span = span), newdata = data.frame(x = x)))
}

# Median and interquartile range across the seasons of each period, LOESS-smoothed
period_summary <- curves |>
  mutate(period = case_when(
    season >= early_window[1]  & season <= early_window[2]  ~ "early",
    season >= recent_window[1] & season <= recent_window[2] ~ "recent")) |>
  filter(!is.na(period)) |>
  group_by(mode, ecoregion, trend, period, time) |>
  summarise(med = median(value),
            lo  = quantile(value, 0.25),
            hi  = quantile(value, 0.75), .groups = "drop") |>
  group_by(mode, ecoregion, trend, period) |>
  arrange(time, .by_group = TRUE) |>
  mutate(across(c(med, lo, hi), ~ smooth_col(time, .x))) |>
  ungroup()

peaks <- period_summary |>
  group_by(mode, ecoregion, trend, period) |>
  slice_max(med, n = 1, with_ties = FALSE) |>
  ungroup() |>
  dplyr::select(mode, ecoregion, trend, period, peak_time = time, peak_med = med)

# Share of all pixels with a significant trend in the basin (per mode)
# that fall in each ecoregion with each sign
shares <- bind_rows(lapply(names(sig_pc), function(v)
  sig_pc[[v]] |>
    count(ecoregion, trend, name = "n") |>
    mutate(mode = v, share = 100 * n / sum(n))))

mean_curve <- data.frame(time  = seq_len(n_steps),
                         value = smooth_col(seq_len(n_steps), as.vector(B %*% pca$mean_coef)))

save_out(list(curves = curves, period_summary = period_summary, peaks = peaks,
              shares = shares, mean_curve = mean_curve,
              early_window = early_window, recent_window = recent_window),
         "two_period_cycles")


