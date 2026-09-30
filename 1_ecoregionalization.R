# =============================================================================
# 01_ecoregionalization.R
#
# Hydrological ecoregionalization of the Southern Ocean from the bivariate
# functional PCA of temperature-salinity profiles (0-1000 m).
#
# Pouly et al. (2026). Multidecadal trends in phytoplankton biomass and
# phenology in Southern Ocean water masses. JGR: Oceans.
#
# Input  : monthly temperature and salinity of the Global Ocean Ensemble Physics
#          Reanalysis (GLOBAL_MULTIYEAR_PHY_ENS_001_031, doi:10.48670/moi-00024),
#          NetCDF files with variables `thetao_mean` and `so_mean`,
#          0-1000 m, 90-30 S, October-March, seasons 1998/99 to 2023/24.
#
# Outputs: outputs/01_ecoregionalization/
#          - ecoregions/ecoregions_k8.shp     ecoregion polygons 
#          - ecoregions/cluster_sd.tif        inter-month SD of cluster number (Figure 1b)
#          - data/                            profiles, fPCA object and scores
#          - figures/                         Figures 2, S3, S4 and S5 of the paper
# =============================================================================


# ---- 0. Packages and parameters ---------------------------------------------

library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(data.table)
library(ncdf4)
library(fda)        # fd(), eval.fd()
library(fda.oce)    # remotes::install_github("EPauthenet/fda.oce")
library(cluster)    # clusGap(), maxSE()
library(terra)
library(sf)
library(here)

# Paths (relative to the repository root)
raw_dir <-  "E:/TS DATA/Raw data 1998 2023 all SO/all_zone"         # downloaded NetCDF files
out_dir <- "G:/papier_phyto_change_2026/test code/ecoregion"
dir_data <- file.path(out_dir, "data")
dir_eco  <- file.path(out_dir, "ecoregions")
dir_fig  <- file.path(out_dir, "figures")
for (d in c(dir_data, dir_eco, dir_fig)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Study domain and period
seasons     <- 1998:2023          # season y = October of year y to March of year y + 1
months_kept <- c(10, 11, 12, 1, 2, 3)
lon_range   <- c(-180, 180)
lat_range   <- c(-90, -30)
depth_range <- c(0, 1000)
var_temp    <- "thetao_mean"
var_sal     <- "so_mean"

# Clustering
n_pc       <- 3      # number of thermohaline modes used for clustering
k_retained <- 8      # number of clusters retained (see gap statistic, Figure S3)

eco_palette <- c("#006837", "#66BD63", "#D9EF8B", "#F0F097",
                 "#ACD1DD", "#CACCE0", "#8073AC", "#542788")


font_reg  <- "C:/Users/etudhenri5/Desktop/lmroman10-regular.otf"
font_bold <- "C:/Users/etudhenri5/Desktop/lmroman10-bold.otf"
if (file.exists(font_reg)) {
  sysfonts::font_add("LM Roman 10", regular = font_reg,
                     bold = if (file.exists(font_bold)) font_bold else NULL)
  showtext::showtext_auto()
  showtext::showtext_opts(dpi = 300)
  base_family <- "LM Roman 10"
} else {
  base_family <- "serif"
}

theme_pub <- theme_bw(base_size = 18, base_family = base_family) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95"))

# Season of a given month: October-December belong to season y,
# January-March to season y - 1
season_of <- function(year, month) ifelse(month >= 10, year, year - 1)

# Decode the NetCDF time axis from its "units" attribute
# (e.g. "hours since 1950-01-01 00:00:00")
decode_time <- function(nc) {
  t     <- ncvar_get(nc, "time")
  units <- ncatt_get(nc, "time", "units")$value
  parts <- regmatches(units, regexec("^(\\w+) since (.+)$", units))[[1]]
  mult  <- switch(parts[2], seconds = 1, minutes = 60, hours = 3600, days = 86400,
                  stop("Unknown time unit: ", units))
  origin <- as.POSIXct(substr(parts[3], 1, 10), tz = "UTC")
  as.POSIXct(t * mult, origin = origin, tz = "UTC")
}


# ---- 1. Read and format the temperature-salinity profiles -------------------

nc_files <- list.files(raw_dir, pattern = "\\.nc$", full.names = TRUE)
if (length(nc_files) == 0) stop("No NetCDF file found in ", raw_dir)

T_list <- S_list <- meta_list <- list()
Pi <- NULL

for (f in nc_files) {
  nc <- nc_open(f)
  lon_all   <- ncvar_get(nc, "longitude")
  lat_all   <- ncvar_get(nc, "latitude")
  depth_all <- ncvar_get(nc, "depth")
  dates     <- decode_time(nc)
  
  i_lon <- which(lon_all   >= lon_range[1]   & lon_all   <= lon_range[2])
  i_lat <- which(lat_all   >= lat_range[1]   & lat_all   <= lat_range[2])
  i_dep <- which(depth_all >= depth_range[1] & depth_all <= depth_range[2])
  if (is.null(Pi)) Pi <- as.numeric(depth_all[i_dep])
  
  yr <- as.integer(format(dates, "%Y"))
  mo <- as.integer(format(dates, "%m"))
  i_time <- which(mo %in% months_kept & season_of(yr, mo) %in% seasons)
  
  # Longitude varies fastest, as in the NetCDF arrays
  grid   <- expand.grid(lon = lon_all[i_lon], lat = lat_all[i_lat])
  n_cell <- nrow(grid)
  n_dep  <- length(i_dep)
  
  for (j in i_time) {
    start <- c(min(i_lon), min(i_lat), min(i_dep), j)
    count <- c(length(i_lon), length(i_lat), n_dep, 1)
    temp  <- ncvar_get(nc, var_temp, start = start, count = count)
    sal   <- ncvar_get(nc, var_sal,  start = start, count = count)
    
    # (lon, lat, depth) array -> depth x profile matrix
    T_m <- t(matrix(temp, nrow = n_cell, ncol = n_dep))
    S_m <- t(matrix(sal,  nrow = n_cell, ncol = n_dep))
    
    # Keep only profiles complete over 0-1000 m in both variables
    # (this excludes regions shallower than 1000 m)
    ok <- colSums(is.na(T_m)) == 0 & colSums(is.na(S_m)) == 0
    if (!any(ok)) next
    
    T_list[[length(T_list) + 1]]       <- T_m[, ok, drop = FALSE]
    S_list[[length(S_list) + 1]]       <- S_m[, ok, drop = FALSE]
    meta_list[[length(meta_list) + 1]] <- data.frame(
      grid[ok, ], year = yr[j], month = mo[j], season = season_of(yr[j], mo[j]))
    
    message(format(dates[j], "%Y-%m"), ": ", sum(ok), " profiles")
  }
  nc_close(nc)
}

Temp <- do.call(cbind, T_list)
Sal  <- do.call(cbind, S_list)
meta <- bind_rows(meta_list)
rm(T_list, S_list, meta_list); gc()

saveRDS(Temp, file.path(dir_data, "temperature_profiles.rds"))
saveRDS(Sal,  file.path(dir_data, "salinity_profiles.rds"))
saveRDS(meta, file.path(dir_data, "profiles_metadata.rds"))
saveRDS(Pi,   file.path(dir_data, "depth_levels.rds"))

message(ncol(Temp), " profiles on ", length(Pi), " depth levels")


# ---- 2. Bivariate functional PCA ---------------------------------------------

# Array levels x profiles x variables, as expected by fda.oce
Xi <- array(NA_real_, dim = c(nrow(Temp), ncol(Temp), 2))
Xi[, , 1] <- Temp
Xi[, , 2] <- Sal

# here, we order the profiles by increasing depth levels, as required by fda.oce
ord <- order(Pi, na.last = NA)
Pi <- ceiling(as.numeric(Pi))
Pi <- Pi[ord]
if (length(dim(Xi)) == 2) {
  Xi <- Xi[ord, , drop = FALSE]
} else if (length(dim(Xi)) == 3) {
  Xi <- Xi[ord, , , drop = FALSE]
}

# Note: the fda.oce functions create their outputs in the global environment
fda.oce::bspl(Pi, Xi)          # B-spline fit (20 cubic B-splines by default) -> `fdobj`
fda.oce::fpca(fdobj)           # bivariate fPCA                              -> `pca`
fda.oce::proj(fdobj, pca)      # scores of each profile                      -> `pc`
rm(Xi); gc()

saveRDS(pca, file.path(dir_data, "fpca_thermohaline.rds"))
saveRDS(pc,  file.path(dir_data, "fpca_scores.rds"))


# ---- 3. Number of clusters: gap statistic (Figure S3) ------------------------

features <- pc[, 1:n_pc]

set.seed(1)
idx_gap <- sample(nrow(features), 1e4)
gap <- clusGap(features[idx_gap, ], FUN = kmeans, K.max = 15, B = 50,
               nstart = 25, iter.max = 50)

gap_tab <- as.data.frame(gap$Tab) |> mutate(k = row_number())
k_firstSEmax <- maxSE(gap_tab$gap, gap_tab$SE.sim, method = "firstSEmax")
message("Gap statistic, firstSEmax rule: k = ", k_firstSEmax,
        " (retained: k = ", k_retained, ")")

saveRDS(gap_tab, file.path(dir_data, "gap_statistic.rds"))

fig_S3 <- ggplot(gap_tab, aes(k, gap)) +
  geom_vline(xintercept = k_retained, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = gap_tab$k) +
  labs(x = "Number of clusters (k)", y = "Gap statistic") +
  theme_pub

ggsave(file.path(dir_fig, "FigS3_gap_statistic.png"), fig_S3,
       width = 20, height = 14, units = "cm", dpi = 300)


# ---- 4. K-means clustering and numbering of the clusters ---------------------

set.seed(1)
km <- kmeans(features, centers = k_retained, nstart = 25,
             iter.max = 20, algorithm = "Lloyd")

# Number the clusters from north (1) to south (k) by median latitude
med_lat <- tapply(meta$lat, km$cluster, median)
new_id  <- rank(-med_lat)
meta$cluster <- as.integer(new_id[as.character(km$cluster)])

profiles <- cbind(meta, setNames(as.data.frame(pc[, 1:n_pc]), paste0("PC", 1:n_pc)))
saveRDS(profiles, file.path(dir_data, "profiles_scores_clusters.rds"))


# ---- 5. Ecoregion layers ------------------------------------------------------

# Mean and inter-month standard deviation of the cluster number per 0.25 deg cell
pix <- as.data.table(meta)[, .(cluster_mean = mean(cluster),
                               cluster_sd   = sd(cluster),
                               n_months     = .N),
                           by = .(lon, lat)]
saveRDS(pix, file.path(dir_data, "cluster_mean_sd_per_cell.rds"))

# SD of the cluster number: high values mark transition zones (Figure 1b)
r_sd <- rast(as.data.frame(pix[, .(lon, lat, cluster_sd)]),
             type = "xyz", crs = "EPSG:4326")
writeRaster(r_sd, file.path(dir_eco, "cluster_sd.tif"), overwrite = TRUE)

# Ecoregions: rounding the mean cluster number places each boundary on the
# median line of the transition zone between two clusters
r_eco <- rast(as.data.frame(pix[, .(lon, lat, ecoregion = round(cluster_mean))]),
              type = "xyz", crs = "EPSG:4326")
eco_poly <- as.polygons(r_eco, dissolve = TRUE) |> st_as_sf()
st_write(eco_poly, file.path(dir_eco, "ecoregions_k8.shp"), delete_dsn = TRUE)


# ---- 6. Figures ----------------------------------------------------------------

# Figure 2: mean temperature and salinity profiles per ecoregion (+/- 1 SD)
profile_stats <- function(mat, cl, var_name) {
  map_dfr(sort(unique(cl)), function(g) {
    x <- mat[, cl == g, drop = FALSE]
    n <- ncol(x)
    m <- rowMeans(x)
    s <- sqrt(pmax(rowMeans(x^2) - m^2, 0) * n / (n - 1))
    data.frame(depth = Pi, ecoregion = g, variable = var_name,
               mean = m, lo = m - s, hi = m + s)
  })
}

df_prof <- bind_rows(profile_stats(Temp, meta$cluster, "temperature"),
                     profile_stats(Sal,  meta$cluster, "salinity"))
saveRDS(df_prof, file.path(dir_data, "mean_profiles_per_ecoregion.rds"))

plot_profile <- function(var_name, x_lab) {
  ggplot(filter(df_prof, variable == var_name),
         aes(y = depth, colour = factor(ecoregion), fill = factor(ecoregion))) +
    geom_ribbon(aes(xmin = lo, xmax = hi), orientation = "y",
                alpha = 0.18, colour = NA) +
    geom_path(aes(x = mean), linewidth = 1.2) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_colour_manual(name = "Ecoregion", values = eco_palette) +
    scale_fill_manual(name = "Ecoregion", values = eco_palette) +
    labs(x = x_lab, y = "Depth (m)") +
    theme_pub
}

fig_2 <- (plot_profile("temperature", "Temperature (°C)") + labs(tag = "(a)")) +
  (plot_profile("salinity", "Practical salinity") + labs(y = NULL, tag = "(b)")) +
  plot_layout(guides = "collect")

ggsave(file.path(dir_fig, "Fig2_TS_profiles.png"), fig_2,
       width = 30, height = 22, units = "cm", dpi = 300)


# Figure S4: effect of the first three modes on the mean profiles
z <- seq(pca$basis$rangeval[1], pca$basis$rangeval[2], length.out = 200)
var_names <- if (!is.null(pca$fdnames)) {
  pca$fdnames[3:(2 + pca$ndim)]
} else {
  c("Temperature", "Salinity")
}

df_eig <- map_dfr(1:3, function(m) {
  map_dfr(seq_len(pca$ndim), function(v) {
    idx   <- ((v - 1) * pca$nbas + 1):(v * pca$nbas)
    c_mu  <- pca$Cm[idx]
    c_eig <- pca$axe[idx, m]
    data.frame(
      depth    = z,
      variable = var_names[[v]],
      mode     = paste0("fPC", m),
      mean     = as.numeric(eval.fd(z, fd(c_mu, pca$basis))),
      plus     = as.numeric(eval.fd(z, fd(c_mu + c_eig, pca$basis))),
      minus    = as.numeric(eval.fd(z, fd(c_mu - c_eig, pca$basis)))
    )
  })
})

fig_S4 <- ggplot(df_eig, aes(y = depth)) +
  geom_path(aes(x = mean),  colour = "black",   linewidth = 1) +
  geom_path(aes(x = plus),  colour = "#D33F6A", linewidth = 1) +
  geom_path(aes(x = minus), colour = "#4A6FE3", linewidth = 1) +
  facet_grid(mode ~ variable, scales = "free_x") +
  scale_y_reverse(expand = c(0, 0)) +
  labs(x = NULL, y = "Depth (m)") +
  theme_pub +
  theme(panel.spacing.x = unit(1.5, "lines"))

ggsave(file.path(dir_fig, "FigS4_eigenfunctions.png"), fig_S4,
       width = 15, height = 25, units = "cm", dpi = 300)


# Figure S5: thermohaline score space colored by ecoregion

df_scores <- profiles |>
  slice_sample(n = min(5e5, nrow(profiles))) |>
  mutate(ecoregion = factor(cluster))

pct <- round(pca$pval[1:3], 1)

plot_scores <- function(x, y, i, j) {
  ggplot(df_scores, aes(.data[[x]], .data[[y]], colour = ecoregion)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_point(size = 0.3, alpha = 0.25) +
    scale_colour_manual(name = "Ecoregion", values = eco_palette) +
    labs(x = paste0("fPC", i, " (", pct[i], "%)"),
         y = paste0("fPC", j, " (", pct[j], "%)")) +
    theme_pub
}

fig_S5 <- (plot_scores("PC1", "PC2", 1, 2) + labs(tag = "(a)")) +
  (plot_scores("PC2", "PC3", 2, 3) + labs(tag = "(b)")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4, alpha = 1)))

ggsave(file.path(dir_fig, "FigS5_score_space.png"), fig_S5,
       width = 40, height = 20, units = "cm", dpi = 300)

message("Done. Outputs in ", out_dir)