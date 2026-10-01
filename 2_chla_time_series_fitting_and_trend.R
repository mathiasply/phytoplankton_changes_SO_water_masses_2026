# =============================================================================
# 02_chl_fda_trends.R
#
# Seasonal Chl-a curves: extraction of the productive seasons, B-spline
# fitting, seasonally integrated biomass, phenological fPCA and multidecadal
# trends (Sen's slope + Mann-Kendall test).
#
# Pouly et al. (2026). Multidecadal trends in phytoplankton biomass and
# phenology in Southern Ocean water masses. JGR: Oceans.
#
# Input  : ESA OC-CCI v6.0 chlorophyll-a, 8-day composites, 4 km
#          (doi:10.5285/5011d22aae5a4671b0cbc7d05c56c4f0): one NetCDF file per
#          composite, variable `chlor_a`, start date YYYYMMDD in the file name.
#
# Outputs: outputs/02_chl_fda_trends/
#          - data/grid_coords.rds          lon/lat of every 4 km cell (cell index)
#          - data/n_valid_<season>.rds     valid composites per cell and season (Figure S1)
#          - data/chunks/                  per-season curves, coefficients and scores
#          - data/biomass_matrix.rds       integrated Chl-a, cells x seasons
#          - data/scores_matrices.rds      fPC1-fPC4 scores, cells x seasons
#          - data/fpca_phenology.rds       mean curve, harmonics, variance explained
#          - data/trends_<variable>.rds    Sen's slope and Mann-Kendall p-value per cell
#          - figures/FigS2_bspline_fits.png
# =============================================================================


# ---- 0. Packages and parameters ---------------------------------------------

library(dplyr)
library(ggplot2)
library(ncdf4)
library(fda)
library(pbapply)
library(trend)
library(here)

# Paths (relative to the repository root)
chl_dir    <- " "
out_dir    <- " "
dir_data   <- file.path(out_dir, "data")
dir_chunks <- file.path(dir_data, "chunks")
dir_fig    <- file.path(out_dir, "figures")
for (d in c(dir_chunks, dir_fig)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Study domain and period
seasons   <- 1998:2023        # season y = October of year y to March of year y + 1
lon_range <- c(-180, 180)
lat_range <- c(-83, -30)

chl_max   <- 1000             # values above are treated as missing

# Seasonal curves
n_steps    <- 23              # 8-day composites per productive season
min_steps  <- 10              # minimum valid composites per season to keep a cell
anchor_run <- 4               # consecutive missing composites at a season edge triggering an anchor
anchor_win <- 4               # composites averaged for the edge climatology of each cell

# B-spline fit and fPCA
nbasis <- 20                  # cubic B-splines
lambda <- 0.01                # roughness penalty on the second derivative
nharm  <- 10                  # harmonics computed
n_keep <- 4                   # fPC scores kept for the outputs (fPC1-fPC4)

# Trends
min_seasons <- 10             # minimum seasons with data to estimate a trend
p_thresh    <- 0.05


font_reg  <- " "
font_bold <- " "
if (file.exists(font_reg)) {
  sysfonts::font_add("LM Roman 10", regular = font_reg,
                     bold = if (file.exists(font_bold)) font_bold else NULL)
  showtext::showtext_auto()
  showtext::showtext_opts(dpi = 300)
  base_family <- "LM Roman 10"
} else {
  base_family <- "serif"
}


season_label <- function(sy) sprintf("%d-%d", sy, sy + 1)

season_dates <- function(sy, origin = as.Date("1998-09-04")) { #this date is the first date of the ESA OC-CCI v6.0 8-day composites
  start_season <- as.Date(sprintf("%d-10-01", sy))
  end_season   <- as.Date(sprintf("%d-03-31", sy + 1))

  all_dates <- seq(
    from = origin,
    to = end_season,
    by = "8 days"
  )

  out <- all_dates[
    all_dates >= start_season &
      all_dates <= end_season
  ]
  # stopifnot(length(out) == n_steps)
  out
}

# First time step of each month, for axis labels
month_breaks <- c(1, 4, 8, 12, 16, 20)
month_labels <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar")


# ---- 1. File index and grid ---------------------------------------------------

files <- list.files(chl_dir, pattern = "\\.nc$", full.names = TRUE)
dates <- as.Date(sub(".*?(\\d{8}).*", "\\1", basename(files)), format = "%Y%m%d")
files <- files[!is.na(dates)]
dates <- dates[!is.na(dates)]

nc  <- nc_open(files[1])
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nc_close(nc)

i_lon <- which(lon >= lon_range[1] & lon <= lon_range[2])
i_lat <- which(lat >= lat_range[1] & lat <= lat_range[2])
o_lon <- order(lon[i_lon])
o_lat <- order(lat[i_lat])

# Cell order used everywhere: longitude varies fastest, latitude increasing
coords <- expand.grid(lon = lon[i_lon][o_lon], lat = lat[i_lat][o_lat])
n_cell <- nrow(coords)
saveRDS(coords, file.path(dir_data, "grid_coords.rds"))

read_chl <- function(f) {
  nc <- nc_open(f)
  on.exit(nc_close(nc))
  x <- ncvar_get(nc, "chlor_a",
                 start = c(min(i_lon), min(i_lat), 1),
                 count = c(length(i_lon), length(i_lat), 1))
  x[!is.finite(x) | x > chl_max] <- NA_real_
  as.vector(x[o_lon, o_lat])
}


# ---- 2. Seasonal time series per cell ------------------------------------------

chunk_file <- function(type, sy) file.path(dir_chunks, sprintf("%s_%s.rds", type, season_label(sy)))

for (sy in seasons) {
  expected <- season_dates(sy)
  sel <- which(dates %in% expected)
  if (length(sel) < n_steps) {
    message(season_label(sy), ": ", n_steps - length(sel), " composite(s) missing")
  }
  
  # Each composite is placed at its own time step, so that a missing file
  # leaves a gap instead of shifting the following composites
  V <- matrix(NA_real_, n_steps, n_cell)
  for (k in sel) V[match(dates[k], expected), ] <- read_chl(files[k])
  
  n_ok <- colSums(is.finite(V))
  saveRDS(n_ok, file.path(dir_data, sprintf("n_valid_%s.rds", season_label(sy))))
  
  keep <- which(n_ok >= min_steps)
  saveRDS(V[, keep, drop = FALSE], chunk_file("Y_raw", sy))
  saveRDS(data.frame(cell = keep, season = sy, n_valid = n_ok[keep]),
          chunk_file("meta", sy))
  
  message(season_label(sy), ": ", length(keep), " cells kept")
  rm(V, n_ok, keep); gc(FALSE)
}


# ---- 3. Anchoring of the season edges ------------------------------------------

# When the first (last) `anchor_run` composites of a season are all missing,
# the first (last) time step is set to the cell's climatological mean of the
# first (last) `anchor_win` composites, to constrain the fit at the edges.

i_head <- seq_len(anchor_win)
i_tail <- (n_steps - anchor_win + 1):n_steps

sum_h <- sum_t <- numeric(n_cell)
n_h   <- n_t   <- integer(n_cell)

for (sy in seasons) {
  Y    <- readRDS(chunk_file("Y_raw", sy))
  cell <- readRDS(chunk_file("meta", sy))$cell
  sum_h[cell] <- sum_h[cell] + colSums(Y[i_head, , drop = FALSE], na.rm = TRUE)
  n_h[cell]   <- n_h[cell]   + colSums(is.finite(Y[i_head, , drop = FALSE]))
  sum_t[cell] <- sum_t[cell] + colSums(Y[i_tail, , drop = FALSE], na.rm = TRUE)
  n_t[cell]   <- n_t[cell]   + colSums(is.finite(Y[i_tail, , drop = FALSE]))
}
clim_h <- ifelse(n_h > 0, sum_h / n_h, NA_real_)
clim_t <- ifelse(n_t > 0, sum_t / n_t, NA_real_)
rm(sum_h, sum_t, n_h, n_t)

j_head <- seq_len(anchor_run)
j_tail <- (n_steps - anchor_run + 1):n_steps

for (sy in seasons) {
  Y    <- readRDS(chunk_file("Y_raw", sy))
  meta <- readRDS(chunk_file("meta", sy))
  
  head_na <- colSums(is.finite(Y[j_head, , drop = FALSE])) == 0
  tail_na <- colSums(is.finite(Y[j_tail, , drop = FALSE])) == 0
  jh <- which(head_na & is.finite(clim_h[meta$cell]))
  jt <- which(tail_na & is.finite(clim_t[meta$cell]))
  Y[1, jh]       <- clim_h[meta$cell[jh]]
  Y[n_steps, jt] <- clim_t[meta$cell[jt]]
  
  meta$anchor_head <- seq_along(meta$cell) %in% jh
  meta$anchor_tail <- seq_along(meta$cell) %in% jt
  saveRDS(Y,    chunk_file("Y_anchored", sy))
  saveRDS(meta, chunk_file("meta", sy))
}


# ---- 4. B-spline fit, integrated biomass and normalization ---------------------

basis <- create.bspline.basis(rangeval = c(1, n_steps), nbasis = nbasis, norder = 4)
B     <- eval.basis(seq_len(n_steps), basis)
R     <- eval.penalty(basis, Lfdobj = 2)
J     <- eval.penalty(basis, Lfdobj = 0)
int_phi <- as.vector(inprod(basis, create.constant.basis(c(1, n_steps))))

# Penalized least squares, one system per pattern of missing composites
# (identical to fitting each curve separately, much faster)
fit_block <- function(Y) {
  na     <- is.na(Y)
  key    <- as.vector(crossprod(2^(seq_len(nrow(Y)) - 1), na))
  groups <- split(seq_len(ncol(Y)), key)
  out    <- matrix(NA_real_, nbasis, ncol(Y))
  for (g in groups) {
    ok <- !na[, g[1]]
    Bk <- B[ok, , drop = FALSE]
    out[, g] <- qr.solve(crossprod(Bk) + lambda * R,
                         crossprod(Bk, Y[ok, g, drop = FALSE]))
  }
  out
}

# Cells present in at least one season (rows of the cells x seasons matrices)
all_cells <- sort(unique(unlist(lapply(seasons, function(sy)
  readRDS(chunk_file("meta", sy))$cell))))
cell_idx  <- function(cell) match(cell, all_cells)

biomass  <- matrix(NA_real_, length(all_cells), length(seasons))
sum_coef <- numeric(nbasis)
n_prof   <- 0


for (s in seq_along(seasons)) {
  sy   <- seasons[s]
  Y    <- readRDS(chunk_file("Y_anchored", sy))
  meta <- readRDS(chunk_file("meta", sy))
  
  coefs <- fit_block(Y)
  coefs[coefs < 0] <- 0      
  
  biomass[cell_idx(meta$cell), s] <- as.vector(crossprod(int_phi, coefs))
  
  mu <- colMeans(coefs)
  sd <- sqrt(colSums(sweep(coefs, 2, mu)^2) / (nbasis - 1))
  ok <- is.finite(sd) & sd > 0
  coefs_std <- sweep(sweep(coefs[, ok, drop = FALSE], 2, mu[ok]), 2, sd[ok], "/")
  
  saveRDS(coefs_std,     chunk_file("coefs_std", sy))
  saveRDS(meta[ok, ],    chunk_file("meta_std", sy))
  
  sum_coef <- sum_coef + rowSums(coefs_std)
  n_prof   <- n_prof + ncol(coefs_std)
  message(season_label(sy), ": ", ncol(coefs_std), " curves fitted")
  rm(Y, coefs, coefs_std); gc(FALSE)
}

saveRDS(list(cells = all_cells, seasons = seasons, biomass = biomass),
        file.path(dir_data, "biomass_matrix.rds"))



# ---- 5. Phenological fPCA of the normalized curves ------------------------------

mean_coef <- sum_coef / n_prof

C <- matrix(0, nbasis, nbasis)
for (sy in seasons) {
  X <- readRDS(chunk_file("coefs_std", sy)) - mean_coef
  C <- C + tcrossprod(X)
}
C <- C / (n_prof - 1)

eig       <- eigen(C %*% J)
values    <- Re(eig$values)
harmonics <- Re(eig$vectors[, seq_len(nharm)])
varprop   <- values[seq_len(nharm)] / sum(values[values > 0])

# Sign convention of the paper: positive fPC1 = early (October) bloom,
# positive fPC2 = bloom apex in November-December
h_eval <- B %*% harmonics
if (h_eval[1, 1] < 0) harmonics[, 1] <- -harmonics[, 1]
if (h_eval[9, 2] < 0) harmonics[, 2] <- -harmonics[, 2]

pca <- list(basis = basis, mean_coef = mean_coef, harmonics = harmonics,
            values = values, varprop = varprop, n_curves = n_prof)
saveRDS(pca, file.path(dir_data, "fpca_phenology.rds"))


# ---- 6. Scores ------------------------------------------------------------------

scores <- replicate(n_keep, matrix(NA_real_, length(all_cells), length(seasons)),
                    simplify = FALSE)
names(scores) <- paste0("PC", seq_len(n_keep))

for (s in seq_along(seasons)) {
  sy   <- seasons[s]
  X    <- readRDS(chunk_file("coefs_std", sy)) - mean_coef
  meta <- readRDS(chunk_file("meta_std", sy))
  sc   <- crossprod(X, J %*% harmonics[, seq_len(n_keep)])
  for (m in seq_len(n_keep)) scores[[m]][cell_idx(meta$cell), s] <- sc[, m]
}

saveRDS(list(cells = all_cells, seasons = seasons, scores = scores),
        file.path(dir_data, "scores_matrices.rds"))


# ---- 7. Multidecadal trends: Sen's slope and Mann-Kendall test -----------------

CB <- lapply(seq_along(seasons), function(k) if (k >= 2) utils::combn(k, 2) else NULL)

calc_sen_coef <- function(y) {
  ok <- !is.na(y)
  n  <- sum(ok)
  if (n < min_seasons) return(c(NA_real_, n, NA_real_))
  yv <- y[ok]
  xv <- seasons[ok]
  if (var(yv) == 0) return(c(0, n, NA_real_))
  ij    <- CB[[n]]
  slope <- median((yv[ij[2, ]] - yv[ij[1, ]]) / (xv[ij[2, ]] - xv[ij[1, ]]))
  pval  <- tryCatch(trend::mk.test(yv)$p.value, error = function(e) NA_real_)
  c(slope, n, pval)
}

# M: cells x seasons matrix (rows = all_cells, columns = seasons)
save_trends <- function(M, name, cl = NULL) {
  res <- pbapply::pbapply(M, 1, calc_sen_coef, cl = cl)      # 3 x cells
  tr  <- data.frame(cell = all_cells, coords[all_cells, ],
                    slope = res[1, ], n_valid = res[2, ], pval = res[3, ])
  tr$signif <- tr$pval < p_thresh
  saveRDS(tr, file.path(dir_data, sprintf("trends_%s.rds", name)))
  message(name, ": ", round(100 * mean(tr$signif, na.rm = TRUE), 1),
          "% of cells with a significant trend")
}

save_trends(biomass,    "biomass")
save_trends(scores$PC1, "fPC1")
save_trends(scores$PC2, "fPC2")


# ---- 8. Figure S2: examples of B-spline fits -------------------------------------

# One random curve for each number of valid composites (10 to 23),
# in the 2012/13 season
set.seed(1)
sy_ex <- 2012
Y_raw <- readRDS(chunk_file("Y_raw", sy_ex))
Y_anc <- readRDS(chunk_file("Y_anchored", sy_ex))
n_val <- colSums(is.finite(Y_raw))

pick <- sapply(min_steps:n_steps, function(n) {
  idx <- which(n_val == n)
  if (length(idx) == 0) NA_integer_ else idx[sample.int(length(idx), 1)]
})
pick <- pick[!is.na(pick)]

coef_ex <- fit_block(Y_anc[, pick, drop = FALSE])
coef_ex[coef_ex < 0] <- 0
t_fine  <- seq(1, n_steps, length.out = 500)

df_obs <- data.frame(t = rep(seq_len(n_steps), length(pick)),
                     chl = as.vector(Y_raw[, pick]),
                     n = rep(n_val[pick], each = n_steps))
df_fit <- data.frame(t = rep(t_fine, length(pick)),
                     chl = as.vector(eval.basis(t_fine, basis) %*% coef_ex),
                     n = rep(n_val[pick], each = length(t_fine)))

fig_S2 <- ggplot() +
  geom_point(data = df_obs, aes(t, chl), colour = "grey20", alpha = 0.7, na.rm = TRUE) +
  geom_line(data = df_fit, aes(t, chl), colour = "blue", linewidth = 0.8) +
  facet_wrap(~ n, scales = "free_y", labeller = labeller(n = function(x) paste("n =", x))) +
  scale_x_continuous(breaks = month_breaks, labels = month_labels) +
  labs(x = "Month", y = expression(paste("Chl-", italic(a), " (mg ", m^{-3}, ")"))) +
  theme_bw(base_size = 16, base_family = base_family) +
  theme(panel.grid = element_blank())

ggsave(file.path(dir_fig, "FigS2_bspline_fits.png"), fig_S2,
       width = 40, height = 30, units = "cm", dpi = 300)

message("Done. Outputs in ", out_dir)
