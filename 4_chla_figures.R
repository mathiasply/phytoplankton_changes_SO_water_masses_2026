# =============================================================================
# 04_figures.R
#
# Figures of the paper and of the Supporting Information built from the
# outputs of scripts 01-03 (Figures 2, S2, S3, S4 and S5 are produced by
# scripts 01 and 02).
#
# Pouly et al. (2026). Multidecadal trends in phytoplankton biomass and
# phenology in Southern Ocean water masses. JGR: Oceans.
#
# Main text : Figures 1, 3, 4, 5, 6, 7, 8
# SI        : Figures S1, S6, S7, S8, S9, S10, S11, S12, S13
# Outputs   : outputs/04_figures/
# =============================================================================


# ---- 0. Packages, paths and parameters ----------------------------------------

library(dplyr)
library(ggplot2)
library(patchwork)
library(sf)
library(terra)
library(ggpattern)
library(ggridges)
library(fda)
library(here)

dir_01  <- "G:/papier_phyto_change_2026/test code/ecoregion"
dir_02  <- "G:/papier_phyto_change_2026/test code/chla_fit"
dir_03 <- "G:/papier_phyto_change_2026/test code/chla_analysis"
fig_dir <- "G:/papier_phyto_change_2026/test code/chla_fig"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

eco_file   <- file.path(dir_01, "ecoregions", "ecoregions_k8.shp")
eco_field  <- "ecoregion"                                       # "zone" for layer_k8.shp
label_file <- here("data", "ancillary", "ecoregion_labels.shp")  # optional label positions

n_steps      <- 23
month_breaks <- c(1, 4, 8, 12, 16, 20)        # first time step of each month
month_labels <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar")

eco_palette <- c("1" = "#006837", "2" = "#66BD63", "3" = "#D9EF8B", "4" = "#F0F097",
                 "5" = "#ACD1DD", "6" = "#CACCE0", "7" = "#8073AC", "8" = "#542788")

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
        strip.background = element_rect(fill = "grey98"))

save_fig <- function(p, name, width = 8, height = 8) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300, bg = "white")
}

colourbar <- function() {
  guides(fill = guide_colourbar(barwidth = unit(12, "cm"), barheight = unit(0.4, "cm"),
                                ticks.colour = "black", frame.colour = "black"))
}

# Inputs
coords <- readRDS(file.path(dir_02, "data","grid_coords.rds"))
pca    <- readRDS(file.path(dir_02, "data","fpca_phenology.rds"))
trends <- readRDS(file.path(dir_03, "trends_with_ecoregion.rds"))
clim   <- readRDS(file.path(dir_03, "climatologies.rds"))
tp     <- readRDS(file.path(dir_03, "two_period_cycles.rds"))


# ---- 1. Polar map layers -------------------------------------------------------

polar_crs <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +datum=WGS84 +units=m +no_defs"

p30      <- st_coordinates(st_transform(st_sfc(st_point(c(0, -30)), crs = 4326), polar_crs))
radius_m <- sqrt(sum(p30^2))                                  # map edge at 30 S
circle   <- st_buffer(st_sfc(st_point(c(0, 0)), crs = polar_crs), radius_m)
box      <- st_as_sfc(st_bbox(c(xmin = -1.05, xmax = 1.05, ymin = -1.05, ymax = 1.05) * radius_m,
                              crs = st_crs(polar_crs)))
mask     <- st_difference(box, circle)

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
  st_union() |> st_transform(polar_crs) |> st_make_valid() |> st_intersection(circle)

ticks <- do.call(rbind, lapply(seq(0, 315, by = 45), function(lon) {
  p <- st_coordinates(st_transform(st_sfc(st_point(c(lon, -45)), crs = 4326), polar_crs))
  u <- p / sqrt(sum(p^2))
  data.frame(x = u[1] * radius_m, y = u[2] * radius_m,
             xend = u[1] * radius_m * 1.03, yend = u[2] * radius_m * 1.03)
}))

# Ecoregions, projected and lightly smoothed for display (modal filter, 10 km grid)
eco_raw  <- st_read(eco_file, quiet = TRUE) |> st_transform(polar_crs)
r_eco    <- rasterize(vect(eco_raw), rast(vect(eco_raw), resolution = 10000), field = eco_field)
r_eco    <- focal(r_eco, w = matrix(1, 7, 7), fun = "modal", na.rm = TRUE)
eco_poly <- as.polygons(r_eco) |> st_as_sf() |> st_make_valid()
names(eco_poly)[1] <- "ecoregion"

# Regions outside the ecoregions (shallower than 1000 m), hatched on the maps
outside <- st_difference(circle, st_union(eco_poly)) |> st_make_valid()

# Ecoregion numbers: positions from a file if available, else inside each ecoregion
if (file.exists(label_file)) {
  lab_pts <- st_read(label_file, quiet = TRUE) |> st_transform(polar_crs)
  labels_df <- data.frame(st_coordinates(lab_pts), label = seq_len(nrow(lab_pts)))
} else {
  lab_pts <- eco_poly |>
    st_cast("POLYGON") |>
    mutate(area = as.numeric(st_area(geometry))) |>
    group_by(ecoregion) |>
    slice_max(area, n = 1) |>
    ungroup() |>
    st_point_on_surface()
  labels_df <- data.frame(st_coordinates(lab_pts), label = lab_pts$ecoregion)
}

# Regular lon/lat template of the 4 km grid (cell index of script 02)
lon_u <- sort(unique(coords$lon))
lat_u <- sort(unique(coords$lat))
dx <- median(diff(lon_u)); dy <- median(diff(lat_u))
tmpl <- rast(xmin = min(lon_u) - dx / 2, xmax = max(lon_u) + dx / 2,
             ymin = min(lat_u) - dy / 2, ymax = max(lat_u) + dy / 2,
             ncols = length(lon_u), nrows = length(lat_u), crs = "EPSG:4326")
raster_cell <- function(cell) {
  i <- (cell - 1) %% length(lon_u) + 1
  j <- (cell - 1) %/% length(lon_u) + 1
  (length(lat_u) - j) * length(lon_u) + i
}

# Values on grid cells -> polar stereographic data frame (x, y, value)
to_polar <- function(cell, value, method = "bilinear") {
  v <- rep(NA_real_, ncell(tmpl))
  v[raster_cell(cell)] <- value
  r <- project(setValues(tmpl, v), polar_crs, method = method, res = 4000)
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "value"
  d
}

# Pixels with fewer than 10 valid seasons around Antarctica
r_na <- setValues(tmpl, NA_real_)
r_na[raster_cell(clim$cell[clim$n_seasons_pc >= 10])] <- 1
r_na <- ifel(is.na(r_na) & init(r_na, "y") < -55, 1, NA)
poly_na <- as.polygons(aggregate(r_na, 3, fun = "max", na.rm = TRUE)) |>
  st_as_sf() |> st_transform(polar_crs) |> st_make_valid() |> st_difference(land)
r_tmp   <- rasterize(vect(poly_na), rast(vect(poly_na), resolution = 10000), field = 1)
poly_na <- as.polygons(focal(r_tmp, w = matrix(1, 3, 3), fun = "max", na.policy = "only")) |>
  st_as_sf() |> st_make_valid()

# Layers common to all maps
map_layers <- function(na_fill = "grey85", show_na = TRUE, boundaries = TRUE,
                       labels = TRUE, extra = NULL) {
  list(
    if (show_na) geom_sf(data = poly_na, fill = na_fill, colour = NA, inherit.aes = FALSE),
    geom_sf_pattern(data = outside, pattern = "stripe", fill = "white", colour = NA,
                    pattern_colour = "grey55", pattern_fill = "grey55", pattern_angle = 45,
                    pattern_density = 0.2, pattern_size = 0.2, pattern_spacing = 0.005,
                    alpha = 0.1, inherit.aes = FALSE),
    geom_sf(data = outside, fill = NA, colour = "grey55", linewidth = 0.4, inherit.aes = FALSE),
    if (boundaries) geom_sf(data = eco_poly, fill = NA, colour = "black",
                            linewidth = 0.4, inherit.aes = FALSE),
    extra,
    if (labels) geom_text(data = labels_df, aes(X, Y, label = label), size = 6,
                          fontface = "bold", family = base_family, inherit.aes = FALSE),
    geom_sf(data = land, fill = "grey40", colour = "grey40", linewidth = 0.2, inherit.aes = FALSE),
    geom_sf(data = mask, fill = "white", colour = NA, inherit.aes = FALSE),
    geom_sf(data = circle, fill = NA, colour = "black", linewidth = 0.6, inherit.aes = FALSE),
    geom_segment(data = ticks, aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.6, inherit.aes = FALSE),
    coord_sf(crs = polar_crs, xlim = c(-1.03, 1.03) * radius_m,
             ylim = c(-1.03, 1.03) * radius_m, expand = FALSE, clip = "off"),
    theme_void(base_family = base_family),
    theme(legend.position = "bottom", legend.title.position = "top",
          legend.title = element_text(size = 18, face = "bold", hjust = 0.5),
          legend.text  = element_text(size = 13))
  )
}

raster_map <- function(d, fill_scale, fill_lab, ...) {
  ggplot() +
    geom_raster(data = d, aes(x, y, fill = value)) +
    fill_scale +
    map_layers(...) +
    labs(fill = fill_lab)
}

# Symmetric, pseudo-log stepped scale for slopes
slope_scale <- function(values, n = 10, sigma = 0.01) {
  q    <- max(abs(quantile(values, c(0.01, 0.99), na.rm = TRUE)))
  plog <- scales::transform_pseudo_log(sigma = sigma)
  br   <- unique(round(plog$inverse(seq(plog$transform(-q), plog$transform(q),
                                        length.out = n)), 2))
  scale_fill_stepsn(colours = rev(hcl.colors(n, "Spectral")), limits = c(-q, q),
                    breaks = br, transform = plog, oob = scales::squish,
                    na.value = "transparent")
}

# Diverging stepped scale for mean fPC scores
score_scale <- function(values, n = 8) {
  q <- max(abs(quantile(values, c(0.01, 0.99), na.rm = TRUE)))
  scale_fill_stepsn(colours = hcl.colors(n - 1, "Blue-Red 2"),
                    breaks = round(seq(-q, q, length.out = n), 2), limits = c(-q, q),
                    oob = scales::squish, na.value = "transparent")
}

# Significance classes of the Mann-Kendall test
sig_levels <- c("Non-significant", "p < 0.05", "p < 0.01", "p < 0.001")
sig_cols   <- setNames(c("white", "#9FE1CB", "#1D9E75", "#0F6E56"), sig_levels)

sig_map <- function(tr) {
  cls <- c(3, 2, 1, 0)[cut(tr$pval, c(-Inf, 0.001, 0.01, 0.05, Inf), labels = FALSE)]
  d   <- to_polar(tr$cell, cls, method = "near")
  d$value <- factor(d$value, levels = 0:3, labels = sig_levels)
  ggplot() +
    geom_raster(data = d, aes(x, y, fill = value)) +
    scale_fill_manual(values = sig_cols, drop = FALSE, na.value = "transparent") +
    map_layers() +
    labs(fill = "Significance") +
    guides(fill = guide_legend(nrow = 1, label.position = "bottom",
                               keywidth = unit(1.5, "cm"), keyheight = unit(0.4, "cm"))) +
    theme(legend.key = element_rect(colour = "black"))
}

# Density of significant slopes per ecoregion (Figures 3b, 6c-d).
# coord_cartesian() only zooms: slopes beyond x_lim still enter densities and medians
ridge_plot <- function(tr, x_lab, x_lim) {
  d <- tr |>
    filter(signif %in% TRUE, !is.na(ecoregion)) |>
    mutate(eco = factor(ecoregion, levels = 8:1))
  ggplot(d, aes(slope, eco, fill = eco)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_density_ridges(scale = 2, rel_min_height = 0.009, quantile_lines = TRUE,
                        quantiles = 2, colour = "black", linewidth = 0.35, alpha = 0.85) +
    scale_fill_manual(values = rev(eco_palette), guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(mult = c(0.02, 0.18))) +
    coord_cartesian(xlim = x_lim) +
    labs(x = x_lab, y = "Ecoregion") +
    theme_pub
}


# ---- 2. Figure 1: ecoregions and inter-month SD of the cluster number ------------

orsi <- st_as_sf(orsifronts::orsifronts) |> st_transform(polar_crs)

fig1a <- ggplot() +
  geom_sf(data = eco_poly, aes(fill = factor(ecoregion)), colour = NA) +
  scale_fill_manual(values = eco_palette, name = "Ecoregion") +
  map_layers(show_na = FALSE, boundaries = FALSE, labels = FALSE,
             extra = geom_sf(data = orsi, colour = "red", linewidth = 0.4,
                             inherit.aes = FALSE)) +
  guides(fill = guide_legend(nrow = 1))

sd_proj <- as.data.frame(project(rast(file.path(dir_01, "ecoregions", "cluster_sd.tif")),
                                 polar_crs, method = "bilinear", res = 10000),
                         xy = TRUE, na.rm = TRUE)
names(sd_proj)[3] <- "value"

fig1b <- raster_map(sd_proj,
                    scale_fill_stepsn(colours = viridisLite::turbo(9), n.breaks = 9,
                                      na.value = "transparent", oob = scales::squish),
                    "Cluster standard deviation",
                    show_na = FALSE, boundaries = FALSE, labels = FALSE) +
  colourbar()

save_fig(fig1a + fig1b + plot_annotation(tag_levels = list(c("(a)", "(b)"))),
         "Fig1_ecoregions", width = 16, height = 9)


# ---- 3. Figure 3: biomass trends ---------------------------------------------------

sig_bio <- filter(trends$biomass, signif %in% TRUE)

fig3a <- raster_map(to_polar(sig_bio$cell, sig_bio$slope), slope_scale(sig_bio$slope),
                    expression("Biomass slope (mg m"^-3~"yr"^-1*")")) +
  colourbar()
fig3b <- ridge_plot(trends$biomass, expression("Biomass slope (mg m"^-3~"yr"^-1*")"),
                    c(-0.5, 0.5))

save_fig(fig3a + fig3b + plot_layout(widths = c(1.2, 1)) +
           plot_annotation(tag_levels = list(c("(a)", "(b)"))),
         "Fig3_biomass_trends", width = 16, height = 8)


# ---- 4. Figure 4: mean trajectory per ecoregion ------------------------------------

f4   <- readRDS(file.path(dir_03, "fig4_biomass_by_ecoregion.rds"))
sen4 <- readRDS(file.path(dir_03, "fig4_sen_by_ecoregion.rds"))
dark <- setNames(colorspace::darken(eco_palette, amount = 0.15), names(eco_palette))

fig4 <- ggplot(f4, aes(season, mean_bio, colour = factor(ecoregion), fill = factor(ecoregion))) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, colour = NA) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_abline(data = sen4, aes(intercept = intercept, slope = slope,
                               colour = factor(ecoregion)), linewidth = 1.1) +
  geom_text(data = sen4, aes(x = -Inf, y = Inf, label = sprintf("Slope: %.3f", slope)),
            hjust = -0.1, vjust = 1.5, size = 4.5, family = base_family, inherit.aes = FALSE) +
  facet_wrap(~ ecoregion, ncol = 4) +
  scale_colour_manual(values = dark, guide = "none") +
  scale_fill_manual(values = dark, guide = "none") +
  scale_x_continuous(breaks = seq(1998, 2022, by = 6)) +
  labs(x = "Season (starting year)",
       y = expression(paste("Mean integrated Chl-", italic(a), " (mg ", m^{-3}, " season)"))) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(fig4, "Fig4_biomass_trajectories", width = 14, height = 8)


# ---- 5. Figures 5 and S9: phenological modes ----------------------------------------

t_fine <- seq(1, n_steps, length.out = 300)
B_fine <- eval.basis(t_fine, pca$basis)
mu     <- as.vector(B_fine %*% pca$mean_coef)
H      <- B_fine %*% pca$harmonics

# Background: random sample of individual normalized curves (season 2012/13)
set.seed(1)
X_bg  <- readRDS(file.path(dir_02, "data","chunks", "coefs_std_2012-2013.rds"))
bg    <- B_fine %*% X_bg[, sample(ncol(X_bg), 100)]
df_bg <- data.frame(t = rep(t_fine, ncol(bg)), value = as.vector(bg),
                    id = rep(seq_len(ncol(bg)), each = length(t_fine)))
rm(X_bg)

# Mean curve +/- 2 SD of each mode
plot_modes <- function(modes, ncol = length(modes), c_sd = 2) {
  d <- bind_rows(lapply(modes, function(k) {
    s <- c_sd * sqrt(pca$values[k])
    data.frame(t = t_fine, mode = sprintf("fPC%d (%.1f%%)", k, 100 * pca$varprop[k]),
               mean = mu, plus = mu + s * H[, k], minus = mu - s * H[, k])
  }))
  d$mode <- factor(d$mode, levels = unique(d$mode))
  ggplot(d, aes(t)) +
    geom_line(data = df_bg, aes(t, value, group = id), colour = "grey88",
              linewidth = 0.4, inherit.aes = FALSE) +
    geom_line(aes(y = plus),  colour = "#D33F6A", linewidth = 1) +
    geom_line(aes(y = minus), colour = "#4A6FE3", linewidth = 1) +
    geom_line(aes(y = mean),  colour = "black", linetype = "dashed", linewidth = 1) +
    facet_wrap(~ mode, ncol = ncol) +
    scale_x_continuous(breaks = month_breaks, labels = month_labels, expand = c(0, 0)) +
    labs(x = NULL, y = "Amplitude") +
    theme_pub
}

clim_map <- function(var, lab) {
  ok <- !is.na(clim[[var]])
  raster_map(to_polar(clim$cell[ok], clim[[var]][ok]), score_scale(clim[[var]]), lab,
             na_fill = "white") +
    colourbar()
}

fig5 <- (plot_modes(1) | plot_modes(2)) /
  (clim_map("mean_PC1", "Mean fPC1") | clim_map("mean_PC2", "Mean fPC2")) +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)")))
save_fig(fig5, "Fig5_phenological_modes", width = 14, height = 14)

save_fig(plot_modes(1:8, ncol = 4), "FigS9_phenological_modes_1-8", width = 16, height = 8)


# ---- 6. Figure 6: phenological trends ----------------------------------------------

sig_pc  <- lapply(trends[c("fPC1", "fPC2")], function(tr) filter(tr, signif %in% TRUE))
pc_scale <- slope_scale(c(sig_pc$fPC1$slope, sig_pc$fPC2$slope))   # common to both maps

fig6_map <- function(v) {
  raster_map(to_polar(sig_pc[[v]]$cell, sig_pc[[v]]$slope), pc_scale,
             bquote(.(v) ~ "slope (yr"^-1*")")) +
    colourbar()
}

fig6 <- (fig6_map("fPC1") | fig6_map("fPC2")) /
  (ridge_plot(trends$fPC1, expression("fPC1 slope (yr"^-1*")"), c(-0.5, 0.5)) |
     ridge_plot(trends$fPC2, expression("fPC2 slope (yr"^-1*")"), c(-0.5, 0.5))) +
  plot_layout(heights = c(1.4, 1)) +
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)")))
save_fig(fig6, "Fig6_phenology_trends", width = 14, height = 14)


# ---- 7. Figures 7, 8, S12 and S13: two-period comparison ----------------------------

period_lab <- c(early  = sprintf("%d–%d", tp$early_window[1],  tp$early_window[2]  + 1),
                recent = sprintf("%d–%d", tp$recent_window[1], tp$recent_window[2] + 1))

# All ecoregions, one mode, positive (top) and negative (bottom) trends
period_facets <- function(mode_name) {
  ps_mode <- filter(tp$period_summary, mode == mode_name)
  y_rng   <- range(c(ps_mode$lo, ps_mode$hi, tp$mean_curve$value), na.rm = TRUE)
  one_sign <- function(trend_label, recent_col) {
    ps <- filter(ps_mode, trend == trend_label)
    sh <- filter(tp$shares, mode == mode_name, trend == trend_label)
    ggplot(ps, aes(time)) +
      geom_line(data = tp$mean_curve, aes(time, value), colour = "grey75",
                linetype = "dashed", linewidth = 0.8, inherit.aes = FALSE) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = period), alpha = 0.35) +
      geom_line(aes(y = med, colour = period), linewidth = 1.1) +
      geom_text(data = sh, aes(label = sprintf("%.1f%%", share)), x = Inf, y = Inf,
                hjust = 1.1, vjust = 1.4, size = 4.5, colour = "grey30",
                family = base_family, inherit.aes = FALSE) +
      facet_wrap(~ ecoregion, ncol = 4) +
      coord_cartesian(ylim = y_rng) +
      scale_colour_manual(values = c(early = "grey55", recent = recent_col),
                          labels = period_lab, name = NULL) +
      scale_fill_manual(values = c(early = "grey55", recent = recent_col), guide = "none") +
      scale_x_continuous(breaks = month_breaks, labels = month_labels, expand = c(0, 0)) +
      labs(x = "Month", y = "Amplitude") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
  }
  one_sign("positive", "#E44748") / one_sign("negative", "#6796BF")
}

save_fig(period_facets("fPC1"), "FigS12_two_periods_fPC1", width = 14, height = 12)
save_fig(period_facets("fPC2"), "FigS13_two_periods_fPC2", width = 14, height = 12)

# Selected panels (Figures 7 and 8)
densify <- function(time, y, n = 500) spline(time, y, n = n)

period_panels <- function(sel, show_peaks = TRUE) {
  ps <- bind_rows(lapply(seq_len(nrow(sel)), function(i) {
    d <- filter(tp$period_summary, mode == sel$mode[i],
                ecoregion == sel$ecoregion[i], trend == sel$trend[i])
    bind_rows(lapply(split(d, d$period), function(p) data.frame(
      period = p$period[1],
      time   = densify(p$time, p$med)$x,
      med    = densify(p$time, p$med)$y,
      lo     = densify(p$time, p$lo)$y,
      hi     = densify(p$time, p$hi)$y))) |>
      mutate(panel = sel$panel[i], colour = sel$colour[i])
  }))
  pk <- bind_rows(lapply(seq_len(nrow(sel)), function(i)
    filter(tp$peaks, mode == sel$mode[i], ecoregion == sel$ecoregion[i],
           trend == sel$trend[i]) |> mutate(panel = sel$panel[i], colour = sel$colour[i])))
  sh <- bind_rows(lapply(seq_len(nrow(sel)), function(i)
    filter(tp$shares, mode == sel$mode[i], ecoregion == sel$ecoregion[i],
           trend == sel$trend[i]) |> mutate(panel = sel$panel[i])))
  for (d in c("ps", "pk", "sh")) {
    x <- get(d); x$panel <- factor(x$panel, levels = sel$panel); assign(d, x)
  }
  ps$key <- ifelse(ps$period == "early", "early", ps$colour)
  pk$key <- ifelse(pk$period == "early", "early", pk$colour)
  cols   <- c(early = "grey25", setNames(unique(sel$colour), unique(sel$colour)))
  mc     <- data.frame(densify(tp$mean_curve$time, tp$mean_curve$value))
  
  p <- ggplot(ps, aes(time)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey75", linewidth = 0.6) +
    geom_line(data = mc, aes(x, y), colour = "grey65", linetype = "longdash",
              linewidth = 0.9, inherit.aes = FALSE) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = key, group = period), alpha = 0.3) +
    geom_line(aes(y = med, colour = key, group = period), linewidth = 1.2) +
    geom_text(data = sh, aes(label = sprintf("%.1f%%", share)), x = -Inf, y = Inf,
              hjust = -0.25, vjust = 1.5, size = 6, colour = "grey30",
              family = base_family, inherit.aes = FALSE)
  if (show_peaks) {
    p <- p +
      geom_vline(data = pk, aes(xintercept = peak_time, colour = key),
                 linetype = "dashed", linewidth = 0.6) +
      geom_point(data = pk, aes(peak_time, peak_med), colour = "grey20", size = 4) +
      geom_point(data = pk, aes(peak_time, peak_med, colour = key), size = 2.4)
  }
  p +
    facet_wrap(~ panel, ncol = nrow(sel)) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_fill_manual(values = cols, guide = "none") +
    scale_x_continuous(breaks = month_breaks, labels = month_labels, expand = c(0, 0)) +
    labs(x = "Month", y = "Amplitude") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

sel7 <- data.frame(mode = "fPC2", ecoregion = c(6, 7), trend = "positive",
                   panel = c("Eco6 · fPC2 +", "Eco7 · fPC2 +"), colour = "#F76D5E")
sel8 <- data.frame(mode = "fPC1", ecoregion = 5, trend = c("positive", "negative"),
                   panel = c("Eco5 · fPC1 +", "Eco5 · fPC1 −"),
                   colour = c("#A50021", "#264CFF"))

save_fig(period_panels(sel7, show_peaks = TRUE),  "Fig7_two_periods_eco6_7", width = 11, height = 5.5)
save_fig(period_panels(sel8, show_peaks = FALSE), "Fig8_two_periods_eco5",   width = 11, height = 5.5)


# ---- 8. Supporting Information maps ------------------------------------------------

# Figure S1: missing eight-day composites, season 2012/13
n_valid <- readRDS(file.path(dir_02, "data","n_valid_2012-2013.rds"))
figS1 <- raster_map(to_polar(seq_along(n_valid), n_steps - n_valid, method = "near"),
                    scale_fill_stepsn(colours = rev(hcl.colors(11, "ag_Sunset")),
                                      breaks = seq(0, 22, by = 2), limits = c(0, 23),
                                      oob = scales::squish, na.value = "transparent"),
                    "Number of missing observations",
                    show_na = FALSE, boundaries = FALSE, labels = FALSE) +
  colourbar()
save_fig(figS1, "FigS1_missing_composites_2012-2013")

# Figures S6 and S11: significance of the trends
save_fig(sig_map(trends$biomass), "FigS6_significance_biomass")
save_fig(sig_map(trends$fPC1) + sig_map(trends$fPC2) +
           plot_annotation(tag_levels = list(c("(a)", "(b)"))),
         "FigS11_significance_phenology", width = 16, height = 8)

# Figure S7: productive-season Chl-a climatology (log scale)
ok <- !is.na(clim$chl_mean)
figS7 <- raster_map(to_polar(clim$cell[ok], clim$chl_mean[ok]),
                    scale_fill_gradientn(colours = rev(hcl.colors(10, "YlGnBu")),
                                         trans = "log10", breaks = c(0.05, 0.1, 0.3, 0.5, 1),
                                         oob = scales::squish, na.value = "transparent"),
                    expression(paste("Chl-", italic(a), " (mg ", m^{-3}, ")")),
                    show_na = FALSE, labels = FALSE) +
  colourbar()
save_fig(figS7, "FigS7_chl_climatology")

# Figure S8: relative change in integrated Chl-a
rel <- readRDS(file.path(dir_03, "relative_change.rds"))
figS8 <- raster_map(to_polar(rel$cell, rel$perc_change),
                    scale_fill_stepsn(colours = rev(hcl.colors(8, "Spectral")),
                                      breaks = seq(-200, 200, by = 50), limits = c(-200, 200),
                                      oob = scales::squish, na.value = "transparent"),
                    "Percentage change in biomass (%)") +
  colourbar()
save_fig(figS8, "FigS8_relative_change")

# Figure S10: mean fPC3 and fPC4 scores
save_fig(clim_map("mean_PC3", "Mean fPC3") + clim_map("mean_PC4", "Mean fPC4") +
           plot_annotation(tag_levels = list(c("(a)", "(b)"))),
         "FigS10_mean_fPC3_fPC4", width = 16, height = 8)

message("Done. Figures in ", fig_dir)