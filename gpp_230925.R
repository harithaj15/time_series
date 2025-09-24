# ============================================================
# GPP Resilience Analysis (Yao et al. 2024 workflow)
# Mean GPP, STL decomposition, AR(1), Mann–Kendall trend
# ============================================================

# -------------------------
# 0. Load packages
# -------------------------
pkgs <- c("terra", "stlplus", "Kendall", "viridis", "ggplot2", "dplyr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(terra)
library(stlplus)
library(Kendall)
library(viridis)
library(ggplot2)
library(dplyr)

# -------------------------
# 1. User inputs
# -------------------------
indir  <- "E:/haritha/PhD/resilience/GPP/gpp_yearly_stack_230925/western_ghats"  
outdir <- "E:/haritha/PhD/resilience/GPP/results_240925/western_ghats"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ncores <- 4
start_year <- 2001
end_year   <- 2023
years      <- start_year:end_year
freq       <- 12
nlayers    <- length(years) * freq

# -------------------------
# 2. Read raster stack
# -------------------------
cat("Reading GPP stack...\n")
gpp_files <- list.files(indir, pattern = "\\.tif$", full.names = TRUE)
gpp_files <- sort(gpp_files)
gpp_stack <- rast(gpp_files)
stopifnot(nlyr(gpp_stack) == nlayers)

time_seq <- seq(as.Date(sprintf("%04d-01-01", start_year)),
                by = "month", length.out = nlayers)
layer_year <- as.integer(format(time_seq, "%Y"))

# -------------------------
# 3. Mean GPP
# -------------------------
cat("Computing mean GPP across years...\n")
gpp_mean <- app(gpp_stack, mean, na.rm = TRUE, cores = ncores)
gpp_mean_wgs84 <- project(gpp_mean, "EPSG:4326")
writeRaster(gpp_mean_wgs84, file.path(outdir, "gpp_mean_2001_2023_wgs84.tif"), overwrite = TRUE)

# Annual mean GPP time series
gpp_yearly <- tapp(gpp_stack, index = layer_year, fun = mean, na.rm = TRUE)
gpp_mean_series <- global(gpp_yearly, mean, na.rm = TRUE)[, 1]
gpp_df <- data.frame(year = years, mean_GPP = gpp_mean_series)

# Plot mean GPP time series
p1 <- ggplot(gpp_df, aes(x = year, y = mean_GPP)) +
  geom_line(color = "forestgreen", linewidth = 1.2) +
  geom_point(color = "darkgreen") +
  labs(title = "Mean Annual GPP (2001–2023)", x = "Year", y = "Mean GPP") +
  theme_minimal()
ggsave(file.path(outdir, "line_mean_gpp_2001_2023.png"), plot = p1, width = 8, height = 5, dpi = 300)

p1



# -------------------------
# 3b. Pixel-wise GPP trends (2001–2023)
# -------------------------
cat("Computing pixel-wise GPP trends (2001–2023)...\n")

# ---- Mann–Kendall per pixel ----
mk_gpp_fun <- function(y) {
  if (all(is.na(y))) return(c(NA_real_, NA_real_))
  if (length(na.omit(y)) < 3) return(c(NA_real_, NA_real_))
  res <- tryCatch(Kendall::MannKendall(na.omit(y)), error = function(e) NULL)
  if (is.null(res)) return(c(NA_real_, NA_real_))
  return(c(as.numeric(res$tau), as.numeric(res$sl)))
}

mk_gpp_rasters <- app(
  gpp_yearly, fun = mk_gpp_fun, cores = ncores,
  filename = file.path(outdir, "gpp_mk_trend_2001_2023.tif"), overwrite = TRUE
)
names(mk_gpp_rasters) <- c("tau", "pval")

mk_gpp_rasters_wgs84 <- project(mk_gpp_rasters, "EPSG:4326")
writeRaster(mk_gpp_rasters_wgs84,
            file.path(outdir, "gpp_mk_trend_2001_2023_wgs84.tif"),
            overwrite = TRUE)

# ---- Linear regression slope per pixel ----
# Linear regression slope per pixel (ΔGPP per year)
lm_gpp_fun <- function(y) {
  if (all(is.na(y))) return(NA_real_)
  if (length(na.omit(y)) < 3) return(NA_real_)
  tt <- seq(start_year, end_year)[1:length(y)]  # build years internally
  res <- tryCatch({
    coef(lm(y ~ tt))[2]   # slope = change per year
  }, error = function(e) NA_real_)
  return(as.numeric(res))
}

slope_gpp_raster <- app(
  gpp_yearly, fun = lm_gpp_fun, cores = ncores,
  filename = file.path(outdir, "gpp_lm_slope_2001_2023.tif"), overwrite = TRUE
)

slope_gpp_raster_wgs84 <- project(slope_gpp_raster, "EPSG:4326")
writeRaster(slope_gpp_raster_wgs84,
            file.path(outdir, "gpp_lm_slope_2001_2023_wgs84.tif"),
            overwrite = TRUE)

# -------------------------
# 3c. Plot pixel-wise GPP trends
# -------------------------
# Convert rasters to data frames
gpp_tau_df  <- as.data.frame(mk_gpp_rasters[["tau"]], xy = TRUE, na.rm = TRUE)
gpp_pval_df <- as.data.frame(mk_gpp_rasters[["pval"]], xy = TRUE, na.rm = TRUE)
gpp_slope_df <- as.data.frame(slope_gpp_raster, xy = TRUE, na.rm = TRUE)
names(gpp_slope_df)[3] <- "slope"

# Tau (non-parametric trend strength/direction)
p_gpp_tau <- ggplot(gpp_tau_df, aes(x = x, y = y, fill = tau)) +
  geom_tile() +
  scale_fill_viridis(option = "D", name = "Tau", limits = c(-1, 1)) +
  coord_equal() +
  labs(
    title = "Pixel-wise GPP Trends (2001–2023, MK Tau)",
    subtitle = "Positive = increasing GPP, Negative = decreasing GPP"
  ) +
  theme_minimal()

# P-values (significance of MK trend)
p_gpp_pval <- ggplot(gpp_pval_df, aes(x = x, y = y, fill = pval)) +
  geom_tile() +
  scale_fill_viridis(option = "C", direction = -1, name = "p-value") +
  coord_equal() +
  labs(
    title = "Pixel-wise GPP Trend Significance (2001–2023, MK test)",
    subtitle = "Darker = more significant"
  ) +
  theme_minimal()

# Regression slope (ΔGPP per year)
p_gpp_slope <- ggplot(gpp_slope_df, aes(x = x, y = y, fill = slope)) +
  geom_tile() +
  scale_fill_viridis(option = "B", name = "Slope\n(GPP/year)") +
  coord_equal() +
  labs(
    title = "Pixel-wise GPP Linear Trend (2001–2023)",
    subtitle = "Units: GPP per year"
  ) +
  theme_minimal()

# Save plots
ggsave(file.path(outdir, "map_gpp_trends_tau.png"),   plot = p_gpp_tau,   width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "map_gpp_trends_pval.png"),  plot = p_gpp_pval,  width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "map_gpp_trends_slope.png"), plot = p_gpp_slope, width = 7, height = 5, dpi = 300)




# -------------------------
# 4. STL decomposition → residuals
# -------------------------
cat("Running STL decomposition per pixel...\n")

stlplus_fun <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  if (sum(!is.na(x)) < 24) return(rep(NA_real_, length(x)))
  res <- tryCatch({
    fit <- stlplus::stlplus(x, n.p = 12, s.window = "periodic")
    as.numeric(fit$data$remainder)
  }, error = function(e) rep(NA_real_, length(x)))
  return(res)
}

resid_stack <- app(gpp_stack, fun = stlplus_fun, cores = ncores,
                   filename = file.path(outdir, "residuals_monthly_2001_2023.tif"),
                   overwrite = TRUE)
resid_stack_wgs84 <- project(resid_stack, "EPSG:4326")
writeRaster(resid_stack_wgs84, file.path(outdir, "residuals_monthly_2001_2023_wgs84.tif"), overwrite = TRUE)

# -------------------------
# 5. AR(1) calculation (5-year sliding windows)
# -------------------------
cat("Computing AR(1) with 5-year sliding window...\n")

ar1_fun <- function(x) {
  win_size <- 60   # 5 years = 60 months
  n_win <- length(x) - win_size + 1
  out <- rep(NA_real_, n_win)
  
  ar1_single <- function(y) {
    if (all(is.na(y))) return(NA_real_)
    if (length(na.omit(y)) < 2) return(NA_real_)
    tryCatch(acf(y, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2],
             error = function(e) NA_real_)
  }
  
  for (i in seq_len(n_win)) {
    seg <- x[i:(i + win_size - 1)]
    out[i] <- ar1_single(seg)
  }
  return(out)
}

ar1_stack <- app(resid_stack, fun = ar1_fun, cores = ncores,
                 filename = file.path(outdir, "ar1_sliding_2001_2023.tif"),
                 overwrite = TRUE)
ar1_stack_wgs84 <- project(ar1_stack, "EPSG:4326")
writeRaster(ar1_stack_wgs84, file.path(outdir, "ar1_sliding_2001_2023_wgs84.tif"), overwrite = TRUE)

# -------------------------
# 6. Mann–Kendall on AR(1) time series (per pixel)
# -------------------------
cat("Running Mann–Kendall test on AR1 series...\n")

mk_fun <- function(y) {
  if (all(is.na(y))) return(c(NA_real_, NA_real_))
  if (length(na.omit(y)) < 3) return(c(NA_real_, NA_real_))
  res <- tryCatch(Kendall::MannKendall(na.omit(y)), error = function(e) NULL)
  if (is.null(res)) return(c(NA_real_, NA_real_))
  return(c(as.numeric(res$tau), as.numeric(res$sl)))
}

mk_rasters <- app(ar1_stack, fun = mk_fun, cores = ncores,
                  filename = file.path(outdir, "ar1_mk_trend.tif"), overwrite = TRUE)
names(mk_rasters) <- c("tau", "pval")
mk_rasters_wgs84 <- project(mk_rasters, "EPSG:4326")
writeRaster(mk_rasters_wgs84, file.path(outdir, "ar1_mk_trend_wgs84.tif"), overwrite = TRUE)

# -------------------------
# 7. Global mean AR1 time series
# -------------------------
cat("Computing mean AR1 time series...\n")
win_size <- 60
mid_dates <- time_seq[seq(win_size/2, nlayers - win_size/2)]
mid_years <- as.integer(format(mid_dates, "%Y"))

ar1_mean_series <- global(ar1_stack, mean, na.rm = TRUE)[, 1]
ar1_df <- data.frame(
  year = mid_years[1:length(ar1_mean_series)],
  mean_AR1 = ar1_mean_series
)

# -------------------------
# 7a. Mann–Kendall test on global mean AR1
# -------------------------
mk_global <- MannKendall(na.omit(ar1_df$mean_AR1))
print(mk_global)

p2 <- ggplot(ar1_df, aes(x = year, y = mean_AR1)) +
  geom_line(color = "firebrick", linewidth = 1.2) +
  geom_point(color = "darkred") +
  labs(
    title = "Global Mean AR1 (5-yr sliding windows)",
    subtitle = sprintf("Mann–Kendall: tau = %.3f, p = %.3f",
                       mk_global$tau, mk_global$sl),
    x = "Year", y = "Mean AR1"
  ) +
  theme_minimal()
ggsave(file.path(outdir, "line_mean_ar1_mk_test.png"),
       plot = p2, width = 8, height = 5, dpi = 300)

# -------------------------
# 7b. Linear regression trend on global AR1
# -------------------------
lm_fit <- lm(mean_AR1 ~ year, data = ar1_df)
summary_fit <- summary(lm_fit)
slope <- coef(lm_fit)[2]
pval  <- summary_fit$coefficients[2, 4]
cat("Linear trend in mean AR1:", slope, "per year; p =", pval, "\n")

p3 <- ggplot(ar1_df, aes(x = year, y = mean_AR1)) +
  geom_line(color = "firebrick", linewidth = 1.2) +
  geom_point(color = "darkred") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    title = "Mean AR1 of Residuals (5-yr sliding windows)",
    subtitle = sprintf("Trend slope = %.4f per year (p = %.3f)", slope, pval),
    x = "Year", y = "Mean AR1"
  ) +
  theme_minimal()
p3
ggsave(file.path(outdir, "line_mean_ar1_with_trend_2001_2023.png"),
       plot = p3, width = 8, height = 5, dpi = 300)

# -------------------------
# 8. Maps of tau and p-values
# -------------------------
cat("Plotting spatial trend maps...\n")
tau_df  <- as.data.frame(mk_rasters[["tau"]], xy = TRUE, na.rm = TRUE)
pval_df <- as.data.frame(mk_rasters[["pval"]], xy = TRUE, na.rm = TRUE)

p_tau <- ggplot(tau_df, aes(x = x, y = y, fill = tau)) +
  geom_tile() +
  scale_fill_viridis(option = "D", name = "Tau", limits = c(-1, 1)) +
  coord_equal() +
  labs(
    title = "Resilience Trends (Mann–Kendall Tau)",
    subtitle = "Negative = increasing resilience (AR1 decreasing), Positive = declining resilience"
  ) +
  theme_minimal()

p_pval <- ggplot(pval_df, aes(x = x, y = y, fill = pval)) +
  geom_tile() +
  scale_fill_viridis(option = "C", direction = -1, name = "p-value") +
  coord_equal() +
  labs(
    title = "Resilience Trend Significance (Mann–Kendall p-values)",
    subtitle = "Darker = more significant"
  ) +
  theme_minimal()

p_tau
p_pval

ggsave(file.path(outdir, "map_resilience_trends_tau.png"), plot = p_tau, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "map_resilience_trends_pval.png"), plot = p_pval, width = 7, height = 5, dpi = 300)

# -------------------------
# Done
# -------------------------
cat("All results saved to:", outdir, "\n")
