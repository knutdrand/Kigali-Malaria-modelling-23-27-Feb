#!/usr/bin/env Rscript
# =============================================================================
# run_final_model.R
# Standalone batch script for the final BYM + RW1 + IID space-time interaction
# malaria model with lagged climate covariates, using R-INLA.
#
# Usage:
#   Rscript run_final_model.R data.csv shapefile.geojson [output_dir]
#
# Arguments:
#   data.csv          Merged CSV with columns: District, Sector, Year, Month,
#                     Malaria_Cases, Population, Max_temperature, Min_temperature,
#                     Precipitation, Relative_humidity
#   shapefile.geojson Rwanda ADM3 sector polygons (GeoJSON)
#   output_dir        Output directory (default: ./output)
# =============================================================================

# ---------------------------------------------------------------------------
# 1) Parse command-line arguments & load packages
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript run_final_model.R data.csv shapefile.geojson [output_dir]\n",
    call. = FALSE
  )
}

data_path <- args[1]
shp_path  <- args[2]
output_dir <- if (length(args) >= 3) args[3] else "./output"

if (!file.exists(data_path)) stop("Data file not found: ", data_path, call. = FALSE)
if (!file.exists(shp_path))  stop("Shapefile not found: ", shp_path, call. = FALSE)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading packages...\n")
pacman::p_load(
  dplyr, sf, tidyr, ggplot2, stringr,
  INLA, spdep, viridis, tibble
)

options(stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# 2) Read inputs
# ---------------------------------------------------------------------------
cat("Reading data from:", data_path, "\n")
csv_data <- read.csv(data_path, stringsAsFactors = FALSE)

cat("Reading shapefile from:", shp_path, "\n")
Rwa_data <- st_read(shp_path, quiet = TRUE) %>%
  st_make_valid()

# ---------------------------------------------------------------------------
# 3) Spatial prep: fix sector names, build adjacency graph
# ---------------------------------------------------------------------------
Rwa_data <- Rwa_data %>%
  mutate(Sec_ID = row_number()) %>%
  mutate(
    ADM3_EN = case_when(
      ADM3_EN == "Mageregere" ~ "Mageragere",
      ADM3_EN == "Shyrongi"   ~ "Shyorongi",
      ADM3_EN == "Ririma"     ~ "Rilima",
      TRUE                    ~ ADM3_EN
    )
  )

graph_file <- file.path(output_dir, "Adj_Map.graph")
nb_obj <- poly2nb(Rwa_data, queen = FALSE)
nb2INLA(graph_file, nb_obj)
Rwa_adj <- normalizePath(graph_file)

cat("Adjacency graph written to:", graph_file, "\n")
cat("Number of sectors:", nrow(Rwa_data), "\n")

# ---------------------------------------------------------------------------
# 4) Data merge: standardize keys, join CSV to geometry
# ---------------------------------------------------------------------------
csv_data <- csv_data %>%
  mutate(
    district = str_to_sentence(str_squish(as.character(District))),
    sector   = str_to_sentence(str_squish(as.character(Sector))),
    Year     = as.integer(Year),
    Month    = as.integer(Month)
  )

Geo_data <- merge(
  Rwa_data, csv_data,
  by.x = c("ADM3_EN", "ADM2_EN"),
  by.y = c("sector", "district"),
  all.x = FALSE
)

n_matched <- nrow(Geo_data)
n_sectors_matched <- length(unique(Geo_data$Sec_ID))
cat("Rows after spatial join:", n_matched, "\n")
cat("Sectors matched:", n_sectors_matched, "of", nrow(Rwa_data), "\n")

if (n_matched == 0) {
  stop("No rows matched during spatial join. Check District/Sector names.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# 5) Feature engineering: indices, lags, standardization
# ---------------------------------------------------------------------------
dat <- Geo_data %>%
  arrange(Sec_ID, Year, Month) %>%
  mutate(
    ID            = Sec_ID,
    Malaria_Cases = replace_na(Malaria_Cases, 0),
    ID.time       = rep(
      seq_len(n_distinct(Year) * 12),
      times = n_distinct(Sec_ID)
    ),
    ID.space.time = row_number(),
    offset_pop    = log(Population / 1000)
  )

# Rename climate columns to safe names (handle both possible column name styles)
rename_map <- c(
  "Precipitation..ERA5.Land." = "Prec",
  "Precipitation (ERA5-Land)" = "Prec",
  "Precipitation"             = "Prec",
  "Min_.temperature"          = "Min_temperature",
  "Min_ temperature"          = "Min_temperature",
  "Relative.humidity..ERA5.Land." = "Relative_humidity",
  "Relative humidity (ERA5-Land)" = "Relative_humidity"
)

for (old_name in names(rename_map)) {
  if (old_name %in% names(dat)) {
    names(dat)[names(dat) == old_name] <- rename_map[[old_name]]
  }
}

# Create lagged climate variables (by sector)
dat <- dat %>%
  arrange(Sec_ID, Year, Month) %>%
  group_by(Sec_ID) %>%
  mutate(
    Max_temp_lag1 = lag(Max_temperature,   1),
    Max_temp_lag2 = lag(Max_temperature,   2),
    Prec_lag1     = lag(Prec,              1),
    Prec_lag2     = lag(Prec,              2),
    Min_temp_lag1 = lag(Min_temperature,   1),
    Min_temp_lag2 = lag(Min_temperature,   2),
    RH_lag1       = lag(Relative_humidity, 1),
    RH_lag2       = lag(Relative_humidity, 2)
  ) %>%
  ungroup()

lag_vars <- c(
  "Max_temp_lag1", "Max_temp_lag2",
  "Min_temp_lag1", "Min_temp_lag2",
  "Prec_lag1",     "Prec_lag2",
  "RH_lag1",       "RH_lag2"
)

# Convert to numeric safely
to_num <- function(x) as.numeric(as.character(x))
dat <- dat %>%
  mutate(across(all_of(lag_vars), to_num))

# Drop rows with missing lags (first 2 months per sector)
n_before <- nrow(dat)
dat <- dat %>%
  filter(if_all(all_of(lag_vars), ~ !is.na(.x)))
cat("Rows dropped due to missing lags:", n_before - nrow(dat), "\n")

# Re-index after filtering
dat <- dat %>%
  arrange(Sec_ID, Year, Month) %>%
  mutate(
    ID.time       = as.integer(factor(paste(Year, Month), levels = unique(paste(Year, Month)))),
    ID.space.time = row_number()
  )

# Standardize lag variables
stdz <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
dat <- dat %>%
  mutate(across(all_of(lag_vars), stdz))

cat("Final modelling dataset:", nrow(dat), "rows,",
    n_distinct(dat$ID), "sectors,",
    max(dat$ID.time), "time points\n")

# ---------------------------------------------------------------------------
# 6) Fit final model: BYM + RW1 + IID space-time + lagged climate
# ---------------------------------------------------------------------------
cat("\nFitting final model (BYM + RW1 + IID space-time + 8 lagged climate covariates)...\n")

hyper_bym <- list(
  prec.unstruct = list(prior = "pc.prec", param = c(1, 0.01)),
  prec.spatial  = list(prior = "pc.prec", param = c(1, 0.01))
)

formula_final_bym <- as.integer(Malaria_Cases) ~ 1 +
  Max_temp_lag1 + Max_temp_lag2 +
  Min_temp_lag1 + Min_temp_lag2 +
  Prec_lag1 + Prec_lag2 +
  RH_lag1 + RH_lag2 +
  f(ID, model = "bym", graph = Rwa_adj, scale.model = TRUE, constr = TRUE, hyper = hyper_bym) +
  f(ID.time, model = "rw1", constr = TRUE, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(ID.space.time, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

fit_final <- inla(
  formula           = formula_final_bym,
  family            = "poisson",
  data              = dat,
  offset            = offset_pop,
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE, link = 1)
)

cat("\n=== Model Fit Summary ===\n")
cat("WAIC:", fit_final$waic$waic, "\n")
cat("DIC: ", fit_final$dic$dic, "\n")

# Fixed effects table
fixef <- fit_final$summary.fixed %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  mutate(
    RR     = exp(mean),
    RR_LCL = exp(`0.025quant`),
    RR_UCL = exp(`0.975quant`)
  )

fixef_climate <- fixef %>%
  filter(term %in% lag_vars) %>%
  select(term, mean, `0.025quant`, `0.975quant`, RR, RR_LCL, RR_UCL)

cat("\n=== Fixed Effects (Relative Risk) ===\n")
print(fixef_climate, row.names = FALSE)

# ---------------------------------------------------------------------------
# 7) Post-processing: spatial RR, exceedance P(RR>1), temporal RR
# ---------------------------------------------------------------------------
cat("\nComputing posterior summaries...\n")

# --- Spatial RR ---
spat_df <- as.data.frame(fit_final$summary.random$ID)
if (!("ID" %in% names(spat_df))) {
  spat_df <- spat_df %>% rownames_to_column("ID") %>% mutate(ID = as.integer(ID))
} else {
  spat_df <- spat_df %>% mutate(ID = as.integer(ID))
}

spat_out <- spat_df %>%
  mutate(
    RR_spatial     = exp(mean),
    RR_spatial_LCL = exp(`0.025quant`),
    RR_spatial_UCL = exp(`0.975quant`)
  ) %>%
  select(ID, RR_spatial, RR_spatial_LCL, RR_spatial_UCL)

# Exceedance probability: P(RR > 1)
spat_out$PP_RR_gt_1 <- sapply(seq_along(fit_final$marginals.random$ID), function(k) {
  1 - inla.pmarginal(q = 0, marginal = fit_final$marginals.random$ID[[k]])
})

# Join to one geometry per sector
geo_sec <- dat %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup() %>%
  select(ID, ADM3_EN, ADM2_EN, geometry)

map_spatial <- geo_sec %>%
  left_join(spat_out, by = "ID")

# --- Temporal RR ---
temp_out <- fit_final$summary.random$ID.time %>%
  as.data.frame() %>%
  rownames_to_column("ID.time") %>%
  mutate(
    ID.time     = as.integer(ID.time),
    RR_time     = exp(mean),
    RR_time_LCL = exp(`0.025quant`),
    RR_time_UCL = exp(`0.975quant`)
  ) %>%
  select(ID.time, RR_time, RR_time_LCL, RR_time_UCL)

time_lookup <- dat %>%
  st_drop_geometry() %>%
  distinct(ID.time, Year, Month) %>%
  arrange(ID.time)

temp_out <- temp_out %>%
  left_join(time_lookup, by = "ID.time") %>%
  arrange(ID.time)

# --- Space-time RR ---
st_out <- NULL
if (!is.null(fit_final$summary.random$ID.space.time)) {
  st_out <- fit_final$summary.random$ID.space.time %>%
    as.data.frame() %>%
    rownames_to_column("ID.space.time") %>%
    mutate(
      ID.space.time = as.integer(ID.space.time),
      RR_st     = exp(mean),
      RR_st_LCL = exp(`0.025quant`),
      RR_st_UCL = exp(`0.975quant`)
    ) %>%
    select(ID.space.time, RR_st, RR_st_LCL, RR_st_UCL)
}

# ---------------------------------------------------------------------------
# 8) Forecasting: 3-month-ahead via INLA missing-data trick
# ---------------------------------------------------------------------------
cat("Generating 3-month-ahead forecast...\n")

T_last <- max(dat$ID.time)
h <- 3

future <- expand.grid(
  ID      = sort(unique(dat$ID)),
  ID.time = (T_last + 1):(T_last + h)
)

# Use last observed population per sector
pop_last <- dat %>%
  st_drop_geometry() %>%
  group_by(ID) %>%
  filter(ID.time == T_last) %>%
  ungroup() %>%
  select(ID, Population, offset_pop)

future <- future %>%
  left_join(pop_last, by = "ID")

# Climate persistence assumption: use last observed lag values
clim_last <- dat %>%
  st_drop_geometry() %>%
  group_by(ID) %>%
  filter(ID.time == T_last) %>%
  ungroup() %>%
  select(ID, all_of(lag_vars))

future <- future %>%
  left_join(clim_last, by = "ID")

# Space-time index for forecast rows
max_st <- max(dat$ID.space.time)
future <- future %>%
  arrange(ID, ID.time) %>%
  mutate(ID.space.time = max_st + row_number())

# Combine observed + future
dat_forecast <- bind_rows(
  dat %>% mutate(forecast = 0L),
  future %>% mutate(
    Malaria_Cases = NA_integer_,
    forecast      = 1L
  )
)

# Re-fit with missing outcomes (INLA prediction trick)
cat("Re-fitting model with forecast horizon...\n")
fit_forecast <- inla(
  formula           = formula_final_bym,
  family            = "poisson",
  data              = dat_forecast,
  offset            = offset_pop,
  control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE, link = 1)
)

# Extract forecasts
forecast_out <- dat_forecast %>%
  mutate(
    mu_mean = fit_forecast$summary.fitted.values$mean,
    mu_lcl  = fit_forecast$summary.fitted.values$`0.025quant`,
    mu_ucl  = fit_forecast$summary.fitted.values$`0.975quant`
  ) %>%
  filter(forecast == 1L) %>%
  arrange(ID, ID.time)

# Add calendar labels
last_cal <- dat %>%
  st_drop_geometry() %>%
  filter(ID.time == T_last) %>%
  distinct(Year, Month)

forecast_out <- forecast_out %>%
  mutate(
    Month = ((last_cal$Month + ID.time - T_last - 1) %% 12) + 1,
    Year  = last_cal$Year + (last_cal$Month + ID.time - T_last - 1) %/% 12
  )

cat("Forecast months:",
    paste(unique(paste(forecast_out$Year, forecast_out$Month, sep = "-")), collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# 9) Save outputs
# ---------------------------------------------------------------------------
cat("\nSaving outputs to:", output_dir, "\n")

saveRDS(spat_out,     file.path(output_dir, "post_spatial_RR.rds"))
saveRDS(temp_out,     file.path(output_dir, "post_temporal_RR.rds"))
if (!is.null(st_out)) saveRDS(st_out, file.path(output_dir, "post_spacetime_RR.rds"))
saveRDS(forecast_out, file.path(output_dir, "forecast_3month.rds"))

cat("  RDS files saved.\n")

# --- Plot: Spatial RR map ---
p_rr <- ggplot(map_spatial) +
  geom_sf(aes(fill = pmin(RR_spatial, 3)), color = NA) +
  scale_fill_viridis_c(option = "C", name = "RR", breaks = c(0.5, 1, 2, 3)) +
  labs(title = "Spatial Relative Risk (RR)") +
  theme_minimal()

ggsave(file.path(output_dir, "spatial_RR_map.png"), p_rr,
       width = 10, height = 8, dpi = 300)

# --- Plot: Exceedance probability map ---
p_pp <- ggplot(map_spatial) +
  geom_sf(aes(fill = PP_RR_gt_1), color = NA) +
  scale_fill_viridis_c(option = "C", name = "P(RR > 1)", limits = c(0, 1)) +
  labs(title = "Posterior Exceedance Probability P(RR > 1)") +
  theme_minimal()

ggsave(file.path(output_dir, "exceedance_probability_map.png"), p_pp,
       width = 10, height = 8, dpi = 300)

# --- Plot: Temporal RR trend ---
p_temp <- ggplot(temp_out, aes(x = ID.time)) +
  geom_ribbon(aes(ymin = RR_time_LCL, ymax = RR_time_UCL), alpha = 0.2) +
  geom_line(aes(y = RR_time), linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "Time (months)", y = "Relative Risk",
    title = "Temporal Relative Risk Trend"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "temporal_RR_trend.png"), p_temp,
       width = 12, height = 5, dpi = 300)

# --- Plot: Forecast vs historical (national) ---
dat0 <- st_drop_geometry(dat)
dat0$mu_hat <- fit_final$summary.fitted.values$mean
dat0$mu_lcl <- fit_final$summary.fitted.values$`0.025quant`
dat0$mu_ucl <- fit_final$summary.fitted.values$`0.975quant`

hist_nat <- dat0 %>%
  group_by(ID.time) %>%
  summarise(
    Year    = first(Year),
    Month   = first(Month),
    obs     = sum(Malaria_Cases, na.rm = TRUE),
    fit     = sum(mu_hat, na.rm = TRUE),
    fit_lcl = sum(mu_lcl, na.rm = TRUE),
    fit_ucl = sum(mu_ucl, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(ID.time)

fc_nat <- forecast_out %>%
  st_drop_geometry() %>%
  group_by(ID.time) %>%
  summarise(
    Year    = first(Year),
    Month   = first(Month),
    fc_mean = sum(mu_mean, na.rm = TRUE),
    fc_lcl  = sum(mu_lcl,  na.rm = TRUE),
    fc_ucl  = sum(mu_ucl,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(ID.time)

all_time <- sort(unique(c(hist_nat$ID.time, fc_nat$ID.time)))
time_map <- data.frame(ID.time = all_time, t_plot = seq_along(all_time))
hist_nat <- hist_nat %>% left_join(time_map, by = "ID.time")
fc_nat   <- fc_nat   %>% left_join(time_map, by = "ID.time")

p_forecast <- ggplot() +
  geom_ribbon(data = hist_nat, aes(x = t_plot, ymin = fit_lcl, ymax = fit_ucl),
              fill = "blue", alpha = 0.18) +
  geom_line(data = hist_nat, aes(x = t_plot, y = fit),
            color = "blue", linewidth = 0.9) +
  geom_line(data = hist_nat, aes(x = t_plot, y = obs),
            color = "blue", linetype = "dashed", linewidth = 0.9) +
  geom_ribbon(data = fc_nat, aes(x = t_plot, ymin = fc_lcl, ymax = fc_ucl),
              fill = "red", alpha = 0.22) +
  geom_line(data = fc_nat, aes(x = t_plot, y = fc_mean),
            color = "red", linewidth = 1.1) +
  geom_vline(xintercept = max(hist_nat$t_plot), linetype = "dotted") +
  labs(
    x = "Time (months)",
    y = "Malaria cases (national total)",
    title = "Historical (blue) and 3-month forecast (red) malaria trends",
    subtitle = "Solid = fitted/forecast mean; Dashed = observed; Shaded = 95% credible interval"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "forecast_vs_historical.png"), p_forecast,
       width = 14, height = 6, dpi = 300)

cat("  PNG plots saved.\n")
cat("\nDone.\n")
