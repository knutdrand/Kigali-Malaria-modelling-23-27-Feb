# ============================================================================
# MALARIA MODELLING PIPELINE (Scientific + Logical Order)
# Rwanda sector-level spatio-temporal malaria modelling with climate covariates
# ============================================================================

# ----------------------------------------------------------------------------
# 1) SETUP: Packages, options, and paths
# ----------------------------------------------------------------------------
pacman::p_load(
  dplyr, sf, tidyr, ggplot2, lubridate, ggtext,
  INLA, spdep, readxl, rio, here, stringr,
  INLAOutputs, viridis, patchwork, rnaturalearthdata,
  plotly, gganimate, leaflet, rnaturalearth,
  RColorBrewer, gifski, forecast, shiny, tibble
  
)

# (Optional) Global options
options(stringsAsFactors = FALSE)

# Key input files
shp_path   <- "shapefile/rwa_adm3_2006_NISR_WGS1984_20181002.shp"
u5_path    <- "U_5_diarrhea_climate.xlsx"
mal_path   <- "malaria_modelling_Oslo.xlsx"
graph_file <- "Adj_Map.graph"

# ----------------------------------------------------------------------------
# 2) SPATIAL DATA: Read polygons, clean names, create IDs, adjacency graph
# ----------------------------------------------------------------------------
# 2.1 Load shapefile
Rwa_data <- st_read(shp_path, quiet = FALSE) %>%
  st_make_valid()

# 2.2 Assign unique sector IDs + fix known spelling issues
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

# 2.3 Build adjacency graph for INLA
nb_obj <- poly2nb(Rwa_data, queen = FALSE)
nb2INLA(graph_file, nb_obj)
Rwa_adj <- file.path(getwd(), graph_file)

# 2.4 Quick check plot
ggplot(Rwa_data) +
  geom_sf(color = "blue", fill = "white") +
  coord_sf() + theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12))


# ----------------------------------------------------------------------------
# 3) DATA INGESTION: malaria datasets, clean keys
# ----------------------------------------------------------------------------
U5_raw      <- read_excel(u5_path)
Malaria_raw <- read_excel(mal_path)

# 3.1 Clean malaria sector/district names
Malaria_raw <- Malaria_raw %>%
  mutate(
    Sector   = str_to_sentence(str_replace(Sector, "\\s*\\([^\\)]+\\)", "")),
    District = str_to_sentence(District)
  )


# ----------------------------------------------------------------------------
# 4) DATA CLEANING & IMPUTATION (U5): handle missing values, standardize keys
# ----------------------------------------------------------------------------
U5_clean <- U5_raw %>%
  group_by(District, Year) %>%
  mutate(
    Population      = if_else(is.na(Population),      median(Population, na.rm = TRUE),      Population),
    U_5_Population  = if_else(is.na(U_5_Population),  median(U_5_Population, na.rm = TRUE),  U_5_Population),
    U_5_Ratio       = if_else(is.na(U_5_Ratio),       median(U_5_Ratio, na.rm = TRUE),       U_5_Ratio),
    Diarrhoea_cases = if_else(is.na(Diarrhoea_cases), median(Diarrhoea_cases, na.rm = TRUE), Diarrhoea_cases)
  ) %>%
  ungroup() %>%
  rename(district = District, sector = Sector) %>%
  mutate(
    district = str_to_lower(str_squish(district)),
    sector   = str_to_lower(str_squish(sector)),
    Year     = as.integer(Year),
    Month    = as.integer(Month)
  )


# ----------------------------------------------------------------------------
# 5) MALARIA SUBSET: keep needed columns + standardize keys
# ----------------------------------------------------------------------------
Malaria_sub <- Malaria_raw %>%
  select(District, Sector, Year, Month, `Malaria Cases`) %>%
  mutate(
    district = str_to_lower(str_squish(District)),
    sector   = str_to_lower(str_squish(Sector))
  ) %>%
  rename(Malaria_Cases = `Malaria Cases`) %>%
  select(district, sector, Year, Month, Malaria_Cases)


# ----------------------------------------------------------------------------
# 6) MERGE (U5 + MALARIA) and JOIN WITH GEOMETRY
# ----------------------------------------------------------------------------
Merged <- left_join(
  U5_clean, Malaria_sub,
  by = c("district", "sector", "Year", "Month")
) %>%
  replace_na(list(
    Malaria_Cases   = 0,
    Diarrhoea_cases = 0
  )) %>%
  mutate(
    district = str_to_sentence(district),
    sector   = str_to_sentence(sector)
  )

# 6.1 Join to polygons (ensure keys match shapefile fields)
Geo_data <- merge(
  Rwa_data, Merged,
  by.x = c("ADM3_EN", "ADM2_EN"),
  by.y = c("sector", "district"),
  all.x = FALSE
)


# ----------------------------------------------------------------------------
# 7) MODELLING DATASET: sort, indices, rename covariates, create lags
# ----------------------------------------------------------------------------
dat <- Geo_data %>%
  arrange(Sec_ID, Year, Month) %>%
  mutate(
    ID             = Sec_ID,
    Malaria_Cases  = replace_na(Malaria_Cases, 0),
    Diarrhoea_cases= replace_na(Diarrhoea_cases, 0),
    ID.time        = rep(1:(10*12), times = n_distinct(Sec_ID)),
    ID.space.time  = row_number()
  )

# 7.1 Rename climate columns to syntactically safe names
dat <- dat %>%
  rename(
    Prec              = `Precipitation (ERA5-Land)`,
    Min_temperature   = `Min_ temperature`,
    Relative_humidity = `Relative humidity (ERA5-Land)`
  )

# 7.2 Create lagged climate variables (by sector)
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
    RH_lag2       = lag(Relative_humidity, 2),
    Air_temp_lag1 = lag(Air_Temperature,   1),
    Air_temp_lag2 = lag(Air_Temperature,   2)
  ) %>%
  ungroup()

# ------------------------- MULTICOLLINEARITY CHECK -----------------------
df_clim <- dat %>% st_drop_geometry() %>%
  select(Max_temperature, Prec, Min_temperature, Relative_humidity) %>%
  na.omit()

#df_clim <- dat %>% st_drop_geometry() %>%
# select(Max_temperature, Prec, Min_temperature, Relative_humidity, Air_Temperature) %>%
#  na.omit()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# adjust margins so the title has room
par(mar = c(0, 0, 2, 0))

# draw the corrplot with a title
corrplot::corrplot(
  cor(df_clim),
  method      = "color",
  addCoef.col = "black",
  tl.cex      = 0.8,
  mar         = c(0, 0, 1, 0)
)
title("Correlation Matrix of Climate Covariates")

# VIF diagnostics
df_num <- as.data.frame(lapply(df_clim, as.numeric))
vif_res <- usdm::vif(df_num)
print(vif_res)

# Du tu high colleration between Air temperature and Max temperature, 
# we need to drop Air temperature
# ----------------------------------------------------------------------------
# 8) MODEL SPECIFICATION: priors + formulas
# ----------------------------------------------------------------------------
# 8.1 Hyperpriors for BYM2 (PC priors)
hyper_bym2 <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)),
  phi  = list(prior = "pc",      param = c(0.5, 0.5))
)
dat <- dat %>%
  mutate(
    offset_pop = log(Population / 1000)
  )

#1.4 BYM2 (recommended modern reparameterisation)
hyper_bym2 <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)),
  phi  = list(prior = "pc",      param = c(0.5, 0.5))
)

# 8.3 BYM + concurrent climate
formula_bym_climate <- as.integer(Malaria_Cases) ~ 1 +
  scale(Max_temperature) + scale(Prec) + scale(Min_temperature) +
  scale(Relative_humidity) +
  f(ID, model = "bym", graph = Rwa_adj, scale.model = TRUE, constr = TRUE) +
  f(ID.time, model = "rw1", constr = TRUE, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(ID.space.time, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

# 8.4 BYM + lagged climate
formula_bym_climate_lags <- as.integer(Malaria_Cases) ~ 1 +
  scale(Max_temp_lag1) + scale(Max_temp_lag2) +
  scale(Prec_lag1) + scale(Prec_lag2) +
  scale(Min_temp_lag1) + scale(Min_temp_lag2) +
  scale(RH_lag1) + scale(RH_lag2)  +
  f(ID, model = "bym", graph = Rwa_adj, scale.model = TRUE, constr = TRUE) +
  f(ID.time, model = "rw1", constr = TRUE, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  f(ID.space.time, model = "iid", constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
# ----------------------------------------------------------------------------
# 9) MODEL FITTING (INLA)


spatial_model_clima <- inla(
  formula           = formula_bym_climate,
  family            = "poisson",
  data              = dat,
  offset            = offset_pop,
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)
)

# (Optional) lag model
spatial_model_clima_lags <- inla(
  formula           = formula_bym_climate_lags,
  family            = "poisson",
  data              = dat,
  offset            = offset_pop,
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE)
)


# ------------------------------------------------------------
# Compare INLA model fits (DIC, WAIC, CPO/LCPO)
# ------------------------------------------------------------
compare_inla_models <- function(models, sort_by = c("waic", "dic", "mlcpo"), decreasing = FALSE) {
  sort_by <- match.arg(sort_by)
  
  get_safe <- function(x, path, default = NA_real_) {
    # safe nested extraction without errors
    for (p in path) {
      if (is.null(x)) return(default)
      x <- x[[p]]
    }
    if (length(x) == 0 || is.null(x)) default else x
  }
  
  res <- lapply(names(models), function(nm) {
    m <- models[[nm]]
    
    dic  <- get_safe(m, c("dic",  "dic"))
    waic <- get_safe(m, c("waic", "waic"))
    
    # CPO: lower "mlcpo" (mean log CPO) is worse; higher is better.
    # We'll compute:
    # - sum_lcpo = sum(log(CPO_i))  (often called LCPO)
    # - mlcpo    = mean(log(CPO_i))
    cpo_vec <- get_safe(m, c("cpo", "cpo"), default = NA_real_)
    if (all(is.na(cpo_vec))) {
      sum_lcpo <- NA_real_
      mlcpo    <- NA_real_
      n_bad_cpo <- NA_integer_
      p_bad_cpo <- NA_real_
    } else {
      lcpo_vec <- log(cpo_vec)
      sum_lcpo <- sum(lcpo_vec, na.rm = TRUE)
      mlcpo    <- mean(lcpo_vec, na.rm = TRUE)
      
      # some CPOs can be NA/0 -> log undefined; count as "bad"
      bad <- is.na(cpo_vec) | cpo_vec <= 0
      n_bad_cpo <- sum(bad)
      p_bad_cpo <- mean(bad)
    }
    
    data.frame(
      model      = nm,
      DIC        = as.numeric(dic),
      WAIC       = as.numeric(waic),
      sum_LCPO   = as.numeric(sum_lcpo),
      mean_LCPO  = as.numeric(mlcpo),
      n_bad_CPO  = as.integer(n_bad_cpo),
      p_bad_CPO  = as.numeric(p_bad_cpo),
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, res)
  
  # Sorting rule:
  # - DIC, WAIC: smaller is better
  # - mean_LCPO (mean log CPO): larger is better
  if (sort_by %in% c("dic", "waic")) {
    out <- out[order(out[[toupper(sort_by)]], decreasing = decreasing, na.last = TRUE), ]
  } else { # "mlcpo"
    out <- out[order(out[["mean_LCPO"]], decreasing = TRUE, na.last = TRUE), ]
  }
  
  # Add deltas vs best (by WAIC by default)
  best_waic <- suppressWarnings(min(out$WAIC, na.rm = TRUE))
  best_dic  <- suppressWarnings(min(out$DIC,  na.rm = TRUE))
  
  out$delta_WAIC <- out$WAIC - best_waic
  out$delta_DIC  <- out$DIC  - best_dic
  
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------
# Example usage: put your fitted models into a named list
# ------------------------------------------------------------
models_list <- list(
  bym_climate      = spatial_model_clima,
  bym_climate_lags = spatial_model_clima_lags
)

cmp <- compare_inla_models(models_list, sort_by = "waic")
print(cmp)


# Keep only the lag covariates you decided to retain
lag_vars_keep <- c("Max_temp_lag1","Max_temp_lag2",
                   "Min_temp_lag1","Min_temp_lag2",
                   "Prec_lag1","Prec_lag2",
                   "RH_lag1","RH_lag2")

# Convert to numeric safely (works if some are factors/characters)
to_num <- function(x) as.numeric(as.character(x))

dat <- dat %>%
  mutate(across(all_of(lag_vars_keep), to_num))

# Drop rows with missing lags (usually first 2 months for each sector)
dat <- dat %>%
  filter(if_all(all_of(lag_vars_keep), ~ !is.na(.x)))

# Standardize (recommended)
stdz <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
dat <- dat %>%
  mutate(across(all_of(lag_vars_keep), stdz))

# PC priors (adjust if you already set different ones)
# PC priors for each BYM precision component
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

fit_final$waic$waic
fit_final$dic$dic
fixef <- fit_final$summary.fixed %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  mutate(
    RR     = exp(mean),
    RR_LCL = exp(`0.025quant`),
    RR_UCL = exp(`0.975quant`)
  )

# Keep only climate lag terms
fixef_climate <- fixef %>%
  filter(term %in% lag_vars_keep) %>%
  select(term, mean, `0.025quant`, `0.975quant`, RR, RR_LCL, RR_UCL)

print(fixef_climate)

# ============================================================
# POST-PROCESSING: Spatial RR, Exceedance PP, Temporal RR + CI
# Requires: fit_final, dat (sf), INLA
# ============================================================

# ============================================================
# 1) SPATIAL: RR_i + 95% CrI + posterior exceedance P(RR>1)
# ============================================================


spat_df <- as.data.frame(fit_final$summary.random$ID)

# If an ID column already exists, keep it; otherwise take rownames
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


# Exceedance probability: PP = P(s_i > 0) = P(RR_i > 1)
spat_out$PP_RR_gt_1 <- sapply(seq_along(fit_final$marginals.random$ID), function(k) {
  1 - inla.pmarginal(q = 0, marginal = fit_final$marginals.random$ID[[k]])
})

# Join to ONE geometry per sector (avoid duplicates across months)
geo_sec <- dat %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup() %>%
  select(ID, ADM3_EN, ADM2_EN, ADM1_EN, geometry)

map_spatial <- geo_sec %>%
  left_join(spat_out, by = "ID")

# ============================================================
# 2) TEMPORAL: RR_t + 95% CrI  (FIXED)
# ============================================================

temp_out <- fit_final$summary.random$ID.time %>%
  as.data.frame() %>%
  rownames_to_column("ID.time") %>%         # <-- FIX: extract from rownames
  mutate(
    ID.time = as.integer(ID.time),
    RR_time     = exp(mean),
    RR_time_LCL = exp(`0.025quant`),
    RR_time_UCL = exp(`0.975quant`)
  ) %>%
  select(ID.time, RR_time, RR_time_LCL, RR_time_UCL)

# Add calendar labels (if present)
time_lookup <- dat %>%
  distinct(ID.time, Year, Month, Months_year) %>%
  arrange(ID.time)

temp_out <- temp_out %>%
  left_join(time_lookup, by = "ID.time") %>%
  arrange(ID.time)

# ============================================================
# 3) (Optional) SPACE–TIME interaction: RR_st + 95% CrI
# ============================================================

# ============================================================
# 5) Quick prints
# ============================================================
# Spatial Relative Risk

p_RR_cap <- ggplot(map_spatial) +
  geom_sf(aes(fill = pmin(RR_spatial, 3)), color = NA) +
  scale_fill_viridis_c(option = "C", name = "RR ",
                       breaks = c(0.5, 1, 2, 3)) +
  labs(title = "Spatial Relative Risk (RR)") +
  theme_minimal()

p_RR_cap


# Extract space–time random effects
st_df <- fit_final$summary.random$ID.space.time %>%
  as.data.frame() %>%
  rownames_to_column("ID.space.time") %>%
  mutate(ID.space.time = as.integer(ID.space.time))

# Attach index back to original data
dat_st <- dat %>%
  st_drop_geometry() %>%
  select(ID.space.time, ID, ID.time, Year)

st_join <- dat_st %>%
  left_join(st_df, by = "ID.space.time")

# Compute RR for interaction only
st_join <- st_join %>%
  mutate(
    RR_st = exp(mean)
  )
# Compute P(RR_st > 1)
st_join$PP_st_gt1 <- sapply(
  seq_along(fit_final$marginals.random$ID.space.time),
  function(k) {
    1 - inla.pmarginal(
      q = 0,
      marginal = fit_final$marginals.random$ID.space.time[[k]]
    )
  }
)


st_join2 <- st_join %>%
  mutate(
    ID = coalesce(ID.x, ID.y)   # pick whichever exists
  ) %>%
  select(-ID.x, -ID.y)

year_pp <- st_join2 %>%
  group_by(ID, Year) %>%
  summarise(
    PP_year = mean(PP_st_gt1, na.rm = TRUE),
    .groups = "drop"
  )

geo_sec <- dat %>%
  group_by(ID) %>% slice(1) %>% ungroup() %>%
  select(ID, geometry)

year_pp_sf <- geo_sec %>%
  left_join(year_pp, by = "ID") %>%
  st_as_sf()

p_PP_yearly <- ggplot(year_pp_sf) +
  geom_sf(aes(fill = PP_year), color = NA) +
  scale_fill_viridis_c(option = "C", name = "P(RR > 1)", limits = c(0, 1)) +
  facet_wrap(~ Year, ncol = 5) +
  labs(
    title = "Yearly exceedance probability maps",
    subtitle = "Posterior probability that RR > 1 (space–time interaction)"
  ) +
  theme_minimal()

p_PP_yearly

ggsave("PP_yearly_facets.png", p_PP_yearly, width = 14, height = 9, dpi = 300)


#Identify last observed month
T_last <- max(dat$ID.time)
T_last
#Build future time index (next 3 months)
# Number of forecast months
h <- 3

future_time <- data.frame(
  ID.time = (T_last + 1):(T_last + h)
)
# Build future sector–time grid

future <- expand.grid(
  ID = sort(unique(dat$ID)),
  ID.time = future_time$ID.time
)
#Add population + offset (use last observed values)
# Get last available population per sector
pop_last <- dat %>%
  group_by(ID) %>%
  filter(ID.time == T_last) %>%
  ungroup() %>%
  select(ID, Population, offset_pop)

future <- future %>%
  left_join(pop_last, by = "ID")
#Supply climate lag covariates (persistence assumption)
#We use the most recent observed lag values per sector.
lag_vars <- c("Max_temp_lag1","Max_temp_lag2",
              "Min_temp_lag1","Min_temp_lag2",
              "Prec_lag1","Prec_lag2",
              "RH_lag1","RH_lag2")

clim_last <- dat %>%
  group_by(ID) %>%
  filter(ID.time == T_last) %>%
  ungroup() %>%
  select(ID, all_of(lag_vars))

future <- future %>%
  left_join(clim_last, by = "ID")
#Create space–time index for forecasting rows
# Next available index
max_st <- max(dat$ID.space.time)

future <- future %>%
  arrange(ID, ID.time) %>%
  mutate(
    ID.space.time = max_st + row_number()
  )
#Combine observed + future data
dat_forecast <- bind_rows(
  dat %>% mutate(forecast = 0),
  future %>% mutate(
    Malaria_Cases = NA,   # missing outcome → prediction
    forecast = 1
  )
)
#Re-fit model with missing outcomes (INLA prediction trick)
fit_forecast <- inla(
  formula           = formula_final_bym,   # or formula_final_bym2 if used
  family            = "poisson",
  data              = dat_forecast,
  offset            = offset_pop,
  control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE),
  control.predictor = list(compute = TRUE, link = 1)
)
#Extract forecasts (next 3 months only)

forecast_out <- dat_forecast %>%
  mutate(
    mu_mean = fit_forecast$summary.fitted.values$mean,
    mu_lcl  = fit_forecast$summary.fitted.values$`0.025quant`,
    mu_ucl  = fit_forecast$summary.fitted.values$`0.975quant`
  ) %>%
  filter(forecast == 1) %>%
  arrange(ID, ID.time)

#Add calendar labels (Year / Month)

last_cal <- dat %>%
  filter(ID.time == T_last) %>%
  distinct(Year, Month)

forecast_out <- forecast_out %>%
  mutate(
    Month = ((last_cal$Month + ID.time - T_last - 1) %% 12) + 1,
    Year  = last_cal$Year + (last_cal$Month + ID.time - T_last - 1) %/% 12
  )
#Final forecast output
forecast_out %>%
  select(ID, Year, Month, mu_mean, mu_lcl, mu_ucl) %>%
  head()
#Save for dashboard / policy use
saveRDS(forecast_out, "forecast_malaria_next_3_months.rds")


# --- Historical: attach fitted values to dat (if not already) ---
if (is.null(dat$mu_hat)) {
  dat$mu_hat <- fit_final$summary.fitted.values$mean
  dat$mu_lcl <- fit_final$summary.fitted.values$`0.025quant`
  dat$mu_ucl <- fit_final$summary.fitted.values$`0.975quant`
}

## --- Create a month label (use Months_year if it exists) ---
hist_nat <- dat %>%
group_by(ID.time) %>%
  summarise(Months_year = dplyr::first(Months_year),
Year  = dplyr::first(Year),
Month = dplyr::first(Month),
obs   = sum(Malaria_Cases, na.rm = TRUE),
fit   = sum(mu_hat, na.rm = TRUE),
fit_lcl = sum(mu_lcl, na.rm = TRUE),
fit_ucl = sum(mu_ucl, na.rm = TRUE),
groups = "drop"
) %>%
arrange(ID.time)

# --- Forecast: aggregate to national ---
fc_nat <- forecast_out %>%
  group_by(ID.time) %>%
  summarise(
    Year  = dplyr::first(Year),
    Month = dplyr::first(Month),
    fc_mean = sum(mu_mean, na.rm = TRUE),
    fc_lcl  = sum(mu_lcl,  na.rm = TRUE),
    fc_ucl  = sum(mu_ucl,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(ID.time)

# --- Build a plotting index (continuous) ---
all_time <- sort(unique(c(hist_nat$ID.time, fc_nat$ID.time)))
time_map <- data.frame(ID.time = all_time, t_plot = seq_along(all_time))

hist_nat <- hist_nat %>% left_join(time_map, by = "ID.time")
fc_nat   <- fc_nat   %>% left_join(time_map, by = "ID.time")

# --- Plot ---
ggplot() +
  # Historical fitted uncertainty
  geom_ribbon(data = hist_nat, aes(x = t_plot, ymin = fit_lcl, ymax = fit_ucl), alpha = 0.15) +
  geom_line(data = hist_nat, aes(x = t_plot, y = fit), linewidth = 0.9) +
  
  # Historical observed
  geom_line(data = hist_nat, aes(x = t_plot, y = obs), linewidth = 0.9, linetype = "dashed") +
  
  # Forecast uncertainty + mean
  geom_ribbon(data = fc_nat, aes(x = t_plot, ymin = fc_lcl, ymax = fc_ucl), alpha = 0.20) +
  geom_line(data = fc_nat, aes(x = t_plot, y = fc_mean), linewidth = 1.1) +
  
  # Vertical separator at last observed time
  geom_vline(xintercept = max(hist_nat$t_plot), linetype = "dotted") +
  
  labs(
    x = "Time (months)",
    y = "Malaria cases (national total)",
    title = "Historical vs Fitted vs 3-month Forecast",
    subtitle = "Dashed = observed; Solid = fitted/forecast mean; Shaded = 95% credible interval"
  ) +
  theme_minimal()



# 1) Drop geometry (critical for memory)
dat0 <- sf::st_drop_geometry(dat)

# 2) Ensure fitted values exist
if (!("mu_hat" %in% names(dat0))) {
  dat0$mu_hat <- fit_final$summary.fitted.values$mean
  dat0$mu_lcl <- fit_final$summary.fitted.values$`0.025quant`
  dat0$mu_ucl <- fit_final$summary.fitted.values$`0.975quant`
}

# 3) National aggregation (much lighter now)
hist_nat <- dat0 %>%
  group_by(ID.time) %>%
  summarise(
    Months_year = dplyr::first(Months_year),
    Year  = dplyr::first(Year),
    Month = dplyr::first(Month),
    obs   = sum(Malaria_Cases, na.rm = TRUE),
    fit   = sum(mu_hat, na.rm = TRUE),
    fit_lcl = sum(mu_lcl, na.rm = TRUE),
    fit_ucl = sum(mu_ucl, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(ID.time)

# Forecast must also be non-sf (it usually is)
fc_nat <- forecast_out %>%
  group_by(ID.time) %>%
  summarise(
    Year  = dplyr::first(Year),
    Month = dplyr::first(Month),
    fc_mean = sum(mu_mean, na.rm = TRUE),
    fc_lcl  = sum(mu_lcl,  na.rm = TRUE),
    fc_ucl  = sum(mu_ucl,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(ID.time)

# Build plotting index
all_time <- sort(unique(c(hist_nat$ID.time, fc_nat$ID.time)))
time_map <- data.frame(ID.time = all_time, t_plot = seq_along(all_time))
hist_nat <- hist_nat %>% left_join(time_map, by = "ID.time")
fc_nat   <- fc_nat   %>% left_join(time_map, by = "ID.time")

# Plot
ggplot() +
  geom_ribbon(data = hist_nat, aes(x = t_plot, ymin = fit_lcl, ymax = fit_ucl), alpha = 0.15) +
  geom_line(data = hist_nat, aes(x = t_plot, y = fit), linewidth = 0.9) +
  geom_line(data = hist_nat, aes(x = t_plot, y = obs), linewidth = 0.9, linetype = "dashed") +
  geom_ribbon(data = fc_nat, aes(x = t_plot, ymin = fc_lcl, ymax = fc_ucl), alpha = 0.20) +
  geom_line(data = fc_nat, aes(x = t_plot, y = fc_mean), linewidth = 1.1) +
  geom_vline(xintercept = max(hist_nat$t_plot), linetype = "dotted") +
  labs(
    x = "Time (months)",
    y = "Malaria cases (national total)",
    title = "Historical vs Fitted vs 3-month Forecast",
    subtitle = "Dashed = observed; Solid = fitted/forecast mean; Shaded = 95% credible interval"
  ) +
  theme_minimal()


# Drop geometry to avoid memory blow-ups
dat0 <- sf::st_drop_geometry(dat) %>%
  select(ADM2_EN, ID.time, Year, Month, Months_year, Malaria_Cases)

# Attach fitted values (historical)
dat0$mu_hat <- fit_final$summary.fitted.values$mean
dat0$mu_lcl <- fit_final$summary.fitted.values$`0.025quant`
dat0$mu_ucl <- fit_final$summary.fitted.values$`0.975quant`

ggplot() +
  # -----------------------------
# Historical fitted (BLUE)
# -----------------------------
geom_ribbon(
  data = hist_nat,
  aes(x = t_plot, ymin = fit_lcl, ymax = fit_ucl),
  alpha = 0.18,
  fill = "blue"
) +
  geom_line(
    data = hist_nat,
    aes(x = t_plot, y = fit),
    linewidth = 0.9,
    color = "blue"
  ) +
  
  # Historical observed (BLUE, dashed)
  geom_line(
    data = hist_nat,
    aes(x = t_plot, y = obs),
    linewidth = 0.9,
    linetype = "dashed",
    color = "blue"
  ) +
  
  # -----------------------------
# Forecast 3 months (RED)
# -----------------------------
geom_ribbon(
  data = fc_nat,
  aes(x = t_plot, ymin = fc_lcl, ymax = fc_ucl),
  alpha = 0.22,
  fill = "red"
) +
  geom_line(
    data = fc_nat,
    aes(x = t_plot, y = fc_mean),
    linewidth = 1.1,
    color = "red"
  ) +
  
  # Separation line
  geom_vline(
    xintercept = max(hist_nat$t_plot),
    linetype = "dotted"
  ) +
  
  labs(
    x = "Time (months)",
    y = "Malaria cases (national total)",
    title = "Historical (blue) and 3-month forecast (red) malaria trends",
    subtitle = "Solid = fitted / forecast mean; Dashed = observed; Shaded = 95% credible interval"
  ) +
  theme_minimal()


######################################

# ------------------------------------------------------------
# 1) Historical (district totals): observed + fitted + CI
# ------------------------------------------------------------
dat0 <- sf::st_drop_geometry(dat) %>%
  select(ADM2_EN, ID.time, Year, Month, Months_year, Malaria_Cases)

# attach fitted values (historical)
dat0$mu_hat <- fit_final$summary.fitted.values$mean
dat0$mu_lcl <- fit_final$summary.fitted.values$`0.025quant`
dat0$mu_ucl <- fit_final$summary.fitted.values$`0.975quant`

hist_dist <- dat0 %>%
  group_by(ADM2_EN, ID.time) %>%
  summarise(
    Year  = first(Year),
    Month = first(Month),
    obs   = sum(Malaria_Cases, na.rm = TRUE),
    fit   = sum(mu_hat, na.rm = TRUE),
    fit_lcl = sum(mu_lcl, na.rm = TRUE),
    fit_ucl = sum(mu_ucl, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(ADM2_EN, ID.time)

# ------------------------------------------------------------
# 2) Forecast (district totals): mean + CI
# ------------------------------------------------------------
fc_dist <- forecast_out %>%
  sf::st_drop_geometry() %>%
  select(ADM2_EN, ID.time, Year, Month, mu_mean, mu_lcl, mu_ucl) %>%
  group_by(ADM2_EN, ID.time) %>%
  summarise(
    Year  = first(Year),
    Month = first(Month),
    fc_mean = sum(mu_mean, na.rm = TRUE),
    fc_lcl  = sum(mu_lcl,  na.rm = TRUE),
    fc_ucl  = sum(mu_ucl,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(ADM2_EN, ID.time)

# ------------------------------------------------------------
# 3) Build a plotting index (continuous within each district)
#    so forecast continues after historical
# ------------------------------------------------------------
hist_dist <- hist_dist %>%
  group_by(ADM2_EN) %>%
  mutate(t_plot = row_number()) %>%
  ungroup()

last_t <- hist_dist %>%
  group_by(ADM2_EN) %>%
  summarise(last_t = max(t_plot), .groups = "drop")

fc_dist <- fc_dist %>%
  left_join(last_t, by = "ADM2_EN") %>%
  group_by(ADM2_EN) %>%
  arrange(ID.time) %>%
  mutate(t_plot = last_t + row_number()) %>%
  ungroup()

# ------------------------------------------------------------
# 4) Plot district-by-district (facets)
# ------------------------------------------------------------
ggplot() +
  # Historical fitted (BLUE)
  geom_ribbon(
    data = hist_dist,
    aes(x = t_plot, ymin = fit_lcl, ymax = fit_ucl),
    fill = "blue", alpha = 0.18
  ) +
  geom_line(
    data = hist_dist,
    aes(x = t_plot, y = fit),
    color = "blue", linewidth = 0.8
  ) +
  geom_line(
    data = hist_dist,
    aes(x = t_plot, y = obs),
    color = "blue", linetype = "dashed", linewidth = 0.8
  ) +
  
  # Forecast (RED)
  geom_ribbon(
    data = fc_dist,
    aes(x = t_plot, ymin = fc_lcl, ymax = fc_ucl),
    fill = "red", alpha = 0.22
  ) +
  geom_line(
    data = fc_dist,
    aes(x = t_plot, y = fc_mean),
    color = "red", linewidth = 1.0
  ) +
  
  # District-specific separation line (end of historical)
  geom_vline(
    data = last_t,
    aes(xintercept = last_t),
    linetype = "dotted"
  ) +
  
  facet_wrap(~ ADM2_EN, scales = "free_y", ncol = 5) +
  labs(
    x = "Time (months)",
    y = "Malaria cases (district total)",
    title = "District-level historical (blue) and 3-month forecast (red) malaria trends",
    subtitle = "Blue dashed = observed; Blue solid = fitted; Red = forecast; Shaded = 95% credible interval"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)
  )


