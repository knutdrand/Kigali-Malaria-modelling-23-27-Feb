# predict.R — CHAP-compatible predict entry point
# BYM + RW1 + IID space-time interaction model with lagged climate covariates
#
# Usage:
#   Rscript predict.R {model} {historic_data} {future_data} {out_file} {model_config} {polygons}
#
# Column mapping (adapters in MLproject):
#   Cases    = disease_cases
#   E        = population
#   month    = month
#   ID_year  = year
#   ID_spat  = location

library(yaml)
library(INLA)
library(sf)
library(spdep)
library(dplyr)
library(tidyr)
library(stringr)

# ---------------------------------------------------------------------------
# Config parsing
# ---------------------------------------------------------------------------
parse_model_configuration <- function(file_path) {
  config <- yaml::yaml.load_file(file_path)

  user_option_values <- if (!is.null(config$user_option_values)) {
    config$user_option_values
  } else {
    list()
  }

  additional_continuous_covariates <- if (!is.null(config$additional_continuous_covariates)) {
    config$additional_continuous_covariates
  } else {
    character()
  }

  list(
    user_option_values = user_option_values,
    additional_continuous_covariates = additional_continuous_covariates
  )
}

# ---------------------------------------------------------------------------
# Main predict function
# ---------------------------------------------------------------------------
predict_chap <- function(model_fn, hist_fn, future_fn, preds_fn, config_fn, polygons_fn) {

  # --- Parse config ---
  covariate_names <- character()
  if (config_fn != "" && file.exists(config_fn)) {
    cat("Loading model configuration from:", config_fn, "\n")
    config <- parse_model_configuration(config_fn)
    covariate_names <- config$additional_continuous_covariates
  }
  cat("Climate covariates:", if (length(covariate_names) == 0) "(none)" else paste(covariate_names, collapse = ", "), "\n")

  # --- Read data ---
  cat("Reading historic data from:", hist_fn, "\n")
  historic_df <- read.csv(hist_fn, stringsAsFactors = FALSE)

  cat("Reading future data from:", future_fn, "\n")
  future_df <- read.csv(future_fn, stringsAsFactors = FALSE)
  future_df$Cases <- NA_integer_

  df <- bind_rows(historic_df, future_df)

  # --- Read polygons and build adjacency ---
  cat("Reading polygons from:", polygons_fn, "\n")
  polygons <- st_read(polygons_fn, quiet = TRUE) %>% st_make_valid()

  # Determine location identifier field in polygons:
  # CHAP convention uses "id", Rwanda shapefiles use "ADM3_EN"
  poly_names <- names(polygons)
  if ("id" %in% poly_names) {
    loc_field <- "id"
  } else if ("ADM3_EN" %in% poly_names) {
    loc_field <- "ADM3_EN"
    # Fix known Rwanda sector name mismatches
    polygons <- polygons %>%
      mutate(
        ADM3_EN = case_when(
          ADM3_EN == "Mageregere" ~ "Mageragere",
          ADM3_EN == "Shyrongi"   ~ "Shyorongi",
          ADM3_EN == "Ririma"     ~ "Rilima",
          TRUE                    ~ ADM3_EN
        )
      )
  } else {
    stop("Cannot find location identifier in polygons. Expected 'id' or 'ADM3_EN'.", call. = FALSE)
  }
  cat("Using polygon location field:", loc_field, "\n")

  polygons <- polygons %>% mutate(Sec_ID = row_number())

  # Build adjacency graph
  graph_file <- tempfile(fileext = ".graph")
  nb_obj <- poly2nb(polygons, queen = FALSE)
  nb2INLA(graph_file, nb_obj)
  Rwa_adj <- normalizePath(graph_file)
  cat("Adjacency graph built. Sectors:", nrow(polygons), "\n")

  # --- Build location → Sec_ID mapping ---
  loc_map <- polygons %>%
    st_drop_geometry() %>%
    select(loc_id = all_of(loc_field), Sec_ID) %>%
    mutate(loc_id_lower = str_to_lower(str_squish(as.character(loc_id))))

  df <- df %>%
    mutate(ID_spat_lower = str_to_lower(str_squish(as.character(ID_spat))))

  df <- df %>%
    left_join(loc_map, by = c("ID_spat_lower" = "loc_id_lower"))

  n_matched <- sum(!is.na(df$Sec_ID))
  cat("Rows matched to polygons:", n_matched, "of", nrow(df), "\n")
  if (n_matched == 0) {
    stop("No rows matched to polygon locations. Check ID_spat values vs polygon '", loc_field, "' field.", call. = FALSE)
  }

  # Drop unmatched rows
  df <- df %>% filter(!is.na(Sec_ID))

  # --- Build indices ---
  df <- df %>%
    arrange(Sec_ID, ID_year, month) %>%
    mutate(
      ID      = Sec_ID,
      ID.time = as.integer(factor(
        paste(ID_year, sprintf("%02d", month)),
        levels = unique(paste(ID_year, sprintf("%02d", month)))
      )),
      ID.space.time = row_number(),
      E = as.numeric(E)
    )

  # --- Create lagged climate covariates ---
  lag_vars <- character()

  if (length(covariate_names) > 0) {
    for (cov in covariate_names) {
      if (!(cov %in% names(df))) {
        cat("Warning: covariate", cov, "not found in data, skipping.\n")
        next
      }
      lag1_name <- paste0(cov, "_lag1")
      lag2_name <- paste0(cov, "_lag2")

      df <- df %>%
        arrange(Sec_ID, ID_year, month) %>%
        group_by(Sec_ID) %>%
        mutate(
          !!lag1_name := lag(as.numeric(.data[[cov]]), 1),
          !!lag2_name := lag(as.numeric(.data[[cov]]), 2)
        ) %>%
        ungroup()

      lag_vars <- c(lag_vars, lag1_name, lag2_name)
    }

    # Drop rows with missing lags (first 2 months per sector in historic data)
    n_before <- nrow(df)
    df <- df %>%
      filter(if_all(all_of(lag_vars), ~ !is.na(.x)))
    cat("Rows dropped due to missing lags:", n_before - nrow(df), "\n")

    # Re-index after filtering
    df <- df %>%
      arrange(Sec_ID, ID_year, month) %>%
      mutate(
        ID.time = as.integer(factor(
          paste(ID_year, sprintf("%02d", month)),
          levels = unique(paste(ID_year, sprintf("%02d", month)))
        )),
        ID.space.time = row_number()
      )

    # Standardize lag variables
    stdz <- function(x) {
      s <- sd(x, na.rm = TRUE)
      if (is.na(s) || s == 0) return(x - mean(x, na.rm = TRUE))
      (x - mean(x, na.rm = TRUE)) / s
    }
    df <- df %>%
      mutate(across(all_of(lag_vars), stdz))
  }

  cat("Modelling dataset:", nrow(df), "rows,",
      n_distinct(df$ID), "sectors,",
      max(df$ID.time), "time points\n")

  # --- Build formula ---
  hyper_bym <- list(
    prec.unstruct = list(prior = "pc.prec", param = c(1, 0.01)),
    prec.spatial  = list(prior = "pc.prec", param = c(1, 0.01))
  )

  # Base formula string
  formula_str <- "Cases ~ 1"

  # Add lagged covariates
  if (length(lag_vars) > 0) {
    formula_str <- paste(formula_str, "+", paste(lag_vars, collapse = " + "))
  }

  # The random effects are added programmatically
  model_formula <- as.formula(formula_str)

  # Add random effects using update()
  model_formula <- update(
    model_formula,
    . ~ . +
      f(ID, model = "bym", graph = Rwa_adj, scale.model = TRUE,
        constr = TRUE, hyper = hyper_bym) +
      f(ID.time, model = "rw1", constr = TRUE, scale.model = TRUE,
        hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
      f(ID.space.time, model = "iid", constr = TRUE,
        hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
  )

  cat("Formula:", deparse(model_formula, width.cutoff = 200), "\n")

  # --- Fit model ---
  cat("Fitting INLA model...\n")
  model <- inla(
    formula           = model_formula,
    family            = "poisson",
    data              = df,
    offset            = log(E),
    control.compute   = list(dic = TRUE, config = TRUE, cpo = TRUE, return.marginals = FALSE),
    control.predictor = list(compute = TRUE, link = 1),
    verbose           = FALSE,
    safe              = FALSE
  )

  cat("Model fitted. DIC:", model$dic$dic, "\n")

  # --- Posterior sampling ---
  cat("Drawing posterior samples...\n")
  idx.pred <- which(is.na(df$Cases))
  mpred <- length(idx.pred)
  cat("Prediction rows:", mpred, "\n")

  s <- 1000
  y.pred <- matrix(NA, mpred, s)

  xx <- inla.posterior.sample(s, model)
  xx.s <- inla.posterior.sample.eval(
    function(idx.pred) Predictor[idx.pred],
    xx,
    idx.pred = idx.pred
  )

  # Draw Poisson counts from each sample's linear predictor
  for (s.idx in 1:s) {
    y.pred[, s.idx] <- rpois(mpred, lambda = exp(xx.s[, s.idx]))
  }

  # --- Build output ---
  new.df <- data.frame(
    time_period = df$time_period[idx.pred],
    location    = df$location[idx.pred],
    y.pred
  )
  colnames(new.df) <- c("time_period", "location", paste0("sample_", 0:(s - 1)))

  # Write predictions
  cat("Writing predictions to:", preds_fn, "\n")
  write.csv(new.df, preds_fn, row.names = FALSE)

  # Save model
  cat("Saving model to:", model_fn, "\n")
  saveRDS(model, file = model_fn)

  cat("Done.\n")
}

# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 4) {
  cat("Running predictions\n")
  print(args)

  model_fn    <- args[1]
  hist_fn     <- args[2]
  future_fn   <- args[3]
  preds_fn    <- args[4]
  config_fn   <- if (length(args) >= 5) args[5] else ""
  polygons_fn <- if (length(args) >= 6) args[6] else ""

  if (polygons_fn == "") {
    stop("polygons path is required for the BYM spatial model", call. = FALSE)
  }

  predict_chap(model_fn, hist_fn, future_fn, preds_fn, config_fn, polygons_fn)
}
