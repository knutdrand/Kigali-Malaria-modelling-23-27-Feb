# train.R — No-op for CHAP compatibility
# The BYM model fits during predict (INLA missing-data trick),
# so no separate training step is needed.
#
# Column mapping (adapters in MLproject):
#   Cases    = disease_cases
#   E        = population
#   month    = month
#   ID_year  = year
#   ID_spat  = location

library(INLA)

train_chap <- function(train_fn, model_fn, config_fn, polygons_fn) {
  # No-op: model is fitted during prediction
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
  train_fn    <- args[1]
  model_fn    <- args[2]
  config_fn   <- if (length(args) >= 3) args[3] else ""
  polygons_fn <- if (length(args) >= 4) args[4] else ""

  train_chap(train_fn, model_fn, config_fn, polygons_fn)
}
