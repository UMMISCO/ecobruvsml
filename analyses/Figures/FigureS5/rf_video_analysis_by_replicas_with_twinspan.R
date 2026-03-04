# =============================================================================
# Script name: rf_video_analysis_by_replicas_with_twinspan.R
# Author: Estephe Kana, Eugeni Belda, Edi Prifti
# Date created: 2025-01-01
# Purpose: Run Random Forest classifiers using IndVal and TWINSPAN indicator
#          species on replica-split data (ABORE train / MBERE test).
# Inputs:  analyses/analysis_outputs/replicas_analysis_results/ (split .Rda files)
# Outputs: analyses/analysis_outputs/rf_analysis_by_replicas/rf_analysis_results_split_by_replicas_with_twinspan.Rda
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script:
# cd /home/estephe/Documents/Papers/2026-aime_key_fish_bruvs/Bruvs-analysis-code/analyses/Figures/FigureS5
# Rscript rf_video_analysis_by_replicas_with_twinspan.R
# -----------------------------------------------------------------------------


# Detect project root automatically (works for any user/machine after cloning)
if (!require("here", quietly = TRUE)) install.packages("here")
base_dir <- here::here()

if (!require("randomForest")) install.packages("randomForest"); library(randomForest)
if (!require("caret"))        install.packages("caret");        library(caret)
if (!require("readr"))        install.packages("readr");        library(readr)
if (!require("ROCR"))         install.packages("ROCR");         library(ROCR)
if (!require("stringr"))      install.packages("stringr");      library(stringr)
if (!require("ggpubr"))       install.packages("ggpubr");       library(ggpubr)
if (!require("purrr"))        install.packages("purrr");        library(purrr)
if (!require("dplyr"))        install.packages("dplyr");        library(dplyr)

# ── Shared RF helper ──────────────────────────────────────────────────────────
# Trains a Random Forest on X.train/y.train using a given set of species,
# evaluates on X.test/y.test, and returns a metrics data frame.
run_rf <- function(X.train, y.train, X.test, y.test,
                   species.names, comp, source, method_label) {
  
  # Subset training data to indicator species
  rf.train_data.X <- X.train[, species.names, drop = FALSE]
  rf.train_data   <- as.data.frame(rf.train_data.X)
  rf.train_data$class <- as.factor(y.train)
  # Remove spaces in species names (RF formula parser requirement)
  colnames(rf.train_data) <- gsub(" ", ".", colnames(rf.train_data))
  
  set.seed(20)
  rf.model <- randomForest(class ~ ., data = rf.train_data, proximity = TRUE)
  
  # Add any training columns missing from test set (set to 0)
  missing_cols <- setdiff(colnames(rf.train_data.X), colnames(X.test))
  for (col in missing_cols) X.test[[col]] <- 0
  
  # Align test set to training columns
  rf.test_data <- X.test[, colnames(rf.train_data.X), drop = FALSE]
  rf.test_data$class <- as.factor(y.test)
  colnames(rf.test_data) <- gsub(" ", ".", colnames(rf.test_data))
  
  # Class predictions
  pred_class <- predict(rf.model, newdata = rf.test_data, type = "class")
  cm         <- confusionMatrix(table(pred_class, rf.test_data$class))
  
  acc <- round(cm$overall["Accuracy"], digits = 3)
  f1  <- round(cm$byClass["F1"],       digits = 3)
  
  # Probability predictions for AUC
  pred_prob       <- predict(rf.model, newdata = rf.test_data, type = "prob")
  pred_obj        <- prediction(pred_prob[, 2], rf.test_data$class)
  auc             <- performance(pred_obj, measure = "auc")@y.values[[1]]
  auc             <- round(auc, digits = 3)
  
  data.frame(
    metrics    = c("auc", "accuracy", "f1"),
    value      = c(auc, acc, f1),
    comparison = comp,
    source     = source,
    method     = method_label,
    featNum    = ncol(rf.train_data.X),
    stringsAsFactors = FALSE
  )
}

# ── RF on Indval results ──────────────────────────────────────────────────────
rf_indval_analysis <- function(data, source, comp) {
  
  data_source <- if (source == "maxN") {
    data$indvalout.maxn[[comp]]
  } else if (source == "presAbs") {
    data$indvalout.bin[[comp]]
  } else {
    stop("Invalid source specified. Use 'maxN' or 'presAbs'.")
  }
  
  indval_species <- data_source$indval_pval[data_source$indval_pval$IsIndSp == 1, ]
  species.names  <- as.character(indval_species$feature)
  
  run_rf(
    X.train       = data_source$Comp_data$X_train,
    y.train       = data_source$Comp_data$y_train,
    X.test        = data_source$Comp_data$X_test,
    y.test        = data_source$Comp_data$y_test,
    species.names = species.names,
    comp          = comp,
    source        = source,
    method_label  = "RF_Indval"
  )
}

# ── RF on Predomics results ───────────────────────────────────────────────────
rf_pred_analysis <- function(data, source, comp, model) {
  
  data_source <- if (source == "maxN") {
    data$predout.maxn[[comp]]
  } else if (source == "presAbs") {
    data$predout.bin[[comp]]
  } else {
    stop("Invalid source specified. Use 'maxN' or 'presAbs'.")
  }
  
  pred_species  <- data_source$pred_out_fbm[data_source$pred_out_fbm$IsIndSp == 1, ]
  species.names <- as.character(pred_species$feature)
  
  run_rf(
    X.train       = data_source$Comp_data$X_train,
    y.train       = data_source$Comp_data$y_train,
    X.test        = data_source$Comp_data$X_test,
    y.test        = data_source$Comp_data$y_test,
    species.names = species.names,
    comp          = comp,
    source        = source,
    method_label  = paste("RF", model, sep = "_")
  )
}

# ── RF on Twinspan results ────────────────────────────────────────────────────
rf_twinspan_analysis <- function(data, source, comp) {
  
  # Twinspan stores results under twinspan.maxn / twinspan.bin (named by comparison)
  data_source <- if (source == "maxN") {
    data$twinspan.maxn[[comp]]
  } else if (source == "presAbs") {
    data$twinspan.bin[[comp]]
  } else {
    stop("Invalid source specified. Use 'maxN' or 'presAbs'.")
  }
  
  if (is.null(data_source)) {
    warning("No Twinspan data found for comp='", comp, "' source='", source, "'. Skipping.")
    return(NULL)
  }
  
  # Indicator species are flagged in indicator.df (IsIndSp == 1)
  tw_species    <- data_source$indicator.df[data_source$indicator.df$IsIndSp == 1, ]
  species.names <- as.character(tw_species$feature)
  species.names <- species.names[!is.na(species.names)]
  
  if (length(species.names) == 0) {
    warning("No indicator species found for comp='", comp, "' source='", source, "'. Skipping.")
    return(NULL)
  }
  
  # Twinspan Comp_data uses the same structure as Indval/Predomics
  run_rf(
    X.train       = data_source$Comp_data$X_train,
    y.train       = data_source$Comp_data$y_train,
    X.test        = data_source$Comp_data$X_test,
    y.test        = data_source$Comp_data$y_test,
    species.names = species.names,
    comp          = comp,
    source        = source,
    method_label  = "RF_Twinspan"
  )
}

# ── Load data ─────────────────────────────────────────────────────────────────
load(file.path(base_dir, "analyses", "Figures", "FigureS5", "results_replicas_analysis", "all_analysis_results_split_replicas_with_twinspan.Rda"))

# ── Parameter grid ────────────────────────────────────────────────────────────
param_grid <- expand.grid(
  replica_comp = names(alldatalist.video_split_replicas),
  prev_rate    = c("prev_0", "prev_10"),
  model        = c("Indval", "bininter", "terinter", "Twinspan"),
  comp         = c("BAR_BAY", "BAR_LAG", "BAY_LAG", "INSHORE_OFFSHORE"),
  stringsAsFactors = FALSE
)

# ── Dispatch table: model name → analysis function ────────────────────────────
analysis_func <- list(
  Indval   = rf_indval_analysis,
  bininter = rf_pred_analysis,
  terinter = rf_pred_analysis,
  Twinspan = rf_twinspan_analysis
)

# ── Main loop ─────────────────────────────────────────────────────────────────
rf.metrics.list <- list()

process_entry <- function(replica_comp, prev_rate, model, comp) {
  
  data_model <- alldatalist.video_split_replicas[[replica_comp]][[prev_rate]][[model]]
  
  if (is.null(data_model)) {
    message("Skipping: replica=", replica_comp, " prev=", prev_rate,
            " model=", model, " (no data loaded)")
    return(invisible(NULL))
  }
  
  analysis_fn <- analysis_func[[model]]
  
  for (source in c("maxN", "presAbs")) {
    if(model=="bininter" | model=="terinter"){
      result <- tryCatch(
        analysis_fn(data_model, source = source, comp = comp, model=model),
        error = function(e) {
          message("Error in ", model, " / ", comp, " / ", source, ": ", conditionMessage(e))
          NULL
        }
      )
    } else{
      result <- tryCatch(
        analysis_fn(data_model, source = source, comp = comp),
        error = function(e) {
          message("Error in ", model, " / ", comp, " / ", source, ": ", conditionMessage(e))
          NULL
        }
      )
    }
    
    if (is.null(result)) next
    
    result$prev_rate    <- prev_rate
    result$replica_comp <- replica_comp
    
    rf.metrics.list[[replica_comp]][[comp]][[source]][[model]][[prev_rate]] <<- result
  }
}

pwalk(param_grid, process_entry)

# ── Flatten results ───────────────────────────────────────────────────────────
flatten_and_bind <- function(nested_list) {
  if (all(sapply(nested_list, is.data.frame))) {
    do.call(rbind, nested_list)
  } else {
    do.call(rbind, lapply(nested_list, flatten_and_bind))
  }
}

rf.metrics.df <- flatten_and_bind(rf.metrics.list)
rownames(rf.metrics.df) <- NULL

# ── Save ──────────────────────────────────────────────────────────────────────
save(rf.metrics.df,
     file = file.path(base_dir, "analyses", "Figures", "FigureS5", "rf_analysis_by_replicas", "rf_analysis_results_split_by_replicas_with_twinspan.Rda"))

message("Done. RF results saved.")