# =============================================================================
# Script name: rf_analysis_by_site_80_20_with_twinspan.R
# Author: Estephe Kana, Eugeni Belda, Edi Prifti
# Date created: 2025-01-01
# Purpose: Evaluate Random Forest classifiers using IndVal and TWINSPAN indicator species, on 80/20 train/test split data.
# Inputs:  analyses/analysis_outputs/ (loaded .Rda files with twinspan)
# Outputs: analyses/analysis_outputs/rf_analysis_results_with_twinspan_80_20.Rda
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script:
# cd /home/estephe/Documents/Papers/2026-aime_key_fish_bruvs/Bruvs-analysis-code/analyses/scripts
# Rscript rf_analysis_by_site_80_20_with_twinspan.R
# -----------------------------------------------------------------------------

# RF analysis
# This script is used to identify the set of key species that characterize sites (Barrier, Bay, Lagoon) with efficiency 

if(!require("randomForest")) install.packages("randomForest") ; library(randomForest)
if(!require("caret")) install.packages("caret") ; library(caret)
if (!require("readr")) install.packages("readr"); library(readr)
if (!require("ROCR")) install.packages("ROCR"); library(ROCR)
if (!require("stringr")) install.packages("stringr"); library(stringr)
if (!require("ggpubr")) install.packages("ggpubr"); library(ggpubr)


# RF analysis on indval results
rf_indval_analysis <-function (data, source, comp){
  
  if(source== 'maxN'){
    # get abundance results data from indval
    data_source= data$indvalout.maxn[[comp]]
  } else if(source=='presAbs'){
    # get abundance results data from indval
    data_source= data$indvalout.bin[[comp]]
  }else {
    stop("Invalid source specified. Use 'maxN' or 'presAbs'.")
  }
  
  if (is.null(data_source)) {
    warning(paste("No data found for comparison:", comp, "source:", source))
    return(NULL)
  }
  
  # get input data
  X.train= data_source$Comp_data$X_train
  y.train= data_source$Comp_data$y_train
  X.test= data_source$Comp_data$X_test
  y.test= data_source$Comp_data$y_test
  
  # get indicator species 
  indval_species= data_source$indval_pval[data_source$indval_pval$IsIndSp==1, ]
  #Get the species in the indval results on train dataset
  indval_species.names <- as.character(indval_species$feature)
  
  if (length(indval_species.names) == 0) {
    warning(paste("No indicator species found for comparison:", comp, "source:", source))
    return(NULL)
  }
  
  if (length(indval_species.names) == 1) {
    warning(paste("Only one indicator species found for comparison:", comp, "source:", source, "- RF requires at least 2 features. Skipping."))
    return(NULL)
  }
  #Subset the train data to indval_species.names
  rf.train_data.X <- X.train[,indval_species.names]
  # bind samples and targets in a dataframe
  rf.train_data <- as.data.frame(rf.train_data.X)
  #Add class variable to the X table
  rf.train_data$class <- as.factor(y.train)
  ##need to remove spaces in colnames (species names) before running the rf model (if not error)
  colnames(rf.train_data) <- gsub(" ",".", colnames(rf.train_data))
  colnames(rf.train_data.X) <- gsub(" ",".", colnames(rf.train_data.X))
  
  set.seed(20)
  rf.train.model <- randomForest(class~., data = rf.train_data, proximity=TRUE) 
  
  # Identify columns in rf.train_data.X that are missing in X.test
  missing_columns <- setdiff(colnames(rf.train_data.X), colnames(X.test))
  # Add missing columns to X.test with values set to zero
  
  for (col in missing_columns) {
    X.test[[col]] <- 0
  }
  
  colnames(X.test) <- gsub(" ",".", colnames(X.test))
  
  # Ensure rf.test_data has the same columns as rf.train_data.X, in the same order
  rf.test_data <- X.test[, colnames(rf.train_data.X), drop = FALSE]
  #Add class variable to the X table
  rf.test_data$class <- as.factor(y.test)
  ##need to remove spaces in colnames (species names) before running the rf model (if not error)
  colnames(rf.test_data) <- gsub(" ",".", colnames(rf.test_data))
  
  # make prediction on test data
  indval_test <- predict(rf.train.model, newdata = rf.test_data, type= "class")
  
  # The prediction to compute the confusion matrix and see the accuracy score 
  indval.test.cm=confusionMatrix(table(indval_test,rf.test_data$class)) 
  
  # get acc, auc and f1 metrics of rf on indval species using test data
  indval.pred_test.acc <- round(indval.test.cm$overall['Accuracy'], digits = 3)
  indval.pred_test.f1 <- round(indval.test.cm$byClass['F1'], digits = 3)
  
  # get the auc of rf model on test data
  indval_rf.test.pred=predict(rf.train.model, newdata = rf.test_data, type = "prob")
  indval_rf.test.prediction = prediction(indval_rf.test.pred[,2], rf.test_data$class)
  indval.pred_test.auc= ROCR::performance(indval_rf.test.prediction, measure = "auc")@y.values[[1]] 
  indval.pred_test.auc= round(indval.pred_test.auc, digits = 3)
  
  #get metric performances of rf with indval species in testing data
  rf.results <- data.frame(metrics= c("auc","accuracy","f1"),
                           value= c(indval.pred_test.auc, indval.pred_test.acc, indval.pred_test.f1),
                           comparison=comp,
                           source= source,
                           featNum= ncol(rf.train_data.X))
  return(rf.results)
}

# RF analysis on predomics results
rf_pred_analysis <-function (data, source, comp){
  
  if(source== 'maxN'){
    # get abundance results data from indval
    data_source= data$predout.maxn[[comp]]
  } else if(source=='presAbs'){
    # get abundance results data from indval
    data_source= data$predout.bin[[comp]]
  }else {
    stop("Invalid  source specified. Use 'maxN' or 'presAbs'.")
  }
  
  if (is.null(data_source)) {
    warning(paste("No data found for comparison:", comp, "source:", source))
    return(NULL)
  }
  
  # get input data
  X.train= data_source$Comp_data$X_train
  y.train= data_source$Comp_data$y_train
  X.test= data_source$Comp_data$X_test
  y.test= data_source$Comp_data$y_test
  
  # get indicator species 
  pred_species= data_source$pred_out_fbm[data_source$pred_out_fbm$IsIndSp==1, ]
  #Get the species in the indval results on train dataset
  pred_species.names <- as.character(pred_species$feature)
  
  if (length(pred_species.names) == 0) {
    warning(paste("No indicator species found for comparison:", comp, "source:", source))
    return(NULL)
  }
  
  if (length(pred_species.names) == 1) {
    warning(paste("Only one indicator species found for comparison:", comp, "source:", source, "- RF requires at least 2 features. Skipping."))
    return(NULL)
  }
  #Subset the train data to pred_species.names
  rf.train_data.X <- X.train[, pred_species.names]
  # bind samples and targets in a dataframe
  rf.train_data <- as.data.frame(rf.train_data.X)
  #Add class variable to the X table
  rf.train_data$class <- as.factor(y.train)
  ##need to remove spaces in colnames (species names) before running the rf model (if not error)
  colnames(rf.train_data) <- gsub(" ",".", colnames(rf.train_data))
  colnames(rf.train_data.X) <- gsub(" ",".", colnames(rf.train_data.X))
  
  set.seed(20)
  rf.train.model <- randomForest(class~., data = rf.train_data, proximity=TRUE) 
  
  # Identify columns in rf.train_data.X that are missing in X.test
  missing_columns <- setdiff(colnames(rf.train_data.X), colnames(X.test))
  # Add missing columns to X.test with values set to zero
  
  for (col in missing_columns) {
    X.test[[col]] <- 0
  }
  
  colnames(X.test) <- gsub(" ",".", colnames(X.test))
  
  # Ensure rf.test_data has the same columns as rf.train_data.X, in the same order
  rf.test_data <- X.test[, colnames(rf.train_data.X), drop = FALSE]
  #Add class variable to the X table
  rf.test_data$class <- as.factor(y.test)
  ##need to remove spaces in colnames (species names) before running the rf model (if not error)
  colnames(rf.test_data) <- gsub(" ",".", colnames(rf.test_data))
  
  # make prediction on test data
  pred_test <- predict(rf.train.model, newdata = rf.test_data, type= "class")
  
  # The prediction to compute the confusion matrix and see the accuracy score 
  pred.test.cm=confusionMatrix(table(pred_test,rf.test_data$class)) 
  
  # get acc, auc and f1 metrics of rf on predomics species using test data
  pred_test.acc <- round(pred.test.cm$overall['Accuracy'], digits = 3)
  pred_test.f1 <- round(pred.test.cm$byClass['F1'], digits = 3)
  
  # get the auc of rf model on test data
  pred_rf.test=predict(rf.train.model, newdata = rf.test_data, type = "prob")
  pred_rf.test.prediction = prediction(pred_rf.test[,2], rf.test_data$class)
  pred_test.auc= ROCR::performance(pred_rf.test.prediction, measure = "auc")@y.values[[1]] 
  pred_test.auc= round(pred_test.auc, digits = 3)
  
  #get metric performances of rf with indval species in testing data
  rf.results <- data.frame(metrics= c("auc","accuracy","f1"),
                           value= c(pred_test.auc, pred_test.acc, pred_test.f1),
                           comparison=comp,
                           source= source,
                           featNum= ncol(rf.train_data.X))
  return(rf.results)
}

# RF analysis on twinspan results
rf_twinspan_analysis <-function (data, source, comp){
  
  if(source== 'maxN'){
    # get abundance results data from twinspan
    data_source= data$twinspan.maxn[[comp]]
  } else if(source=='presAbs'){
    # get abundance results data from twinspan
    data_source= data$twinspan.bin[[comp]]
  }else {
    stop("Invalid source specified. Use 'maxN' or 'presAbs'.")
  }
  
  if (is.null(data_source)) {
    warning(paste("No data found for comparison:", comp, "source:", source))
    return(NULL)
  }
  
  # get input data
  X.train= data_source$Comp_data$X_train
  y.train= data_source$Comp_data$y_train
  X.test= data_source$Comp_data$X_test
  y.test= data_source$Comp_data$y_test
  
  # get indicator species 
  twinspan_species= data_source$indicator.df[data_source$indicator.df$IsIndSp==1, ]
  #Get the species in the twinspan results on train dataset
  twinspan_species.names <- as.character(twinspan_species$feature)
  
  if (length(twinspan_species.names) == 0) {
    warning(paste("No indicator species found for comparison:", comp, "source:", source))
    return(NULL)
  }
  
  if (length(twinspan_species.names) == 1) {
    warning(paste("Only one indicator species found for comparison:", comp, "source:", source, "- RF requires at least 2 features. Skipping."))
    return(NULL)
  }
  
  #Subset the train data to twinspan_species.names
  rf.train_data.X <- X.train[,twinspan_species.names]
  # bind samples and targets in a dataframe
  rf.train_data <- as.data.frame(rf.train_data.X)
  #Add class variable to the X table
  rf.train_data$class <- as.factor(y.train)
  ##need to remove spaces in colnames (species names) before running the rf model (if not error)
  colnames(rf.train_data) <- gsub(" ",".", colnames(rf.train_data))
  colnames(rf.train_data.X) <- gsub(" ",".", colnames(rf.train_data.X))
  
  set.seed(20)
  rf.train.model <- randomForest(class~., data = rf.train_data, proximity=TRUE) 
  
  # Identify columns in rf.train_data.X that are missing in X.test
  missing_columns <- setdiff(colnames(rf.train_data.X), colnames(X.test))
  # Add missing columns to X.test with values set to zero
  
  for (col in missing_columns) {
    X.test[[col]] <- 0
  }
  
  colnames(X.test) <- gsub(" ",".", colnames(X.test))
  
  # Ensure rf.test_data has the same columns as rf.train_data.X, in the same order
  rf.test_data <- X.test[, colnames(rf.train_data.X), drop = FALSE]
  #Add class variable to the X table
  rf.test_data$class <- as.factor(y.test)
  ##need to remove spaces in colnames (species names) before running the rf model (if not error)
  colnames(rf.test_data) <- gsub(" ",".", colnames(rf.test_data))
  
  # make prediction on test data
  twinspan_test <- predict(rf.train.model, newdata = rf.test_data, type= "class")
  
  # The prediction to compute the confusion matrix and see the accuracy score 
  twinspan.test.cm=confusionMatrix(table(twinspan_test,rf.test_data$class)) 
  
  # get acc, auc and f1 metrics of rf on twinspan species using test data
  twinspan.pred_test.acc <- round(twinspan.test.cm$overall['Accuracy'], digits = 3)
  twinspan.pred_test.f1 <- round(twinspan.test.cm$byClass['F1'], digits = 3)
  
  # get the auc of rf model on test data
  twinspan_rf.test.pred=predict(rf.train.model, newdata = rf.test_data, type = "prob")
  twinspan_rf.test.prediction = prediction(twinspan_rf.test.pred[,2], rf.test_data$class)
  twinspan.pred_test.auc= ROCR::performance(twinspan_rf.test.prediction, measure = "auc")@y.values[[1]] 
  twinspan.pred_test.auc= round(twinspan.pred_test.auc, digits = 3)
  
  #get metric performances of rf with twinspan species in testing data
  rf.results <- data.frame(metrics= c("auc","accuracy","f1"),
                           value= c(twinspan.pred_test.auc, twinspan.pred_test.acc, twinspan.pred_test.f1),
                           comparison=comp,
                           source= source,
                           featNum= ncol(rf.train_data.X))
  return(rf.results)
}


# Detect project root automatically (works for any user/machine after cloning)
if (!require("here", quietly = TRUE)) install.packages("here")
base_dir <- here::here()

root_path <- file.path(base_dir, "analyses", "analysis_outputs")

# Get a list of all directories within the root_path recursively
all_dirs <- list.dirs(path = root_path, full.names = TRUE, recursive = TRUE)

# Filter directories that contain "output_data" in their names
output_data_dirs <- grep("output_data/analysis_80_20", all_dirs, value = TRUE)

files <- list.files(output_data_dirs, recursive = TRUE)

alldatalist.video_80_20 <- list()

for (f in files) {
  matches          <- str_match(f, "^([^_]+).*_data_([^.]+)")
  method           <- matches[, 2]
  prev_seed        <- matches[, 3]
  prevalence_rate  <- sub("^(prev_\\d+)_.*", "\\1", prev_seed)
  seed             <- sub(".*_(seed_\\d+)$", "\\1", prev_seed)
  file_dir_path    <- paste(root_path, paste(method, "output_data", sep = "_"), "analysis_80_20", sep = "/")
  fobj             <- load(paste(file_dir_path, f, sep = "/"))
  fobj             <- mget(fobj)
  alldatalist.video_80_20[[seed]][[prevalence_rate]][[method]] <- fobj
}

# save Indval/Predomics/twinspan analysis results data
save(alldatalist.video_80_20, file = paste(root_path, "all_analysis_results_with_twinspan_80_20.Rda", sep="/"))


# Map each model to its analysis function
model_fun <- list(
  Indval  = rf_indval_analysis,
  bininter = rf_pred_analysis,
  terinter = rf_pred_analysis,
  Twinspan = rf_twinspan_analysis
)

# Build all combinations
combinations <- expand.grid(
  seed      = names(alldatalist.video_80_20),
  prev_rate = c("prev_0", "prev_10"),
  model     = c("Indval", "bininter", "terinter", "Twinspan"),
  comp      = c("BAR_BAY", "BAR_LAG", "BAY_LAG", "INSHORE_OFFSHORE"),
  source    = c("maxN", "presAbs"),
  stringsAsFactors = FALSE
)

# Single function to run one combination
run_rf <- function(seed, prev_rate, model, comp, source) {
  
  results_80_20 <- alldatalist.video_80_20[[seed]][[prev_rate]][[model]]
  rf_fun        <- model_fun[[model]]
  
  rf.results <- tryCatch(
    rf_fun(results_80_20, source = source, comp = comp),
    error = function(e) {
      warning(paste("Failed for", seed, prev_rate, model, comp, source, ":", e$message))
      return(NULL)
    }
  )
  
  if (is.null(rf.results)) return(NULL)
  
  rf.results$prev_rate <- prev_rate
  rf.results$seed      <- seed
  rf.results$method    <- paste("RF", model, sep = "_")
  rf.results$source    <- source
  
  return(rf.results)
}

# Run all combinations
results_list <- mapply(
  run_rf,
  seed      = combinations$seed,
  prev_rate = combinations$prev_rate,
  model     = combinations$model,
  comp      = combinations$comp,
  source    = combinations$source,
  SIMPLIFY  = FALSE
)

# Rebuild the same nested list structure: [[comp]][[source]][[model]][[prev_rate]][[seed]]
rf.metrics.list <- list()

for (i in seq_along(results_list)) {
  res <- results_list[[i]]
  if (is.null(res)) next
  
  comp      <- combinations$comp[i]
  source    <- combinations$source[i]
  model     <- combinations$model[i]
  prev_rate <- combinations$prev_rate[i]
  seed      <- combinations$seed[i]
  
  rf.metrics.list[[comp]][[source]][[model]][[prev_rate]][[seed]] <- res
}

# Function to recursively flatten and bind nested lists
flatten_and_bind <- function(nested_list) {
  if (all(sapply(nested_list, is.data.frame))) {
    # If all elements are data frames, bind them
    do.call(rbind, nested_list)
  } else {
    # Otherwise, recursively process
    do.call(rbind, lapply(nested_list, flatten_and_bind))
  }
}

# bind the results of RF
rf.metrics.df <- flatten_and_bind(rf.metrics.list)

# remove rownames of metrics df
rownames(rf.metrics.df) <- NULL

# Factorize all columns of a dataframe
#rf.metrics.df[] <- lapply(rf.metrics.df, as.factor)

# reorder seed from 01, 02, to 10
rf.metrics.df$seed= ifelse(rf.metrics.df$seed == "seed_10", "seed_10", gsub("seed_", "seed_0", rf.metrics.df$seed))
rf.metrics.df= rf.metrics.df[order(rf.metrics.df$seed), ]

# replace bininter by predomics_bin and terinter by predomics_ter
# rf.metrics.df$method= ifelse(rf.metrics.df$method == "Indval", "RF_Indval", ifelse(rf.metrics.df$method == "bininter", "RF_predomics_bin", ifelse(rf.metrics.df$method == "terinter", "RF_predomics_ter", NA)))
# any(is.na(rf.metrics.df$method))

# save rf results
save(rf.metrics.df, file = paste(root_path, "rf_analysis_results_with_twinspan_80_20.Rda", sep="/"))

# # plotting results of Random Forest for each metric
# #filter auc values from rf results
# rf.auc.df <- rf.metrics.df[rf.metrics.df$metrics=='auc', ]
# 
# #filter acc values from rf results
# rf.acc.df <- rf.metrics.df[rf.metrics.df$metrics=='accuracy', ]
# 
# #filter f1 values from rf results
# rf.f1.df <- rf.metrics.df[rf.metrics.df$metrics=='f1', ]
# 
# 
# rf.auc.plot <- ggplot(data = rf.auc.df, aes(x = seed, y = value, group=method)) +
#   geom_line(aes(color=method), position = position_jitterdodge(dodge.width = 0.3), linewidth = 1, alpha = 0.5) +
#   geom_point(size = 1, alpha = 0.5, aes(color=method)) +
#   #geom_text(data=ibest, aes(x=x,y=y, label=label, colour=algorithm), size=3, inherit.aes = FALSE) +
#   ylab("auc") +
#   xlab("seed") +
#   ggtitle("Generalization performance") +
#   theme_bw() +
#   facet_grid(prev_rate+source~comparison) +
#   theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1, size = 8))
# 
# rf.acc.plot <-  ggplot(data = rf.acc.df, aes(x = seed, y = value, group=method)) +
#   geom_line(aes(color=method), position = position_jitterdodge(dodge.width = 0.3), linewidth = 1, alpha = 0.5) +
#   geom_point(size = 1, alpha = 0.5, aes(color=method)) +
#   #geom_text(data=ibest, aes(x=x,y=y, label=label, colour=algorithm), size=3, inherit.aes = FALSE) +
#   ylab("acc") +
#   xlab("seed") +
#   ggtitle("Generalization performance") +
#   theme_bw() +
#   facet_grid(prev_rate+source~comparison) +
#   theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1, size = 8))
# 
# rf.f1.plot <-  ggplot(data = rf.f1.df, aes(x = seed, y = value, group=method)) +
#   geom_line(aes(color=method), position = position_jitterdodge(dodge.width = 0.3), linewidth = 1, alpha = 0.5) +
#   geom_point(size = 1, alpha = 0.5, aes(color=method)) +
#   #geom_text(data=ibest, aes(x=x,y=y, label=label, colour=algorithm), size=3, inherit.aes = FALSE) +
#   ylab("f1") +
#   xlab("seed") +
#   ggtitle("Generalization performance") +
#   theme_bw() +
#   facet_grid(prev_rate+source~comparison) +
#   theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1, size = 8))
# 
# # plot all metrics performances of RF
# # my_comparisons <- list(c("Indval", "predomics_bin"), c("Indval", "predomics_ter"), c("predomics_bin", "predomics_ter"))
# 
# ggplot(data=rf.auc.df, aes(x= comparison, y= value, color=method))+
#   geom_boxplot(alpha = 0.6)+
#   geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
#   labs(x = "Comparison", y = "AUC") +
#   facet_grid(prev_rate~source) +
#   ylim(c(0.5,1.5))+
#   theme_bw() +
#   stat_compare_means(aes(group = method, label = paste0("p = ", after_stat(p.format))), label.y = c(1.1, 1.2, 1.3, 1.4)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12))+
#   scale_size(range = c(1,5), limits = c(1,80))
# 
# ggplot(data=rf.acc.df, aes(x= comparison, y= value, color=method, label=featNum))+
#   geom_boxplot(alpha = 0.6)+
#   geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
#   labs(x = "Comparison", y = "ACC") +
#   facet_grid(prev_rate~source) +
#   ylim(c(0.5,1.5))+
#   theme_bw() +
#   stat_compare_means(aes(group = method, label = paste0("p = ", after_stat(p.format))), label.y = c(1.1, 1.2, 1.3, 1.4)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12))+
#   scale_size(range = c(1,5), limits = c(1,80))
# 
# ggplot(data=rf.f1.df, aes(x= comparison, y= value, color=method))+
#   geom_boxplot(alpha = 0.6)+
#   geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
#   labs(x = "Comparison", y = "f1") +
#   facet_grid(prev_rate~source) +
#   ylim(c(0.5,1.5))+
#   theme_bw() +
#   stat_compare_means(aes(group = method, label = paste0("p = ", after_stat(p.format))), label.y = c(1.1, 1.2, 1.3, 1.4)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12))+
#   scale_size(range = c(1,5), limits = c(1,80))
