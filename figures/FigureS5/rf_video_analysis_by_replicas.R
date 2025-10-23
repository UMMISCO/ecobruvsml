# RF analysis
# This script is used to identify the set of key species that characterize sites (Barrier, Bay, Lagoon) with efficiency 
if(!require("randomForest")) install.packages("randomForest") ; library(randomForest)
if(!require("caret")) install.packages("caret") ; library(caret)
if (!require("readr")) install.packages("readr"); library(readr)
if (!require("ROCR")) install.packages("ROCR"); library(ROCR)
if (!require("stringr")) install.packages("stringr"); library(stringr)
if (!require("ggpubr")) install.packages("ggpubr"); library(ggpubr)
if (!require("purrr")) install.packages("purrr"); library(purrr)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)

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
  # get input data
  X.train= data_source$Comp_data$X_train
  y.train= data_source$Comp_data$y_train
  X.test= data_source$Comp_data$X_test
  y.test= data_source$Comp_data$y_test
  # get indicator species 
  indval_species= data_source$indval_pval[data_source$indval_pval$IsIndSp==1, ]
  #Get the species in the indval results on train dataset
  indval_species.names <- as.character(indval_species$feature)
  #Subset the train data to indval_species.names
  rf.train_data.X <- X.train[,indval_species.names, drop=FALSE]
  # bind samples and targets in a dataframe
  rf.train_data <- as.data.frame(rf.train_data.X)
  #Add class variable to the X table
  rf.train_data$class <- as.factor(y.train)
  ##need to remove spaces in colnames (species names) before running the rf model (if not error)
  colnames(rf.train_data) <- gsub(" ",".", colnames(rf.train_data))
  set.seed(20)
  rf.train.model <- randomForest(class~., data = rf.train_data, proximity=TRUE) 
  
  # Identify columns in rf.train_data.X that are missing in X.test
  missing_columns <- setdiff(colnames(rf.train_data.X), colnames(X.test))
  # Add missing columns to X.test with values set to zero
  
  for (col in missing_columns) {
    X.test[[col]] <- 0
  }
  
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
  indval.pred_test.auc= performance(indval_rf.test.prediction, measure = "auc")@y.values[[1]] 
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
  # get input data
  X.train= data_source$Comp_data$X_train
  y.train= data_source$Comp_data$y_train
  X.test= data_source$Comp_data$X_test
  y.test= data_source$Comp_data$y_test
  # get indicator species 
  pred_species= data_source$pred_out_fbm[data_source$pred_out_fbm$IsIndSp==1, ]
  #Get the species in the indval results on train dataset
  pred_species.names <- as.character(pred_species$feature)
  #Subset the train data to indval_species.names
  rf.train_data.X <- X.train[, pred_species.names, drop=FALSE]
  # bind samples and targets in a dataframe
  rf.train_data <- as.data.frame(rf.train_data.X)
  #Add class variable to the X table
  rf.train_data$class <- as.factor(y.train)
  ##need to remove spaces in colnames (species names) before running the rf model (if not error)
  colnames(rf.train_data) <- gsub(" ",".", colnames(rf.train_data))
  set.seed(20)
  rf.train.model <- randomForest(class~., data = rf.train_data, proximity=TRUE) 
  
  # Identify columns in rf.train_data.X that are missing in X.test
  missing_columns <- setdiff(colnames(rf.train_data.X), colnames(X.test))
  # Add missing columns to X.test with values set to zero
  
  for (col in missing_columns) {
    X.test[[col]] <- 0
  }
  
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
  pred_test.auc= performance(pred_rf.test.prediction, measure = "auc")@y.values[[1]] 
  pred_test.auc= round(pred_test.auc, digits = 3)
  
  #get metric performances of rf with indval species in testing data
  rf.results <- data.frame(metrics= c("auc","accuracy","f1"),
                           value= c(pred_test.auc, pred_test.acc, pred_test.f1),
                           comparison=comp,
                           source= source,
                           featNum= ncol(rf.train_data.X))
  return(rf.results)
}

# load result list of split by replicas analysis (Indval/terinter/bininter)
load("/data/projects/aime/analyses/bruvs/results_replicas_analysis/all_analysis_results_split_replicas.Rda")

# Create the parameter combinations
param_grid <- expand.grid(
  replica_comp = names(alldatalist.video_split_replicas),
  prev_rate = c("prev_0", "prev_10"),
  model = c("Indval", "bininter", "terinter"),
  comp = c("BAR_BAY", "BAR_LAG", "BAY_LAG", "INSHORE_OFFSHORE"),
  stringsAsFactors = FALSE
)

# Initialize empty list
rf.metrics.list <- list()

# Main processing function
process_entry <- function(replica_comp, prev_rate, model, comp) {
  data_model <- alldatalist.video_split_replicas[[replica_comp]][[prev_rate]][[model]]
  method_label <- paste("RF", model, sep = "_")
  
  analysis_fn <- if (model == "Indval") rf_indval_analysis else rf_pred_analysis
  
  for (source in c("maxN", "presAbs")) {
    result <- analysis_fn(data_model, source = source, comp = comp)
    result$prev_rate <- prev_rate
    result$method <- method_label
    result$replica_comp <- replica_comp
    
    # Safe nested list assignment
    rf.metrics.list[[replica_comp]][[comp]][[source]][[model]][[prev_rate]] <<- result
  }
}

# Apply function using pwalk
pwalk(param_grid, process_entry)

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

# save rf results
save(rf.metrics.df, file = "/data/projects/aime/analyses/bruvs/rf_analysis_by_replicas/rf_analysis_results_split_by_replicas.Rda")
