library(ggplot2) # for plots
library(gridExtra) # for multiple plots
library(patchwork) # to display multiples graphs in the same figure
library(dplyr)
library(ggVennDiagram)
library(reshape2)
library(ggpubr)
library(stringr)
if(!require("randomForest")) install.packages("randomForest") ; library(randomForest)
if(!require("caret")) install.packages("caret") ; library(caret)
if (!require("readr")) install.packages("readr"); library(readr)
if (!require("ROCR")) install.packages("ROCR"); library(ROCR)
if (!require("stringr")) install.packages("stringr"); library(stringr)
if (!require("ggpubr")) install.packages("ggpubr"); library(ggpubr)
if (!require("purrr")) install.packages("purrr"); library(purrr)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)

root_path="/data/projects/aime/analyses/bruvs/FigureS5/results_replicas_analysis"

# Get a list of all directories within the root_path recursively
all_dirs <- list.dirs(path = root_path, full.names = TRUE, recursive = TRUE)

# Filter directories that contain "output_data" in their names
output_data_dirs <- grep("inter_output_data", all_dirs, value = TRUE)

files <- list.files(output_data_dirs, recursive = TRUE)

# Initialize a named list to store all result lists
alldatalist.video_split_replicas_predomics <- list()

for(f in files)
{
  # Extract the comparison and prevalence_rate from the filename
  matches <- str_match(f, "^([^_]+).*_replicas_([^_]+_vs_[^_]+)_prev_(\\d+)\\.Rda$")
  
  if (!is.na(matches[1, 1])) { 
    method <- matches[1, 2]
    replicas_comp <- matches[1, 3]  # e.g., "ABORE_vs_MBERE"
    prevalence_rate <- paste0("prev_", matches[1, 4])  # e.g., "prev_0"
    
    # Ensure the sublist exists
    if (is.null(alldatalist.video_split_replicas_predomics[[replicas_comp]])) {
      alldatalist.video_split_replicas_predomics[[replicas_comp]] <- list()
    }
    
    # Get the path of each result
    file_dir_path <- paste(root_path, paste(method, "output_data", sep="_"), sep="/")
    
    # Load data
    fobj_name <- load(paste(file_dir_path, f, sep = "/"))
    fobj <- mget(fobj_name)
    
    # Save to the dynamically-named sublist
    alldatalist.video_split_replicas_predomics[[replicas_comp]][[prevalence_rate]][[method]] <- fobj
  }
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

# define the combinations
param_grid <- expand.grid(
  i = names(alldatalist.video_split_replicas_predomics),
  j = c("prev_0", "prev_10"),
  k = c("bininter", "terinter"),
  l = c("BAR_BAY", "BAR_LAG", "BAY_LAG", "INSHORE_OFFSHORE"),
  stringsAsFactors = FALSE
)

all_species.list <- list()
 
# Loop using pmap
all_species.list <- pmap(param_grid, function(i, j, k, l) {

  # Retrieve model outputs
  bin_model <- alldatalist.video_split_replicas_predomics[[i]][[j]][[k]][["predout.bin"]][[l]]
  maxn_model <- alldatalist.video_split_replicas_predomics[[i]][[j]][[k]][["predout.maxn"]][[l]]

  # best model species for maxn
  best_model.species.maxn <- data.frame(
    feature = maxn_model$digest$best$model$names_,
    comparison = l,
    stringsAsFactors = FALSE
  )

  # best model species for bin
  best_model.species.bin <- data.frame(
    feature = bin_model$digest$best$model$names_,
    comparison = l,
    stringsAsFactors = FALSE
  )

  # update presAbs
  presAbs <- bin_model[["pred_out_fbm"]]
  presAbs$IsIndSp <- as.integer(presAbs$feature %in% best_model.species.bin$feature)
  presAbs$comparison = l
  presAbs$source= paste("BestModel", k, sep = "_")
  presAbs$data= "pres/abs"
  bin_model$pred_out_bestModel <- presAbs

  # update maxN
  maxN <- maxn_model$pred_out_fbm
  maxN$IsIndSp <- as.integer(maxN$feature %in% best_model.species.maxn$feature)
  maxN$comparison = l
  maxN$source= paste("BestModel", k, sep = "_")
  maxN$data= "maxN"
  maxn_model$pred_out_bestModel <- maxN
  
  # Save updated models back into alldatalist
  alldatalist.video_split_replicas_predomics[[i]][[j]][[k]][["predout.bin"]][[l]] <<- bin_model
  alldatalist.video_split_replicas_predomics[[i]][[j]][[k]][["predout.maxn"]][[l]] <<- maxn_model
  

  # bind and return combined result
  bind_rows(presAbs, maxN) %>%
    dplyr::mutate(
      prev_rate = j,
      replica_comp = i
    )
})

# combine into a single data frame if needed
all_species.df <- bind_rows(all_species.list)

# get only indicator species
indics_species.df= all_species.df[all_species.df$IsIndSp==1,]

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
  pred_species= data_source$pred_out_bestModel[data_source$pred_out_bestModel$IsIndSp==1, ]
  #Get the species in the indval results on train dataset
  pred_species.names <- as.character(pred_species$feature)
  #Subset the train data to pred_species.names
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

# run random forest on each with indicator species of each method 
rf.metrics.list <- list()

# Main processing function
process_entry <- function(i, j, k, l) {
  data_model <- alldatalist.video_split_replicas_predomics[[i]][[j]][[k]]
  method_label <- paste("RF_bestModel", k, sep = "_")
  
  for (source in c("maxN", "presAbs")) {
    result <- rf_pred_analysis(data_model, source = source, comp = l)
    result$prev_rate <- j
    result$method <- method_label
    result$replica_comp <- i
    
    rf.metrics.list[[i]][[l]][[source]][[k]][[j]] <<- result
  }
}

# Apply function using pwalk
pwalk(param_grid, process_entry)

# bind the results of RF
rf.metrics.bestModel.df <- flatten_and_bind(rf.metrics.list)

# remove rownames of metrics df
rownames(rf.metrics.bestModel.df) <- NULL

# save the data object of Best models Predomics
save(alldatalist.video_split_replicas_predomics, file = paste(root_path,"all_analysis_results_split_replicas_BestModels_Predomics.Rda", sep = "/"))

# # save Indval/Predomics analysis results data
# save(alldatalist.video_split_replicas, file = paste(root_path,"all_analysis_results_split_replicas.Rda", sep = "/"))

# save analysis results for all species or only indics species with Indval, Predomics_bininter, Predomics_terinter
write.csv(all_species.df, file = paste(root_path, "all_species_analysis_split_by_replicas_BestModels_Predomics.csv", sep="/"), row.names = FALSE)
write.csv(indics_species.df, file = paste(root_path, "indics_species_analysis_split_by_replicas_BestModels_Predomics.csv", sep="/"), row.names = FALSE)

# save rf results on Best Models
save(rf.metrics.bestModel.df, file = "/data/projects/aime/analyses/bruvs/rf_analysis_by_replicas/rf_analysis_results_split_by_replicas_BestModel_Predomics.Rda")
