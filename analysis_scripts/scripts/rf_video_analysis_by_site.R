# RF analysis
# This script is used to identify the set of key species that characterize sites (Barrier, Bay, Lagoon) and zones (Inshore, Offshore) with efficiency 

if(!require("randomForest")) install.packages("randomForest") ; library(randomForest)
if(!require("caret")) install.packages("caret") ; library(caret)
if (!require("readr")) install.packages("readr"); library(readr)
if (!require("ROCR")) install.packages("ROCR"); library(ROCR)
if (!require("stringr")) install.packages("stringr"); library(ROCR)


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
  rf.train_data.X <- X.train[,indval_species.names]
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
  rf.results <- data.frame(rf_metric= c("auc","accuracy","f1"),
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
  rf.train_data.X <- X.train[, pred_species.names]
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
  rf.results <- data.frame(rf_metric= c("auc","accuracy","f1"),
                           value= c(pred_test.auc, pred_test.acc, pred_test.f1),
                           comparison=comp,
                           source= source,
                           featNum= ncol(rf.train_data.X))
  return(rf.results)
}


root_path="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F"

# Get a list of all directories within the root_path recursively
all_dirs <- list.dirs(path = root_path, full.names = TRUE, recursive = TRUE)

# Filter directories that contain "output_data" in their names
output_data_dirs <- grep("output_data/analysis_80_20", all_dirs, value = TRUE)

files <- list.files(output_data_dirs, recursive = TRUE)

alldatalist.video_80_20 <- list()

for(f in files)
{
  #print(f)
  # Extract method and seed from result file name
  matches <- str_match(f, "^([^_]+).*_data_([^.]+)")
  # Store into separate variables
  method <- matches[, 2]  # First part (before the first "_")
  prev_seed <- matches[, 3]  # Part after "data_"
  # Extract prev_rate and seed
  prevalence_rate <- sub("^(prev_\\d+)_.*", "\\1", prev_seed)
  seed <- sub(".*_(seed_\\d+)$", "\\1", prev_seed)
  # get the path of each result
  file_dir_path <- paste(root_path, paste(method, "output_data", sep="_"), "analysis_80_20", sep="/")
  # load results data
  fobj <- load(paste(file_dir_path,f,sep = "/"))
  # get result object
  fobj <- mget(fobj)
  # save result object
  alldatalist.video_80_20[[seed]][[prevalence_rate]][[method]] <- fobj
}

# save Indval/Predomics analysis results data
save(alldatalist.video_80_20, file = paste(root_path, "all_analysis_results_80_20.Rda", sep="/"))

# run random forest on each with indicator species of each method 
rf.metrics.list <- list()
#results_80_20 <- list()
for (seed in names(alldatalist.video_80_20)) {
  for (prev_rate in c("prev_0", "prev_10")) {
    for (model in c("Indval","bininter","terinter")) {
      for (comp in c("BAR_BAY", "BAR_LAG", "BAY_LAG", "INSHORE_OFFSHORE")){
        if (model=="Indval"){
          # Run random Forest in abundance for each comparison
          results_80_20 = alldatalist.video_80_20[[seed]][[prev_rate]][[model]]
          rf.results= rf_indval_analysis(results_80_20, source="maxN", comp)
          rf.results$prev_rate= prev_rate
          rf.results$seed = seed
          rf.results$method= model
          ## save results
          rf.metrics.list[[comp]][["maxN"]][[model]][[prev_rate]][[seed]] <- rf.results
          
          # Run random Forest in presence/absence for each comparison
          rf.results= rf_indval_analysis(results_80_20, source="presAbs", comp)
          rf.results$prev_rate= prev_rate
          rf.results$seed = seed
          rf.results$method= model
          ## save results
          rf.metrics.list[[comp]][["presAbs"]][[model]][[prev_rate]][[seed]] <- rf.results
          
        } else if (model %in% c("bininter", "terinter")){
          
          # Run random Forest in abundance for each comparison
          results_80_20= alldatalist.video_80_20[[seed]][[prev_rate]][[model]]
          rf.results= rf_pred_analysis(results_80_20, source="maxN", comp)
          rf.results$prev_rate= prev_rate
          rf.results$seed = seed
          rf.results$method= model
          ## save results
          rf.metrics.list[[comp]][["maxN"]][[model]][[prev_rate]][[seed]] <- rf.results
          
          # Run random Forest in presence/absence for each comparison
          rf.results= rf_pred_analysis(results_80_20, source="presAbs", comp)
          rf.results$prev_rate= prev_rate
          rf.results$seed = seed
          rf.results$method= model
          ## save results
          rf.metrics.list[[comp]][["presAbs"]][[model]][[prev_rate]][[seed]] <- rf.results
          
        }
      }
    }
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
rf.metrics.df$method= ifelse(rf.metrics.df$method == "Indval", "Indval", ifelse(rf.metrics.df$method == "bininter", "predomics_bin", ifelse(rf.metrics.df$method == "terinter", "predomics_ter", NA)))
any(is.na(rf.metrics.df$method))

# plotting results of Random Forest for each metric
#filter auc values from rf results
rf.auc.df <- rf.metrics.df[rf.metrics.df$rf_metric=='auc', ]

#filter acc values from rf results
rf.acc.df <- rf.metrics.df[rf.metrics.df$rf_metric=='accuracy', ]

#filter f1 values from rf results
rf.f1.df <- rf.metrics.df[rf.metrics.df$rf_metric=='f1', ]


rf.auc.plot <- ggplot(data = rf.auc.df, aes(x = seed, y = value, group=method)) +
  geom_line(aes(color=method), position = position_jitterdodge(dodge.width = 0.3), size = 1, alpha = 0.5) +
  geom_point(size = 1, alpha = 0.5, aes(color=method)) +
  #geom_text(data=ibest, aes(x=x,y=y, label=label, colour=algorithm), size=3, inherit.aes = FALSE) +
  ylab("auc") +
  xlab("seed") +
  ggtitle("Generalization performance") +
  theme_bw() +
  facet_grid(prev_rate+source~comparison) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1, size = 8))

rf.acc.plot <-  ggplot(data = rf.acc.df, aes(x = seed, y = value, group=method)) +
  geom_line(aes(color=method), position = position_jitterdodge(dodge.width = 0.3), size = 1, alpha = 0.5) +
  geom_point(size = 1, alpha = 0.5, aes(color=method)) +
  #geom_text(data=ibest, aes(x=x,y=y, label=label, colour=algorithm), size=3, inherit.aes = FALSE) +
  ylab("acc") +
  xlab("seed") +
  ggtitle("Generalization performance") +
  theme_bw() +
  facet_grid(prev_rate+source~comparison) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1, size = 8))

rf.f1.plot <-  ggplot(data = rf.f1.df, aes(x = seed, y = value, group=method)) +
  geom_line(aes(color=method), position = position_jitterdodge(dodge.width = 0.3), size = 1, alpha = 0.5) +
  geom_point(size = 1, alpha = 0.5, aes(color=method)) +
  #geom_text(data=ibest, aes(x=x,y=y, label=label, colour=algorithm), size=3, inherit.aes = FALSE) +
  ylab("f1") +
  xlab("seed") +
  ggtitle("Generalization performance") +
  theme_bw() +
  facet_grid(prev_rate+source~comparison) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1, size = 8))

# plot all metrics performances of RF
# my_comparisons <- list(c("Indval", "predomics_bin"), c("Indval", "predomics_ter"), c("predomics_bin", "predomics_ter"))

ggplot(data=rf.auc.df, aes(x= comparison, y= value, color=method))+
  geom_boxplot(alpha = 0.6)+
  labs(x = "Comparison", y = "AUC") +
  facet_grid(prev_rate~source) +
  ylim(c(0,1.5))+
  theme_bw() +
  stat_compare_means(aes(group = method), label.y = c(1.1, 1.2, 1.3, 1.4)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12))

ggplot(data=rf.acc.df, aes(x= comparison, y= value, color=method))+
  geom_boxplot(alpha = 0.6)+
  labs(x = "Comparison", y = "ACC") +
  facet_grid(prev_rate~source) +
  ylim(c(0,1.5))+
  theme_bw() +
  stat_compare_means(aes(group = method), label.y = c(1.1, 1.2, 1.3, 1.4)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12))

ggplot(data=rf.f1.df, aes(x= comparison, y= value, color=method))+
  geom_boxplot(alpha = 0.6)+
  labs(x = "Comparison", y = "f1") +
  facet_grid(prev_rate~source) +
  ylim(c(0,1.5))+
  theme_bw() +
  stat_compare_means(aes(group = method), label.y = c(1.1, 1.2, 1.3, 1.4)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12))
