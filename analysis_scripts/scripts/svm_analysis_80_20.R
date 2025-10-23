# SVM analysis
# This script is used to identify the set of key species that characterize sites (Barrier, Bay, Lagoon) with efficiency 

if(!require("predomics")) devtools::install_github("predomics/predomicspkg", dependencies = TRUE) ; library(predomics)
if(!require("e1071")) install.packages("e1071") ; library(e1071) # will allow us to run SVM model
if(!require("caret")) install.packages("caret") ; library(caret) 
if (!require("readr")) install.packages("readr"); library(readr)
if (!require("ROCR")) install.packages("ROCR"); library(ROCR)
if (!require("pROC")) install.packages("pROC"); library(pROC)
if (!require("stringr")) install.packages("stringr"); library(stringr)
if (!require("ggpubr")) install.packages("ggpubr"); library(ggpubr)


# svm analysis on indval results
svm_indval_analysis <-function (data, source, comp){
  
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
  svm.train_data.X <- X.train[,indval_species.names]
  # bind samples and targets in a dataframe
  svm.train_data <- as.data.frame(svm.train_data.X)
  #Add class variable to the X table
  svm.train_data$class <- as.factor(y.train)
  ##need to remove spaces in colnames (species names) before running the svm model (if not error)
  colnames(svm.train_data) <- gsub(" ",".", colnames(svm.train_data))
  set.seed(20)
  
  # sota svm
  # clf_sota <- sota.svm(c(1:30),
  #                     max.nb.features = 1000,
  #                     language="svm",
  #                     nCores = 2,
  #                     seed=20,
  #                     evalToFit = "auc_",
  #                     objective = "auc",
  #                     experiment.save = "nothing")
  # 
  # res_clf_svm <- fit(X = t(svm.train_data.X), y = svm.train_data$class, clf = clf_sota, cross.validate = TRUE, nfolds = 10)
  
  svm.train.model <- svm(class~., data = svm.train_data, scale = FALSE, kernel = "radial", cost = 5) 
  
  # Identify columns in svm.train_data.X that are missing in X.test
  missing_columns <- setdiff(colnames(svm.train_data.X), colnames(X.test))
  # Add missing columns to X.test with values set to zero
  
  for (col in missing_columns) {
    X.test[[col]] <- 0
  }
  
  # Ensure svm.test_data has the same columns as svm.train_data.X, in the same order
  svm.test_data <- X.test[, colnames(svm.train_data.X), drop = FALSE]
  #Add class variable to the X table
  svm.test_data$class <- as.factor(y.test)
  ##need to remove spaces in colnames (species names) before running the svm model (if not error)
  colnames(svm.test_data) <- gsub(" ",".", colnames(svm.test_data))
  
  # make prediction on test data
  indval_test <- predict(svm.train.model, newdata=svm.test_data, type= "class")
  
  # The prediction to compute the confusion matrix and see the accuracy score 
  indval.test.cm=confusionMatrix(table(indval_test,svm.test_data$class))
  
  # get acc, auc and f1 metrics of svm on indval species using test data
  indval.pred_test.acc <- round(indval.test.cm$overall['Accuracy'], digits = 3)
  indval.pred_test.f1 <- round(indval.test.cm$byClass['F1'], digits = 3)
  
  # get the auc of svm model on test data
  indval_svm.test.pred=predict(svm.train.model, newdata = svm.test_data, type = "prob")
  indval_svm.test.prediction = prediction(as.numeric(indval_svm.test.pred), as.numeric(svm.test_data$class))
  indval.pred_test.auc= performance(indval_svm.test.prediction, measure = "auc")@y.values[[1]] 
  indval.pred_test.auc= round(indval.pred_test.auc, digits = 3)
  
  #get metric performances of svm with indval species in testing data
  svm.results <- data.frame(metrics= c("auc","accuracy","f1"),
                           value= c(indval.pred_test.auc, indval.pred_test.acc, indval.pred_test.f1),
                           comparison=comp,
                           source= source,
                           featNum= ncol(svm.train_data.X))
  return(svm.results)
}

# svm analysis on predomics results
svm_pred_analysis <-function (data, source, comp){
  
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
  svm.train_data.X <- X.train[, pred_species.names]
  # bind samples and targets in a dataframe
  svm.train_data <- as.data.frame(svm.train_data.X)
  #Add class variable to the X table
  svm.train_data$class <- as.factor(y.train)
  ##need to remove spaces in colnames (species names) before running the svm model (if not error)
  colnames(svm.train_data) <- gsub(" ",".", colnames(svm.train_data))
  set.seed(20)
  svm.train.model <- svm(class~., data = svm.train_data, scale = FALSE, kernel = "radial", cost = 5)
  
  # Identify columns in svm.train_data.X that are missing in X.test
  missing_columns <- setdiff(colnames(svm.train_data.X), colnames(X.test))
  # Add missing columns to X.test with values set to zero
  
  for (col in missing_columns) {
    X.test[[col]] <- 0
  }
  
  # Ensure svm.test_data has the same columns as svm.train_data.X, in the same order
  svm.test_data <- X.test[, colnames(svm.train_data.X), drop = FALSE]
  #Add class variable to the X table
  svm.test_data$class <- as.factor(y.test)
  ##need to remove spaces in colnames (species names) before running the svm model (if not error)
  colnames(svm.test_data) <- gsub(" ",".", colnames(svm.test_data))
  
  # make prediction on test data
  pred_test <- predict(svm.train.model, newdata = svm.test_data, type= "class")
  
  # The prediction to compute the confusion matrix and see the accuracy score 
  pred.test.cm=confusionMatrix(table(pred_test,svm.test_data$class)) 
  
  # get acc, auc and f1 metrics of svm on predomics species using test data
  pred_test.acc <- round(pred.test.cm$overall['Accuracy'], digits = 3)
  pred_test.f1 <- round(pred.test.cm$byClass['F1'], digits = 3)
  
  # get the auc of svm model on test data
  pred_svm.test=predict(svm.train.model, newdata = svm.test_data, type = "prob")
  pred_svm.test.prediction = prediction(as.numeric(pred_svm.test), as.numeric(svm.test_data$class))
  pred_test.auc= performance(pred_svm.test.prediction, measure = "auc")@y.values[[1]] 
  pred_test.auc= round(pred_test.auc, digits = 3)
  
  #get metric performances of svm with indval species in testing data
  svm.results <- data.frame(metrics= c("auc","accuracy","f1"),
                           value= c(pred_test.auc, pred_test.acc, pred_test.f1),
                           comparison=comp,
                           source= source,
                           featNum= ncol(svm.train_data.X))
  return(svm.results)
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
# save(alldatalist.video_80_20, file = paste(root_path, "all_analysis_results_80_20.Rda", sep="/"))

# run SVM with indicator species of each method 
svm.metrics.list <- list()
#results_80_20 <- list()
for (seed in names(alldatalist.video_80_20)) {
  for (prev_rate in c("prev_0", "prev_10")) {
    for (model in c("Indval","bininter","terinter")) {
      for (comp in c("BAR_BAY", "BAR_LAG", "BAY_LAG", "INSHORE_OFFSHORE")){
        if (model=="Indval"){
          # Run random Forest in abundance for each comparison
          results_80_20 = alldatalist.video_80_20[[seed]][[prev_rate]][[model]]
          svm.results= svm_indval_analysis(results_80_20, source="maxN", comp)
          svm.results$prev_rate= prev_rate
          svm.results$seed = seed
          svm.results$method= paste("SVM", model, sep = "_")
          ## save results
          svm.metrics.list[[comp]][["maxN"]][[model]][[prev_rate]][[seed]] <- svm.results
          
          # Run random Forest in presence/absence for each comparison
          svm.results= svm_indval_analysis(results_80_20, source="presAbs", comp)
          svm.results$prev_rate= prev_rate
          svm.results$seed = seed
          svm.results$method= paste("SVM", model, sep = "_")
          ## save results
          svm.metrics.list[[comp]][["presAbs"]][[model]][[prev_rate]][[seed]] <- svm.results
          
        } else if (model %in% c("bininter", "terinter")){
          
          # Run random Forest in abundance for each comparison
          results_80_20= alldatalist.video_80_20[[seed]][[prev_rate]][[model]]
          svm.results= svm_pred_analysis(results_80_20, source="maxN", comp)
          svm.results$prev_rate= prev_rate
          svm.results$seed = seed
          svm.results$method= paste("SVM", model, sep = "_")
          ## save results
          svm.metrics.list[[comp]][["maxN"]][[model]][[prev_rate]][[seed]] <- svm.results
          
          # Run random Forest in presence/absence for each comparison
          svm.results= svm_pred_analysis(results_80_20, source="presAbs", comp)
          svm.results$prev_rate= prev_rate
          svm.results$seed = seed
          svm.results$method= paste("SVM", model, sep = "_")
          ## save results
          svm.metrics.list[[comp]][["presAbs"]][[model]][[prev_rate]][[seed]] <- svm.results
          
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

# bind the results of svm
svm.metrics.df <- flatten_and_bind(svm.metrics.list)

# remove rownames of metrics df
rownames(svm.metrics.df) <- NULL

# Factorize all columns of a dataframe
#svm.metrics.df[] <- lapply(svm.metrics.df, as.factor)

# reorder seed from 01, 02, to 10
svm.metrics.df$seed= ifelse(svm.metrics.df$seed == "seed_10", "seed_10", gsub("seed_", "seed_0", svm.metrics.df$seed))
svm.metrics.df= svm.metrics.df[order(svm.metrics.df$seed), ]

# replace bininter by predomics_bin and terinter by predomics_ter
# svm.metrics.df$method= ifelse(svm.metrics.df$method == "Indval", "SVM_Indval", ifelse(svm.metrics.df$method == "bininter", "SVM_predomics_bin", ifelse(svm.metrics.df$method == "terinter", "SVM_predomics_ter", NA)))
# any(is.na(svm.metrics.df$method))

# save svm results
save(svm.metrics.df, file = paste(root_path, "svm_analysis_results_80_20.Rda", sep="/"))

# plotting results of Random Forest for each metric
#filter auc values from svm results
svm.auc.df <- svm.metrics.df[svm.metrics.df$metrics=='auc', ]

#filter acc values from svm results
svm.acc.df <- svm.metrics.df[svm.metrics.df$metrics=='accuracy', ]

#filter f1 values from svm results
svm.f1.df <- svm.metrics.df[svm.metrics.df$metrics=='f1', ]


svm.auc.plot <- ggplot(data = svm.auc.df, aes(x = seed, y = value, group=method)) +
  geom_line(aes(color=method), position = position_jitterdodge(dodge.width = 0.3), size = 1, alpha = 0.5) +
  geom_point(size = 1, alpha = 0.5, aes(color=method)) +
  #geom_text(data=ibest, aes(x=x,y=y, label=label, colour=algorithm), size=3, inherit.aes = FALSE) +
  ylab("auc") +
  xlab("seed") +
  ggtitle("Generalization performance") +
  theme_bw() +
  facet_grid(prev_rate+source~comparison) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 9))

svm.acc.plot <-  ggplot(data = svm.acc.df, aes(x = seed, y = value, group=method)) +
  geom_line(aes(color=method), position = position_jitterdodge(dodge.width = 0.3), size = 1, alpha = 0.5) +
  geom_point(size = 1, alpha = 0.5, aes(color=method)) +
  #geom_text(data=ibest, aes(x=x,y=y, label=label, colour=algorithm), size=3, inherit.aes = FALSE) +
  ylab("acc") +
  xlab("seed") +
  ggtitle("Generalization performance") +
  theme_bw() +
  facet_grid(prev_rate+source~comparison) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 9))

svm.f1.plot <-  ggplot(data = svm.f1.df, aes(x = seed, y = value, group=method)) +
  geom_line(aes(color=method), position = position_jitterdodge(dodge.width = 0.3), size = 1, alpha = 0.5) +
  geom_point(size = 1, alpha = 0.5, aes(color=method)) +
  #geom_text(data=ibest, aes(x=x,y=y, label=label, colour=algorithm), size=3, inherit.aes = FALSE) +
  ylab("f1") +
  xlab("seed") +
  ggtitle("Generalization performance") +
  theme_bw() +
  facet_grid(prev_rate+source~comparison) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 9))

# plot all metrics performances of svm
# my_comparisons <- list(c("Indval", "predomics_bin"), c("Indval", "predomics_ter"), c("predomics_bin", "predomics_ter"))

ggplot(data=svm.auc.df, aes(x= comparison, y= value, color=method))+
  geom_boxplot(alpha = 0.6)+
  geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
  labs(x = "Comparison", y = "AUC") +
  facet_grid(prev_rate~source) +
  ylim(c(0.35,1.5))+
  theme_bw() +
  stat_compare_means(aes(group = method, label = paste0("p = ", after_stat(p.format))), label.y = c(1.1, 1.2, 1.3, 1.4)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12))+
  scale_size(range = c(1,5), limits = c(1,80))

ggplot(data=svm.acc.df, aes(x= comparison, y= value, color=method))+
  geom_boxplot(alpha = 0.6)+
  geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
  labs(x = "Comparison", y = "ACC") +
  facet_grid(prev_rate~source) +
  ylim(c(0.35,1.5))+
  theme_bw() +
  stat_compare_means(aes(group = method, label = paste0("p = ", after_stat(p.format))), label.y = c(1.1, 1.2, 1.3, 1.4)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12))+
  scale_size(range = c(1,5), limits = c(1,80))

ggplot(data=svm.f1.df, aes(x= comparison, y= value, color=method))+
  geom_boxplot(alpha = 0.6)+
  geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
  labs(x = "Comparison", y = "f1") +
  facet_grid(prev_rate~source) +
  ylim(c(0.35,1.5))+
  theme_bw() +
  stat_compare_means(aes(group = method, label = paste0("p = ", after_stat(p.format))), label.y = c(1.1, 1.2, 1.3, 1.4)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12))+
  scale_size(range = c(1,5), limits = c(1,80))
