# To run this script, use this command: Rscript codePredomics_all_analyses_split_data.R terinter 10

# import libraries
library(readr) # for data import
library(DataExplorer) # for summary statistics
library(data.table) # for dcast => long to wide
library(momr) # for normalizing, visualizing and microbiome analyses
library(kableExtra) # for markdown table
library(ggplot2) # for plots
library(gridExtra) # for multiple plots
library(vegan) # for ecological analyses
library(viridis) # for colour maps
library(stats) # for euclidean distance computation
library(predomics) # loading predomics
library(dplyr)
library(caTools)
library(labdsv)

# define arguments for the script
args <- commandArgs(trailingOnly = TRUE)
algorithm_language <- args[1] 
prevalence_rate <- as.numeric(args[2])
#maxseed <- as.numeric(args[3])


dataset <- "/data/projects/aime/data/seamount/video/AllDataexportefixed.txt"
mydata=read.delim(dataset, header = TRUE, sep = "\t", quote="\"", dec=".",fill = TRUE)
head(mydata)

# Changing names in data
# mydata$G_SP=ess$SP
mydata$G_SP=gsub("_", " ",mydata$G_SP)
mydata$Transect=paste(substr(mydata$Transect,3,7))
mydata$Site=factor(mydata$Site,levels=c("ABAR","ALAG","ABAI","MBAR","MLAG","MBAI"))
mydata$Site <- do.call(rbind, lapply(mydata$Site, gsub, pattern = "BAI", replacement = "BAY"))[,1] 
mydata$Station <- do.call(rbind, lapply(mydata$Station, gsub, pattern = "BAI", replacement = "BAY"))[,1] 
head(mydata)

#summarize data by station and species
matablepivot=mydata %>%
  group_by(Site,G_SP)
matablepivot=as.data.frame(matablepivot)
matablepivot=matablepivot[c("Site","Station","G_SP","Transect","MaxN")]
head(matablepivot)

# Species community matrix
abund = reshape2::dcast(matablepivot, Station+ Transect + Site ~ G_SP, value.var = "MaxN")
abund[is.na(abund)] = 0    # fill with zeroes for summarise below
rownames(abund)=abund$Station
# Recreate habitat because lazy
abund$hab <- substring(abund$Site, 2)
acc <- abund

table(acc$hab)

##build list of comparisons 2x2
comp <- combn(x = unique(acc$hab), m = 2, simplify = FALSE)

##Add inshore/offshore
acc$hab_inoff <- ifelse(acc$hab %in% c("BAY","LAG"),"INSHORE",
                        ifelse(acc$hab %in% "BAR","OFFSHORE",NA))
which(is.na(acc$hab_inoff))
table(acc$hab, acc$hab_inoff)

# define a list for habitat INSHORE/OFFSHORE
ilevels_hab <- list("INSHORE", "OFFSHORE")

# Function to get sample data by prevalence
get_sample_by_prevalence <- function(abundance_matrix, prevalence_rate) {
  
  # Print the prevalence rate message
  message(paste("You set the prevalence rate to:", paste0(prevalence_rate, "%")))
  
  # Scale the prevalence rate to the range [0, 1]
  prevalence_rate <- prevalence_rate / 100
  
  # Check if prevalence rate is between 0 and 1
  if (prevalence_rate == 0) {
    message("You are using the overall data without any filtering.")
  }
  if (prevalence_rate == 1) {
    message("You are using species that appear in all samples.")
  }
  if (prevalence_rate < 0 || prevalence_rate > 1) {
    stop("Prevalence rate must be between 0 and 1.")
  }
  
  # Get the number of samples (rows) - fixed from columns to rows
  num_samples <- nrow(abundance_matrix)
  
  # Calculate the minimum number of non-zero occurrences required
  min_occurrences <- ceiling(prevalence_rate * num_samples)
  
  # Select columns (features) that meet the prevalence threshold
  valid_features <- colSums(abundance_matrix > 0) >= min_occurrences
  
  # Subset the abundance matrix based on valid features
  sampled_matrix <- abundance_matrix[, valid_features, drop = FALSE]
  
  # Remove samples (rows) that don't have any non-zero values
  sampled_matrix <- sampled_matrix[rowSums(sampled_matrix > 0) > 0, ]
  
  # Return the filtered matrix
  return(sampled_matrix)
}

# Function to prepare data for training
prepare_data_for_analysis <- function(sample_data, habitat_type, hab1, hab2){
  # Make sure hab1 and hab2 are different
  ifelse(hab1==hab2, 
         stop("It seems that the habitats you entered are the same."),
         {
           ##Exclude species absents in all subset samples
           sample_data <- sample_data[, colSums(sample_data>0)>0]
           
           ##Exclude unnecesary columns (keep abundance table only samples x species)
           X <- sample_data[,-match(c("Station","Transect", "Site","hab","hab_inoff"), colnames(sample_data))]
           sample_X <- get_sample_by_prevalence(X, prevalence_rate)
           
           # get the list of species that doesn't meet the prevalence rate
           unprevalent_species<- setdiff(colnames(X), colnames(sample_X))
           # get only species that meeting the prevalence rate
           sample_data <- sample_data[, !(colnames(sample_data) %in% unprevalent_species)]
           
           # Retained only samples from the sample filtered by prevalence 
           sample_data <- sample_data[rownames(sample_data) %in% rownames(sample_X),]
           # Exclude samples with no species
           sample_data <- sample_data[rowSums(sample_data>0)>0, ]
           
           ##Build the class vector to do the predomics comparison
           sample_data.class <- ifelse(sample_data[[habitat_type]]==hab1, -1,
                                       ifelse(sample_data[[habitat_type]]==hab2, 1, NA))
           which(is.na(sample_data.class))
           table(sample_data$hab, sample_data.class)
           ##Exclude unnecesary columns (keep abundance table only samples x species)
           sample_data.X <- sample_data[,-match(c("Station","Transect","Site","hab","hab_inoff"), colnames(sample_data))]
           
           # Return sample.info, sample.X, sample.class
           return(list(sample.info = sample_data, sample.X = sample_data.X, sample.class = sample_data.class))
         }
  )
  
}

# Function for Predomics analysis
predomics_analyis <- function(ilevels.df.train,ilevels.df.train.class,ilevels.df.test,ilevels.df.test.class, algorithm){
  
  ##Do the predomics analyses
  # set the learner
  ifelse(algorithm == "terinter",
         clf <- terBeam(language = 'terinter',
                        nCores = 1,
                        objective = "auc",
                        seed = 20,
                        plot = FALSE),
         ifelse(algorithm == "bininter",
                clf <- terBeam(language = 'bininter',
                               nCores = 1,
                               objective = "auc",
                               seed = 20,
                               plot = FALSE),
                clf <- 'NA'))
  
  #### Run the learner (training) ####
  runit = TRUE
  if(runit)
  {
    res_clf <- fit(X = t(ilevels.df.train), y = ilevels.df.train.class, clf = clf, cross.validate = TRUE, nfolds = 5);
    # save results
    # save(res_clf, clf, 
    #      file = paste("res_clf_terbeam", comparison_var, class1, paste0(class2, ".rda"), sep = "_"), 
    #      compression_level = 9)
  }
  
  # Build a master list to save predomics analysis results
  results_list <- list()
  
  # store data for training and testing
  results_list$Comp_data <- list(X_train= ilevels.df.train, y_train= ilevels.df.train.class,
                                 X_test= ilevels.df.test, y_test= ilevels.df.test.class)
  
  # save clf and classification results in the master list
  results_list$fit <- res_clf
  results_list$clf <- clf
  
  # save digest of results
  res_clf.dig <- digest(obj = res_clf, penalty = 0.75/100, plot = FALSE)
  results_list$digest <- res_clf.dig
  
  ### Testing the model in the testing dataset (20%)
  
  best.model.test <- evaluateModel(mod = res_clf.dig$best$model, X = ilevels.df.test, y = ilevels.df.test.class, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE, mode = "test")
  results_list$test_eval <- best.model.test
  
  #### Family of Best Models (FBM) ####
  
  # get the population of models scrambled by model size
  pop <- modelCollectionToPopulation(res_clf$classifier$models)
  pop.df <- populationToDataFrame(pop)
  results_list$model_pop <- pop.df
  
  # select the best
  fbm <- selectBestPopulation(pop)
  fbm.df <- populationToDataFrame(fbm)
  results_list$fbm <- fbm.df
  
  #### Family of Best Models (FBM) - Feature Annotation ####
  
  fa <- makeFeatureAnnot(pop = fbm,
                         X = t(ilevels.df.train),
                         y = ilevels.df.train.class,
                         clf = clf)
  
  results_list$featureAnnotFBM <- fa
  
  #### Feature importance ####
  
  feat1.import <- mergeMeltImportanceCV(list.results = list(terBeam = res_clf),
                                        filter.cv.prev = 0,
                                        min.kfold.nb = FALSE,
                                        learner.grep.pattern = "*",
                                        #nb.top.features = 50,
                                        feature.selection = rownames(fa$pop.noz),
                                        scaled.importance = TRUE,
                                        make.plot = FALSE,
                                        cv.prevalence = FALSE)
  
  
  feat2.import <- mergeMeltImportanceCV(list.results = list(terBeam = res_clf),
                                        filter.cv.prev = 0,
                                        min.kfold.nb = FALSE,
                                        learner.grep.pattern = "*",
                                        nb.top.features = 148,
                                        #feature.selection = rownames(fa$pop.noz),
                                        scaled.importance = TRUE,
                                        make.plot = FALSE,
                                        cv.prevalence = FALSE)
  
  
  results_list$FI_fmbFeats <- feat1.import
  results_list$FI_allfeat <- feat2.import
  return(results_list)
}

# Function to save predomics results object
save_pred_results <- function(results_pred, X, ilevels){
  
  # get only species from FBM
  subdf <- results_pred$FI_allfeat$summary[, c('feature', 'value','sign')]
  # replace NA values by 0
  subdf$value[is.nan(subdf$value)] <- 0
  # round feature importance value
  subdf$value = round(subdf$value, digits = 3)
  # get the class vector of results
  subdf$sign = ifelse(subdf$sign == -1, ilevels[[1]], ilevels[[2]])
  # get the class vector of results
  fbm.species = results_pred$FI_fmbFeats$summary$feature
  # add a column to indicate if a species is an indicator or not
  subdf$IsIndsp = as.integer(subdf$feature %in% fbm.species)
  colnames(subdf) <- c("feature","featureImportance", "class","IsIndSp")
  #compute prevalence of each species
  subdf$prevalence <- (colSums(X> 0)/nrow(X))* 100
  results_pred$pred_out_fbm <- subdf
  
  save_results <- results_pred
  return(save_results)
}

# compute permanova analysis
permanova_analysis_predomics <- function(sample.info, sample.X, sample.class, results_pred, hab_type, ilevels, algorithm_language, data_type) {
  
  # Initialize a list to store results
  adonis_pred <- list()
  
  if (data_type == "maxN") {
    # Get features from FBM
    fbm.species <- results_pred$FI_fmbFeats$summary$feature
    
    # Subset the abundance table to these species
    subdf.features.df <- sample.X[, fbm.species]
    
    # Exclude samples with no species
    subdf.features.df <- subdf.features.df[rowSums(subdf.features.df) > 0, ]
    
    # Get class according to rows of subdf.features.df
    subdf.features.df.class <- sample.info[rownames(subdf.features.df),]
    subdf.features.df.class <- ifelse(subdf.features.df.class[[hab_type]] == ilevels[[1]], -1,
                                      ifelse(subdf.features.df.class[[hab_type]] == ilevels[[2]], 1, NA))
    
    # Check for any NA classes
    na_indices <- which(is.na(subdf.features.df.class))
    if (length(na_indices) > 0) {
      warning(paste("NA classes found at indices:", paste(na_indices, collapse = ", ")))
    }
    
    # Compute Bray-Curtis distances
    subdf.features.df.bray <- vegdist(subdf.features.df, method = "bray")
    
    # Do the PERMANOVA with subset species
    subdf.meta <- data.frame(sample = rownames(subdf.features.df), class = subdf.features.df.class)
    set.seed(100)
    subdf.features.df.bray.adonis <- adonis2(subdf.features.df.bray ~ class, data = subdf.meta)
    
    # Extract results
    subdf.features.df.bray.adonis <- data.frame(subdf.features.df.bray.adonis)["class", , drop = FALSE]
    subdf.features.df.bray.adonis$comparison <- paste(ilevels, collapse = "_")
    subdf.features.df.bray.adonis$data <- "maxN"
    subdf.features.df.bray.adonis$source <- paste0(algorithm_language, "Species")
    subdf.features.df.bray.adonis$features <- ncol(subdf.features.df)
    
    # Compute Bray-Curtis distances for the entire community
    ilevels.df.bray <- vegdist(sample.X, method = "bray")
    alldf.meta <- data.frame(sample = rownames(sample.X), class = sample.class)
    set.seed(100)
    ilevels.df.bray.adonis <- adonis2(ilevels.df.bray ~ class, data = alldf.meta)
    
    # Extract results
    ilevels.df.bray.adonis <- data.frame(ilevels.df.bray.adonis)["class", , drop = FALSE]
    ilevels.df.bray.adonis$comparison <- paste(ilevels, collapse = "_")
    ilevels.df.bray.adonis$data <- "maxN"
    ilevels.df.bray.adonis$source <- paste("allSpecies", algorithm_language, sep = '_')
    ilevels.df.bray.adonis$features <- ncol(sample.X)
    
    # Combine results
    adonis_pred <- rbind(ilevels.df.bray.adonis, subdf.features.df.bray.adonis)
    
  } else if (data_type == "pres/abs") {
    # Similar steps for presence/absence data
    fbm.species <- results_pred$FI_fmbFeats$summary$feature
    subdf.features.df <- sample.X[, fbm.species]
    
    # Exclude samples with no species
    subdf.features.df <- subdf.features.df[rowSums(subdf.features.df) > 0, ]
    
    # Get class according to rows of subdf.features.df
    subdf.features.df.class <- sample.info[rownames(subdf.features.df),]
    subdf.features.df.class <- ifelse(subdf.features.df.class[[hab_type]] == ilevels[[1]], -1,
                                      ifelse(subdf.features.df.class[[hab_type]] == ilevels[[2]], 1, NA))
    
    # Check for any NA classes
    na_indices <- which(is.na(subdf.features.df.class))
    if (length(na_indices) > 0) {
      warning(paste("NA classes found at indices:", paste(na_indices, collapse = ", ")))
    }
    
    # Compute Jaccard distances
    subdf.features.df.jaccard <- vegdist(subdf.features.df, method = "jaccard", binary = TRUE)
    
    # Do the PERMANOVA with subset species
    subdf.meta <- data.frame(sample = rownames(subdf.features.df), class = subdf.features.df.class)
    set.seed(100)
    subdf.features.df.jaccard.adonis <- adonis2(subdf.features.df.jaccard ~ class, data = subdf.meta)
    
    
    # Extract results
    subdf.features.df.jaccard.adonis <- data.frame(subdf.features.df.jaccard.adonis)["class", , drop = FALSE]
    subdf.features.df.jaccard.adonis$comparison <- paste(ilevels, collapse = "_")
    subdf.features.df.jaccard.adonis$data <- "pres/abs"
    subdf.features.df.jaccard.adonis$source <- paste0(algorithm_language, "Species")
    subdf.features.df.jaccard.adonis$features <- ncol(subdf.features.df)
    
    # Compute Jaccard distances for the entire community
    ilevels.df.jaccard <- vegdist(sample.X, method = "jaccard", binary = TRUE)
    alldf.meta <- data.frame(sample = rownames(sample.X), class = sample.class) 
    set.seed(100)
    ilevels.df.jaccard.adonis <- adonis2(ilevels.df.jaccard ~ class, data = alldf.meta)
    
    # Extract results
    ilevels.df.jaccard.adonis <- data.frame(ilevels.df.jaccard.adonis)["class", , drop = FALSE]
    ilevels.df.jaccard.adonis$comparison <- paste(ilevels, collapse = "_")
    ilevels.df.jaccard.adonis$data <- "pres/abs"
    ilevels.df.jaccard.adonis$source <- paste("allSpecies", algorithm_language, sep = '_')
    ilevels.df.jaccard.adonis$features <- ncol(sample.X)
    
    # Combine results
    adonis_pred <- rbind(ilevels.df.jaccard.adonis, subdf.features.df.jaccard.adonis)
  } else {
    stop("Invalid data_type specified. Use 'maxN' or 'pres/abs'.")
  }
  
  return(adonis_pred)
}


# for (seed in 1:maxseed) {
#   
#   
# }

print("##################### starting predomics analyses on data in abundance ###############################")

####################
## Analyses on abundance data (MaxN)
####################

predout.maxn <- list()
adonis_pred.maxn <- list()

# compute predomics and PERMANOVA analysis on BAY,LAG and BAR
for(i in 1:length(comp))
{
  # print the seed value
  # print(paste0("You are spliting data with seed = ",seed, " / ", maxseed)) 
  
  # print comparison
  print(comp[[i]])
  ##Get the levels to compare
  ilevels <- comp[[i]]
  ##Get the data limited to the levels to compare
  ilevels.df <- acc[acc$hab %in% ilevels,]
  
  # split data into train and test
  # set.seed(seed)  # For reproducibility
  # split_ind <- caTools::sample.split(Y = ilevels.df$hab, SplitRatio = 0.8)
  # table(split_ind)
  # ## train data
  # ilevels.df.train <- ilevels.df[split_ind,]
  # Train data on ABORE replicas
  ilevels.df.train <- ilevels.df[ilevels.df$Transect=="ABORE",]
  
  # prepare data for analysis
  data.train <- prepare_data_for_analysis(ilevels.df.train, habitat_type = "hab", ilevels[[1]], ilevels[[2]])
  # get elements returned
  ilevels.df.train = data.train$sample.info
  ilevels.df.train.X= data.train$sample.X
  ilevels.df.train.class= data.train$sample.class
  
  # ## test data
  # ilevels.df.test <- ilevels.df[!split_ind,]
  ## test data on MBERE replicas
  ilevels.df.test <- ilevels.df[ilevels.df$Transect=="MBERE",]
  # prepare data for analysis
  data.test <- prepare_data_for_analysis(ilevels.df.test, habitat_type = "hab", ilevels[[1]], ilevels[[2]])
  # get elements returned
  ilevels.df.test = data.test$sample.info
  ilevels.df.test.X= data.test$sample.X
  ilevels.df.test.class= data.test$sample.class
  
  ##Do the predomics analyses
  predomics_res_list= predomics_analyis(ilevels.df.train.X, ilevels.df.train.class, ilevels.df.test.X, ilevels.df.test.class, algorithm_language)
  
  #Merge the indval and pvalues output in single table + save
  if(identical(names(predomics_res_list$Comp_data$X_train), as.character(predomics_res_list$FI_allfeat$summary$feature)))
  {
    # save predomics results
    predout.maxn[[paste(ilevels, collapse = "_")]] <- save_pred_results(predomics_res_list, ilevels.df.train.X, ilevels)
    ##Analyse permanova; get the species to include
    adonis_pred.maxn[[paste(ilevels, collapse = "_")]] <- permanova_analysis_predomics(ilevels.df.train,ilevels.df.train.X, ilevels.df.train.class, predomics_res_list, hab_type = "hab", ilevels, algorithm_language, data_type = "maxN")
  }
}

#Build a comparison vector for INSHORE_OFFSHORE  
ilevels.df <- acc

# print the seed value
#print(paste0("You are spliting data with seed = ",seed, " / ", maxseed)) 

# print comparison
print(ilevels_hab)

# split data into train and test
# set.seed(seed)  # For reproducibility
# split_ind_inoff <- caTools::sample.split(Y = ilevels.df$hab_inoff, SplitRatio = 0.8)
# 
# ## train data
# ilevels.df.train <- ilevels.df[split_ind_inoff,]
# Train data on ABORE replicas
ilevels.df.train <- ilevels.df[ilevels.df$Transect=="ABORE",]

# prepare data for analysis
data.train <- prepare_data_for_analysis(ilevels.df.train, habitat_type = "hab_inoff", ilevels_hab[[1]], ilevels_hab[[2]])
# get elements returned
ilevels.df.train= data.train$sample.info
ilevels.df.train.X= data.train$sample.X
ilevels.df.train.class= data.train$sample.class

# ## test data
# ilevels.df.test <- ilevels.df[!split_ind_inoff,]
## test data on MBERE replicas
ilevels.df.test <- ilevels.df[ilevels.df$Transect=="MBERE",]

# prepare data for analysis
data.test <- prepare_data_for_analysis(ilevels.df.test, habitat_type = "hab_inoff", ilevels_hab[[1]], ilevels_hab[[2]])
# get elements returned
ilevels.df.test= data.test$sample.info
ilevels.df.test.X= data.test$sample.X
ilevels.df.test.class= data.test$sample.class

##Do the predomics analyses
predomics_res_list= predomics_analyis(ilevels.df.train.X, ilevels.df.train.class, ilevels.df.test.X, ilevels.df.test.class, algorithm_language)

#Merge the featImport and pvalues output in single table + save
if(identical(names(predomics_res_list$Comp_data$X_train), as.character(predomics_res_list$FI_allfeat$summary$feature)))
{
  # save predomics results
  predout.maxn[[paste(ilevels_hab, collapse = "_")]] <- save_pred_results(predomics_res_list, ilevels.df.train.X, ilevels_hab)
  ##Analyse permanova; get the species to include
  adonis_pred.maxn[[paste(ilevels_hab, collapse = "_")]] <- permanova_analysis_predomics(ilevels.df.train,ilevels.df.train.X, ilevels.df.train.class, predomics_res_list, hab_type = "hab_inoff", ilevels_hab, algorithm_language, data_type = "maxN")
}

##Get the pred_out table
predout.maxn.sub <- lapply(predout.maxn, function(x){x[["pred_out_fbm"]]})
for(i in names(predout.maxn.sub))
{
  predout.maxn.sub[[i]][,"comparison"] <- i
}
predout.maxn.sub <- do.call("rbind", predout.maxn.sub)
predout.maxn.sub$source <- algorithm_language
predout.maxn.sub$data <- "maxN"
#predout.maxn.sub$seed <- seed

print("##################### starting predomics analyses on data in presence/absence ###############################")

####################
## Analyses on presence/absence data
####################

predout.bin <- list()
adonis_pred.bin <- list()

# compute predomics and PERMANOVA analysis on BAY,LAG and BAR

for(i in 1:length(comp))
{
  # print the seed value
  #print(paste0("You are spliting data with seed = ",seed, " / ", maxseed)) 
  
  # print comparison
  print(comp[[i]])
  ##Get the levels to compare
  ilevels <- comp[[i]]
  ##Get the data limited to the levels to compare
  ilevels.df <- acc[acc$hab %in% ilevels,]
  
  # split data into train and test
  # set.seed(seed)  # For reproducibility
  # split_ind <- caTools::sample.split(Y = ilevels.df$hab, SplitRatio = 0.8)
  # table(split_ind)
  
  # ## train data
  # ilevels.df.train <- ilevels.df[split_ind,]
  # Train data on ABORE replicas
  ilevels.df.train <- ilevels.df[ilevels.df$Transect=="ABORE",]
  
  # prepare data for analysis
  data.train <- prepare_data_for_analysis(ilevels.df.train, habitat_type = "hab", ilevels[[1]], ilevels[[2]])
  # get elements returned
  ilevels.df.train = data.train$sample.info
  ilevels.df.train.X= data.train$sample.X
  ##Transform on binary
  ilevels.df.train.X <- as.data.frame(apply(ilevels.df.train.X, 2, function(x){ifelse(x==0,0,1)}))
  ilevels.df.train.class= data.train$sample.class
  
  # ## test data
  # ilevels.df.test <- ilevels.df[!split_ind,]
  ## test data on MBERE replicas
  ilevels.df.test <- ilevels.df[ilevels.df$Transect=="MBERE",]
  
  # prepare data for analysis
  data.test <- prepare_data_for_analysis(ilevels.df.test, habitat_type = "hab", ilevels[[1]], ilevels[[2]])
  # get elements returned
  ilevels.df.test = data.test$sample.info
  ilevels.df.test.X= data.test$sample.X
  ##Transform on binary
  ilevels.df.test.X <- as.data.frame(apply(ilevels.df.test.X, 2, function(x){ifelse(x==0,0,1)}))
  ilevels.df.test.class= data.test$sample.class
  
  ##Do the predomics analyses
  predomics_res_list= predomics_analyis(ilevels.df.train.X, ilevels.df.train.class, ilevels.df.test.X, ilevels.df.test.class, algorithm_language)
  #Merge the indval and pvalues output in single table + save
  if(identical(names(predomics_res_list$Comp_data$X_train), as.character(predomics_res_list$FI_allfeat$summary$feature)))
  {
    # save predomics results
    predout.bin[[paste(ilevels, collapse = "_")]] <- save_pred_results(predomics_res_list, ilevels.df.train.X, ilevels)
    ##Analyse permanova; get the species to include
    adonis_pred.bin[[paste(ilevels, collapse = "_")]] <- permanova_analysis_predomics(ilevels.df.train,ilevels.df.train.X, ilevels.df.train.class, predomics_res_list, hab_type = "hab", ilevels, algorithm_language, data_type = "pres/abs")
    
  }
}

#Build a comparison vector for INSHORE_OFFSHORE  
ilevels.df <- acc

# print the seed value
# print(paste0("You are spliting data with seed = ",seed, " / ", maxseed)) 

# print comparison
print(ilevels_hab)

# split data into train and test
# set.seed(seed)  # For reproducibility
# split_ind_inoff <- caTools::sample.split(Y = ilevels.df$hab_inoff, SplitRatio = 0.8)

# ## train data
# ilevels.df.train <- ilevels.df[split_ind_inoff,]
# Train data on ABORE replicas
ilevels.df.train <- ilevels.df[ilevels.df$Transect=="ABORE",]

# prepare data for analysis
data.train <- prepare_data_for_analysis(ilevels.df.train, habitat_type = "hab_inoff", ilevels_hab[[1]], ilevels_hab[[2]])
# get elements returned
ilevels.df.train= data.train$sample.info
ilevels.df.train.X= data.train$sample.X
##Transform on binary
ilevels.df.train.X <- as.data.frame(apply(ilevels.df.train.X, 2, function(x){ifelse(x==0,0,1)}))
ilevels.df.train.class= data.train$sample.class

# ## test data
# ilevels.df.test <- ilevels.df[!split_ind_inoff,]
## test data on MBERE replicas
ilevels.df.test <- ilevels.df[ilevels.df$Transect=="MBERE",]

# prepare data for analysis
data.test <- prepare_data_for_analysis(ilevels.df.test, habitat_type = "hab_inoff", ilevels_hab[[1]], ilevels_hab[[2]])
# get elements returned
ilevels.df.test= data.test$sample.info
ilevels.df.test.X= data.test$sample.X
##Transform on binary
ilevels.df.test.X <- as.data.frame(apply(ilevels.df.test.X, 2, function(x){ifelse(x==0,0,1)}))
ilevels.df.test.class= data.test$sample.class

##Do the predomics analyses
predomics_res_list= predomics_analyis(ilevels.df.train.X, ilevels.df.train.class, ilevels.df.test.X, ilevels.df.test.class, algorithm_language)

#Merge the featImport and pvalues output in single table + save
if(identical(names(predomics_res_list$Comp_data$X_train), as.character(predomics_res_list$FI_allfeat$summary$feature)))
{
  # save predomics results
  predout.bin[[paste(ilevels_hab, collapse = "_")]] <- save_pred_results(predomics_res_list, ilevels.df.train.X, ilevels_hab)
  ##Analyse permanova; get the species to include
  adonis_pred.bin[[paste(ilevels_hab, collapse = "_")]] <- permanova_analysis_predomics(ilevels.df.train,ilevels.df.train.X, ilevels.df.train.class, predomics_res_list, hab_type = "hab_inoff", ilevels_hab, algorithm_language, data_type = "pres/abs")
  
}

##Get the pred_out table
predout.bin.sub <- lapply(predout.bin, function(x){x[["pred_out_fbm"]]})
for(i in names(predout.bin.sub))
{
  predout.bin.sub[[i]][,"comparison"] <- i
}
predout.bin.sub <- do.call("rbind", predout.bin.sub)
predout.bin.sub$source <- algorithm_language
predout.bin.sub$data <- "pres/abs"
#predout.bin.sub$seed <-seed

# Define the save paths
save_path_pred <- paste0('/data/projects/aime/analyses/bruvs/results_replicas_analysis/', algorithm_language, "_output_data/")

# Create the directory if it doesn't exist
if (!dir.exists(save_path_pred)) {
  dir.create(save_path_pred, recursive = TRUE)
}

# Define the full file path for the .Rda file
file_path <- paste0(save_path_pred, algorithm_language, "_Predomics_all_comparison_replicas_prev_", prevalence_rate, ".Rda")

# Check if the file already exists before saving
if (!file.exists(file_path)) {
  save(
    adonis_pred.bin,
    adonis_pred.maxn,
    predout.bin,
    predout.maxn,
    predout.bin.sub,
    predout.maxn.sub,
    file = file_path
  )
  message("File saved: ", file_path)
} else {
  message("File already exists, skipping save: ", file_path)
}


