# To run this script, use this command: Rscript codePredomics_all_analyses.R terinter

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

# define arguments for the script
args <- commandArgs(trailingOnly = TRUE)
algorithm_language <- args[1]
prevalence_rate <- args[2]
max_seed <- args[3]

dataset <- "/data/projects/aime/data/metabarcoding/AllDataexportefixed.txt"
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
matablepivot=matablepivot[c("Site","Station","G_SP","MaxN")]
head(matablepivot)

# Species community matrix
abund = reshape2::dcast(matablepivot, Station + Site ~ G_SP, value.var = "MaxN")
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

## Running analyses from a range of seed

# function to get sample data by prevalence
get_sample_by_prevalence <- function(abundance_matrix, prevalence_rate) {
  
  print(paste("You set the prevalence rate to:", prevalence_rate, sep = " "))
  # scale the prevalence rate in [0,1]
  prevalence_rate <- prevalence_rate/100
  # Check if prevalence rate is between 0 and 1
  if(prevalence_rate==0){
    print("You are using the overall data without any filtering")
  }
  if(prevalence_rate==1){
    print("You are using species that appears in all samples")
  }
  if (prevalence_rate < 0 || prevalence_rate > 1) {
    stop("Prevalence rate must be between 0 and 1.")
  }
  
  # Get the number of samples (columns)
  num_samples <- nrow(abundance_matrix)
  
  # Calculate the minimum number of non-zero occurrences required
  min_occurrences <- ceiling(prevalence_rate * num_samples)
  
  # Select rows that meet the prevalence threshold
  valid_features <- colSums(abundance_matrix > 0) >= min_occurrences
  
  # Return the subset of the abundance matrix
  sampled_matrix <- abundance_matrix[ ,valid_features, drop = FALSE]
  
  return(sampled_matrix)
}

for (seed in 1:max_seed) {

  ####################
  ## Analyses on abundance data (MaxN)
  ####################
  
  predout.maxn <- list()
  adonis_pred.maxn <- list()
  for(i in 1:length(comp))
  {
    print(i)
    ##Get the levels to compare
    ilevels <- comp[[i]]
    ##Get the data limited to the levels to compare
    ilevels.df <- acc[acc$hab %in% ilevels,]
    
    # 
    ilevels.df <- get_sample_by_prevalence(ilevels.df,)
    
    # split data into train and test
    set.seed(seed)  # For reproducibility
    split_ind <- caTools::sample.split(Y = ilevels.df$hab, SplitRatio = 0.8)
    table(split_ind)
    ## train data
    ilevels.df.train <- ilevels.df[split_ind,]
    ilevels.df.train.class <- ifelse(ilevels.df.train$hab==ilevels[[1]], -1,
                               ifelse(ilevels.df.train$hab==ilevels[[2]], 1, NA))
    which(is.na(ilevels.df.train.class))
    table(ilevels.df.train$hab, ilevels.df.train.class)
    
    ## test data
    ilevels.df.test <- ilevels.df[!split_ind,]
    ilevels.df.test.class <- ifelse(ilevels.df.test$hab==ilevels[[1]], -1,
                                     ifelse(ilevels.df.test$hab==ilevels[[2]], 1, NA))
    which(is.na(ilevels.df.test.class))
    table(ilevels.df.test$hab, ilevels.df.test.class)
    
    ##Exclude unnecesary columns (keep abundance table only samples x species)
    ilevels.df.train <- ilevels.df.train[,-match(c("Station","Site","hab","hab_inoff"), colnames(acc))]
    ilevels.df.test <- ilevels.df.test[,-match(c("Station","Site","hab","hab_inoff"), colnames(acc))]
    
    ##Do the predomics analyses
   
    # set the learner
    ifelse(algorithm_language == "terinter",
              clf <- terBeam(language = 'terinter',
                             nCores = 1,
                             objective = "auc",
                             seed = 20,
                             plot = TRUE),
           ifelse(algorithm_language == "bininter",
              clf <- terBeam(language = 'bininter',
                                         nCores = 1,
                                         objective = "auc",
                                         seed = 20,
                                         plot = TRUE),
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
    predomics_res_list <- list()
    
    # store data for training and testing
    predomics_res_list$Comp_data <- list(X_train= ilevels.df.train, y_train= ilevels.df.train.class, 
                                         X_test= ilevels.df.test, y_test= ilevels.df.test.class)
    
    # save clf and classification results in the master list
    predomics_res_list$fit <- res_clf
    predomics_res_list$clf <- clf
    
    # save digest of results
    res_clf.dig <- digest(obj = res_clf, penalty = 0.75/100, plot = TRUE)
    predomics_res_list$digest <- res_clf.dig
    
    ### Testing the model in the testing dataset (20%)
    
    best.model.test <- evaluateModel(mod = res_clf.dig$best$model, X = ilevels.df.test, y = ilevels.df.test.class, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE, mode = "test")
    predomics_res_list$test_eval <- best.model.test
    
    #### Family of Best Models (FBM) ####
    
    # get the population of models scrambled by model size
    pop <- modelCollectionToPopulation(res_clf$classifier$models)
    pop.df <- populationToDataFrame(pop)
    predomics_res_list$model_pop <- pop.df
    
    # select the best
    fbm <- selectBestPopulation(pop)
    fbm.df <- populationToDataFrame(fbm)
    predomics_res_list$fbm <- fbm.df
    
    #### Family of Best Models (FBM) - Feature Annotation ####
    
    fa <- makeFeatureAnnot(pop = fbm,
                           X = t(ilevels.df.train),
                           y = ilevels.df.train.class,
                           clf = clf)
    
    predomics_res_list$featureAnnotFBM <- fa
    
    #### Feature importance ####
    
    feat1.import <- mergeMeltImportanceCV(list.results = list(terBeam = res_clf),
                                          filter.cv.prev = 0,
                                          min.kfold.nb = FALSE,
                                          learner.grep.pattern = "*",
                                          #nb.top.features = 50,
                                          feature.selection = rownames(fa$pop.noz),
                                          scaled.importance = TRUE,
                                          make.plot = TRUE,
                                          cv.prevalence = FALSE)
    
    
    feat2.import <- mergeMeltImportanceCV(list.results = list(terBeam = res_clf),
                                          filter.cv.prev = 0,
                                          min.kfold.nb = FALSE,
                                          learner.grep.pattern = "*",
                                          nb.top.features = 148,
                                          #feature.selection = rownames(fa$pop.noz),
                                          scaled.importance = TRUE,
                                          make.plot = TRUE,
                                          cv.prevalence = FALSE)
    
    
    predomics_res_list$FI_fmbFeats <- feat1.import
    predomics_res_list$FI_allfeat <- feat2.import
    
    
    #Merge the indval and pvalues output in single table + save
    if(identical(names(predomics_res_list$Comp_data$X_train), as.character(predomics_res_list$FI_allfeat$summary$feature)))
    {
      # get all species and specify if there are from FBM or not with IsIndSp
      subdf <- predomics_res_list$FI_allfeat$summary[, c('feature', 'value','sign')]
      
      subdf$value[is.nan(subdf$value)] <- 0
      subdf$value = round(subdf$value, digits = 3)
      subdf$sign = ifelse(subdf$sign == -1, ilevels[[1]], ilevels[[2]])
      fbm.species = predomics_res_list$FI_fmbFeats$summary$feature
      subdf$IsIndsp = as.integer(subdf$feature %in% fbm.species)
      colnames(subdf) <- c("feature","featureImportance", "class","IsIndSp")
      #compute prevalence of each species
      subdf$prevalence <- (colSums(ilevels.df.train > 0)/nrow(ilevels.df.train))* 100
      predomics_res_list$pred_out_fbm <- subdf
      predout.maxn[[paste(ilevels, collapse = "_")]] <- predomics_res_list
      ##Analyse permanova; get the species to include (from FBM)
      ##Exclude species absents in all subset samples
      ilevels.df.train <- ilevels.df.train[,colSums(ilevels.df.train>0)>0]
      ilevels.df.test <- ilevels.df.test[,colSums(ilevels.df.test>0)>0]
      #Subset the abundance table to these species
      subdf.features.df <- ilevels.df.train[, fbm.species]
      ## Exclude sample with no species
      subdf.features.df <- subdf.features.df[rowSums(subdf.features.df) != 0, ]
      # get class according to rows of subdf.features.df
      subdf.features.df.class <- ilevels.df[rownames(subdf.features.df),]
      subdf.features.df.class <- ifelse(subdf.features.df.class$hab==ilevels[[1]], -1,
                                        ifelse(subdf.features.df.class$hab==ilevels[[2]], 1, NA))
      which(is.na(subdf.features.df.class))
      #Compute bray curtis distances
      subdf.features.df.bray <- vegdist(subdf.features.df, method="bray")
      #Do the permanova with subset species
      subdf.meta <- data.frame(sample=rownames(subdf.features.df), class=subdf.features.df.class)
      set.seed(100)
      subdf.features.df.bray.adonis <- adonis2(subdf.features.df.bray~class, data=subdf.meta)
      subdf.features.df.bray.adonis <- data.frame(subdf.features.df.bray.adonis)
      subdf.features.df.bray.adonis <- subdf.features.df.bray.adonis["class",,drop=FALSE]
      subdf.features.df.bray.adonis$comparison <- paste(ilevels, collapse = "_")
      subdf.features.df.bray.adonis$data <- "maxN"
      subdf.features.df.bray.adonis$source <- paste0(algorithm_language,"Species")
      subdf.features.df.bray.adonis$features <- length(subdf$feature)
      
      #Compute bray curtis distances; all community
      ilevels.df.bray <- vegdist(ilevels.df.train, method="bray")
      alldf.meta <- data.frame(sample=rownames(ilevels.df.train), class=ilevels.df.train.class)
      set.seed(100)
      ilevels.df.bray.adonis <- adonis2(ilevels.df.bray~class, data=alldf.meta)
      ilevels.df.bray.adonis <- data.frame(ilevels.df.bray.adonis)
      ilevels.df.bray.adonis <- ilevels.df.bray.adonis["class",,drop=FALSE]
      ilevels.df.bray.adonis$comparison <- paste(ilevels, collapse = "_")
      ilevels.df.bray.adonis$data <- "maxN"
      ilevels.df.bray.adonis$source <- paste("allSpecies", algorithm_language, sep = '_')
      ilevels.df.bray.adonis$features <- ncol(ilevels.df.train)
      #put all together + save
      ilevels.df.bray.adonis.all <- rbind(ilevels.df.bray.adonis, subdf.features.df.bray.adonis)
      adonis_pred.maxn[[paste(ilevels, collapse = "_")]] <- ilevels.df.bray.adonis.all
    }
  }
  
  #Build a comparison vector for INSHORE_OFFSHORE  
  ilevels.df <- acc
  
  # split data into train and test
  set.seed(seed)  # For reproducibility
  split_ind_inoff <- caTools::sample.split(Y = ilevels.df$hab_inoff, SplitRatio = 0.8)
  
  ## train data
  ilevels.df.train <- ilevels.df[split_ind_inoff,]
  ilevels.df.train.class <- ifelse(ilevels.df.train$hab_inoff=="INSHORE", -1,
                            ifelse(ilevels.df.train$hab_inoff=="OFFSHORE", 1, NA))
  
  which(is.na(ilevels.df.train.class))
  table(ilevels.df.train$hab_inoff, ilevels.df.train.class)
  
  ## test data
  ilevels.df.test <- ilevels.df[!split_ind_inoff,]
  ilevels.df.test.class <- ifelse(ilevels.df.test$hab_inoff=="INSHORE", -1,
                            ifelse(ilevels.df.test$hab_inoff=="OFFSHORE", 1, NA))
  
  which(is.na(ilevels.df.test.class))
  table(ilevels.df.test$hab_inoff, ilevels.df.test.class)
  
  ##Exclude unnecesary columns (keep abundance table only samples x species)
  ilevels.df.train <- ilevels.df.train[,-match(c("Station","Site","hab","hab_inoff"), colnames(acc))]
  ilevels.df.test <- ilevels.df.test[,-match(c("Station","Site","hab","hab_inoff"), colnames(acc))]
  
  ##Do the predomics analyses
  
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
  predomics_res_list <- list()
  
  # store data for training and testing
  predomics_res_list$Comp_data <- list(X_train= ilevels.df.train, y_train= ilevels.df.train.class, 
                                       X_test= ilevels.df.test, y_test= ilevels.df.test.class)
  
  # save clf and classification results in the master list
  predomics_res_list$fit <- res_clf
  predomics_res_list$clf <- clf
  
  # save digest of results
  res_clf.dig <- digest(obj = res_clf, penalty = 0.75/100, plot = TRUE)
  predomics_res_list$digest <- res_clf.dig
  
  ### Testing the model in the testing dataset (20%)
  
  best.model.test <- evaluateModel(mod = res_clf.dig$best$model, X = t(ilevels.df.test), y = ilevels.df.test.class, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE, mode = "test")
  predomics_res_list$test_eval <- best.model.test
  
  #### Family of Best Models (FBM) ####
  
  # get the population of models scrambled by model size
  pop <- modelCollectionToPopulation(res_clf$classifier$models)
  pop.df <- populationToDataFrame(pop)
  predomics_res_list$model_pop <- pop.df
  
  # select the best
  fbm <- selectBestPopulation(pop)
  fbm.df <- populationToDataFrame(fbm)
  predomics_res_list$fbm <- fbm.df
  
  #### Family of Best Models (FBM) - Feature Annotation ####
  
  fa <- makeFeatureAnnot(pop = fbm,
                         X = t(ilevels.df.train),
                         y = ilevels.df.train.class,
                         clf = clf)
  
  predomics_res_list$featureAnnotFBM <- fa
  
  #### Feature importance ####
  
  feat1.import <- mergeMeltImportanceCV(list.results = list(terBeam = res_clf),
                                        filter.cv.prev = 0,
                                        min.kfold.nb = FALSE,
                                        learner.grep.pattern = "*",
                                        #nb.top.features = 50,
                                        feature.selection = rownames(fa$pop.noz),
                                        scaled.importance = TRUE,
                                        make.plot = TRUE,
                                        cv.prevalence = FALSE)
  
  
  feat2.import <- mergeMeltImportanceCV(list.results = list(terBeam = res_clf),
                                        filter.cv.prev = 0,
                                        min.kfold.nb = FALSE,
                                        learner.grep.pattern = "*",
                                        nb.top.features = 148,
                                        #feature.selection = rownames(fa$pop.noz),
                                        scaled.importance = TRUE,
                                        make.plot = TRUE,
                                        cv.prevalence = FALSE)
  
  
  predomics_res_list$FI_fmbFeats <- feat1.import
  predomics_res_list$FI_allfeat <- feat2.import
  
  
  #Merge the featImport and pvalues output in single table + save
  if(identical(names(predomics_res_list$Comp_data$X_train), as.character(predomics_res_list$FI_allfeat$summary$feature)))
  {
    # get all species and specify if there are from FBM or not with IsIndSp
    subdf <- predomics_res_list$FI_allfeat$summary[, c('feature', 'value','sign')]
    
    subdf$value[is.nan(subdf$value)] <- 0
    subdf$value = round(subdf$value, digits = 3)
    subdf$sign = ifelse(subdf$sign == -1, "INSHORE", "OFFSHORE")
    fbm.species = predomics_res_list$FI_fmbFeats$summary$feature
    subdf$IsIndsp = as.integer(subdf$feature %in% fbm.species)
    colnames(subdf) <- c("feature","featureImportance", "class","IsIndSp")
    #compute prevalence of each species
    subdf$prevalence <- (colSums(ilevels.df.train > 0)/nrow(ilevels.df.train))* 100
    predomics_res_list$pred_out_fbm <- subdf
    predout.maxn[["INSHORE_OFFSHORE"]] <- predomics_res_list
    ##Analyse permanova; get the species to include (from FBM)
    ##Exclude species absents in all subset samples
    ilevels.df.train <- ilevels.df.train[,colSums(ilevels.df.train>0)>0]
    ilevels.df.test <- ilevels.df.test[,colSums(ilevels.df.test>0)>0]
    #Subset the abundance table to these species
    subdf.features.df <- ilevels.df.train[, fbm.species]
    ## Exclude sample with no species
    subdf.features.df <- subdf.features.df[rowSums(subdf.features.df) != 0, ]
    # get class according to rows of subdf.features.df
    subdf.features.df.class <- ilevels.df[rownames(subdf.features.df),]
    subdf.features.df.class <- ifelse(subdf.features.df.class$hab_inoff=='INSHORE', -1,
                                      ifelse(subdf.features.df.class$hab_inoff=='OFFSHORE', 1, NA))
    which(is.na(subdf.features.df.class))
    #Compute bray curtis distances
    subdf.features.df.bray <- vegdist(subdf.features.df, method="bray")
    #Do the permanova with subset species
    subdf.meta <- data.frame(sample=rownames(subdf.features.df), class=subdf.features.df.class)
    set.seed(100)
    subdf.features.df.bray.adonis <- adonis2(subdf.features.df.bray~class, data=subdf.meta)
    subdf.features.df.bray.adonis <- data.frame(subdf.features.df.bray.adonis)
    subdf.features.df.bray.adonis <- subdf.features.df.bray.adonis["class",,drop=FALSE]
    subdf.features.df.bray.adonis$comparison <- "INSHORE_OFFSHORE"
    subdf.features.df.bray.adonis$data <- "maxN"
    subdf.features.df.bray.adonis$source <- paste0(algorithm_language,"Species")
    subdf.features.df.bray.adonis$features <- length(subdf$feature)
    
    #Compute bray curtis distances; all community
    ilevels.df.bray <- vegdist(ilevels.df.train, method="bray")
    alldf.meta <- data.frame(sample=rownames(ilevels.df.train), class=ilevels.df.train.class)
    set.seed(100)
    ilevels.df.bray.adonis <- adonis2(ilevels.df.bray~class, data=alldf.meta)
    ilevels.df.bray.adonis <- data.frame(ilevels.df.bray.adonis)
    ilevels.df.bray.adonis <- ilevels.df.bray.adonis["class",,drop=FALSE]
    ilevels.df.bray.adonis$comparison <- "INSHORE_OFFSHORE"
    ilevels.df.bray.adonis$data <- "maxN"
    ilevels.df.bray.adonis$source <- paste("allSpecies", algorithm_language, sep = '_')
    ilevels.df.bray.adonis$features <- ncol(ilevels.df.train)
    #put all together + save
    ilevels.df.bray.adonis.all <- rbind(ilevels.df.bray.adonis, subdf.features.df.bray.adonis)
    adonis_pred.maxn[["INSHORE_OFFSHORE"]] <- ilevels.df.bray.adonis.all
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
  
  
  ####################
  ## Analyses on presence/absence data
  ####################
  
  predout.bin <- list()
  adonis_pred.bin <- list()
  
  for(i in 1:length(comp))
  {
    print(i)
    ##Get the levels to compare
    ilevels <- comp[[i]]
    ##Get the data limited to the levels to compare
    ilevels.df <- acc[acc$hab %in% ilevels,]
    
    # split data into train and test
    set.seed(seed)  # For reproducibility
    split_ind <- caTools::sample.split(Y = ilevels.df$hab, SplitRatio = 0.8)
    table(split_ind)
    ## train data
    ilevels.df.train <- ilevels.df[split_ind,]
    ilevels.df.train.class <- ifelse(ilevels.df.train$hab==ilevels[[1]], -1,
                                     ifelse(ilevels.df.train$hab==ilevels[[2]], 1, NA))
    which(is.na(ilevels.df.train.class))
    table(ilevels.df.train$hab, ilevels.df.train.class)
    
    ## test data
    ilevels.df.test <- ilevels.df[!split_ind,]
    ilevels.df.test.class <- ifelse(ilevels.df.test$hab==ilevels[[1]], -1,
                                    ifelse(ilevels.df.test$hab==ilevels[[2]], 1, NA))
    which(is.na(ilevels.df.test.class))
    table(ilevels.df.test$hab, ilevels.df.test.class)
    
    ##Exclude unnecesary columns (keep abundance table only samples x species)
    ilevels.df.train <- ilevels.df.train[,-match(c("Station","Site","hab","hab_inoff"), colnames(acc))]
    ilevels.df.test <- ilevels.df.test[,-match(c("Station","Site","hab","hab_inoff"), colnames(acc))]
    
    ##Transform on binary
    ilevels.df.train <- as.data.frame(apply(ilevels.df.train, 2, function(x){ifelse(x==0,0,1)}))
    ilevels.df.test <- as.data.frame(apply(ilevels.df.test, 2, function(x){ifelse(x==0,0,1)}))
    
    ##Do the predomics analyses
    
    # set the learner
    ifelse(algorithm_language == "terinter",
           clf <- terBeam(language = 'terinter',
                          nCores = 1,
                          objective = "auc",
                          seed = 20,
                          plot = TRUE),
           ifelse(algorithm_language == "bininter",
                  clf <- terBeam(language = 'bininter',
                                 nCores = 1,
                                 objective = "auc",
                                 seed = 20,
                                 plot = TRUE),
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
    predomics_res_list <- list()
    
    # store data for training and testing
    predomics_res_list$Comp_data <- list(X_train= ilevels.df.train, y_train= ilevels.df.train.class, 
                                         X_test= ilevels.df.test, y_test= ilevels.df.test.class)
    
    # save clf and classification results in the master list
    predomics_res_list$fit <- res_clf
    predomics_res_list$clf <- clf
    
    # save digest of results
    res_clf.dig <- digest(obj = res_clf, penalty = 0.75/100, plot = TRUE)
    predomics_res_list$digest <- res_clf.dig
    
    ### Testing the model in the testing dataset (20%)
    
    best.model.test <- evaluateModel(mod = res_clf.dig$best$model, X = t(ilevels.df.test), y = ilevels.df.test.class, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE, mode = "test")
    predomics_res_list$test_eval <- best.model.test
    
    #### Family of Best Models (FBM) ####
    
    # get the population of models scrambled by model size
    pop <- modelCollectionToPopulation(res_clf$classifier$models)
    pop.df <- populationToDataFrame(pop)
    predomics_res_list$model_pop <- pop.df
    
    # select the best
    fbm <- selectBestPopulation(pop)
    fbm.df <- populationToDataFrame(fbm)
    predomics_res_list$fbm <- fbm.df
    
    #### Family of Best Models (FBM) - Feature Annotation ####
    
    fa <- makeFeatureAnnot(pop = fbm,
                           X = t(ilevels.df.train),
                           y = ilevels.df.train.class,
                           clf = clf)
    
    predomics_res_list$featureAnnotFBM <- fa
    
    #### Feature importance ####
    
    feat1.import <- mergeMeltImportanceCV(list.results = list(terBeam = res_clf),
                                          filter.cv.prev = 0,
                                          min.kfold.nb = FALSE,
                                          learner.grep.pattern = "*",
                                          #nb.top.features = 50,
                                          feature.selection = rownames(fa$pop.noz),
                                          scaled.importance = TRUE,
                                          make.plot = TRUE,
                                          cv.prevalence = FALSE)
    
    
    feat2.import <- mergeMeltImportanceCV(list.results = list(terBeam = res_clf),
                                          filter.cv.prev = 0,
                                          min.kfold.nb = FALSE,
                                          learner.grep.pattern = "*",
                                          nb.top.features = 148,
                                          #feature.selection = rownames(fa$pop.noz),
                                          scaled.importance = TRUE,
                                          make.plot = TRUE,
                                          cv.prevalence = FALSE)
    
    
    predomics_res_list$FI_fmbFeats <- feat1.import
    predomics_res_list$FI_allfeat <- feat2.import
    
    
    #Merge the indval and pvalues output in single table + save
    if(identical(names(predomics_res_list$Comp_data$X_train), as.character(predomics_res_list$FI_allfeat$summary$feature)))
    {
      # get all species and specify if there are from FBM or not with IsIndSp
      subdf <- predomics_res_list$FI_allfeat$summary[, c('feature', 'value','sign')]
      
      subdf$value[is.nan(subdf$value)] <- 0
      subdf$value = round(subdf$value, digits = 3)
      subdf$sign = ifelse(subdf$sign == -1, ilevels[[1]], ilevels[[2]])
      fbm.species = predomics_res_list$FI_fmbFeats$summary$feature
      subdf$IsIndsp = as.integer(subdf$feature %in% fbm.species)
      colnames(subdf) <- c("feature","featureImportance", "class","IsIndSp")
      #compute prevalence of each species
      subdf$prevalence <- (colSums(ilevels.df.train > 0)/nrow(ilevels.df.train))* 100
      predomics_res_list$pred_out_fbm <- subdf
      predout.bin[[paste(ilevels, collapse = "_")]] <- predomics_res_list
      ##Analyse permanova; get the species to include (from FBM)
      ##Exclude species absents in all subset samples
      ilevels.df.train <- ilevels.df.train[,colSums(ilevels.df.train>0)>0]
      ilevels.df.test <- ilevels.df.test[,colSums(ilevels.df.test>0)>0]
      #Subset the abundance table to these species
      subdf.features.df <- ilevels.df.train[, fbm.species]
      ## Exclude sample with no species
      subdf.features.df <- subdf.features.df[rowSums(subdf.features.df) != 0, ]
      # get class according to rows of subdf.features.df
      subdf.features.df.class <- ilevels.df[rownames(subdf.features.df),]
      subdf.features.df.class <- ifelse(subdf.features.df.class$hab==ilevels[[1]], -1,
                                        ifelse(subdf.features.df.class$hab==ilevels[[2]], 1, NA))
      which(is.na(subdf.features.df.class))
      #Compute jaccard distances
      subdf.features.df.jaccard <- vegdist(subdf.features.df, method="jaccard", binary=TRUE)
      #Do the permanova with subset species
      subdf.meta <- data.frame(sample=rownames(subdf.features.df), class=subdf.features.df.class)
      set.seed(100)
      subdf.features.df.jaccard.adonis <- adonis2(subdf.features.df.jaccard~class, data=subdf.meta)
      subdf.features.df.jaccard.adonis <- data.frame(subdf.features.df.jaccard.adonis)
      subdf.features.df.jaccard.adonis <-subdf.features.df.jaccard.adonis["class",,drop=FALSE]
      subdf.features.df.jaccard.adonis$comparison <- paste(ilevels, collapse = "_")
      subdf.features.df.jaccard.adonis$data <- "bin"
      subdf.features.df.jaccard.adonis$source <- paste0(algorithm_language,"Species")
      subdf.features.df.jaccard.adonis$features <- length(subdf$feature)
      
      #Compute bray curtis distances; all community
      ilevels.df.jaccard <- vegdist(ilevels.df.train, method="jaccard", binary=TRUE)
      alldf.meta <- data.frame(sample=rownames(ilevels.df.train), class=ilevels.df.train.class) 
      set.seed(100)
      ilevels.df.jaccard.adonis <- adonis2(ilevels.df.jaccard~class, data=alldf.meta)
      ilevels.df.jaccard.adonis <- data.frame(ilevels.df.jaccard.adonis)
      ilevels.df.jaccard.adonis <- ilevels.df.jaccard.adonis["class",,drop=FALSE]
      ilevels.df.jaccard.adonis$comparison <- paste(ilevels, collapse = "_")
      ilevels.df.jaccard.adonis$data <- "bin"
      ilevels.df.jaccard.adonis$source <- paste("allSpecies", algorithm_language, sep = '_')
      ilevels.df.jaccard.adonis$features <- ncol(ilevels.df.train)
      #put all together + save
      ilevels.df.jaccard.adonis.all <- rbind(ilevels.df.jaccard.adonis,subdf.features.df.jaccard.adonis)
      adonis_pred.bin[[paste(ilevels, collapse = "_")]] <- ilevels.df.jaccard.adonis.all
    }
  }
  
  #Build a comparison vector for INSHORE_OFFSHORE  
  ilevels.df <- acc
  
  # split data into train and test
  set.seed(seed)  # For reproducibility
  split_ind_inoff <- caTools::sample.split(Y = ilevels.df$hab_inoff, SplitRatio = 0.8)
  
  ## train data
  ilevels.df.train <- ilevels.df[split_ind_inoff,]
  ilevels.df.train.class <- ifelse(ilevels.df.train$hab_inoff=="INSHORE", -1,
                                   ifelse(ilevels.df.train$hab_inoff=="OFFSHORE", 1, NA))
  
  which(is.na(ilevels.df.train.class))
  table(ilevels.df.train$hab_inoff, ilevels.df.train.class)
  
  ## test data
  ilevels.df.test <- ilevels.df[!split_ind_inoff,]
  ilevels.df.test.class <- ifelse(ilevels.df.test$hab_inoff=="INSHORE", -1,
                                  ifelse(ilevels.df.test$hab_inoff=="OFFSHORE", 1, NA))
  
  which(is.na(ilevels.df.test.class))
  table(ilevels.df.test$hab_inoff, ilevels.df.test.class)
  
  ##Exclude unnecesary columns (keep abundance table only samples x species)
  ilevels.df.train <- ilevels.df.train[,-match(c("Station","Site","hab","hab_inoff"), colnames(acc))]
  ilevels.df.test <- ilevels.df.test[,-match(c("Station","Site","hab","hab_inoff"), colnames(acc))]
  
  ##Transform on binary
  ilevels.df.train <- as.data.frame(apply(ilevels.df.train, 2, function(x){ifelse(x==0,0,1)}))
  ilevels.df.test <- as.data.frame(apply(ilevels.df.test, 2, function(x){ifelse(x==0,0,1)}))
  
  ##Do the predomics analyses
  
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
  predomics_res_list <- list()
  
  # store data for training and testing
  predomics_res_list$Comp_data <- list(X_train= ilevels.df.train, y_train= ilevels.df.train.class, 
                                       X_test= ilevels.df.test, y_test= ilevels.df.test.class)
  
  # save clf and classification results in the master list
  predomics_res_list$fit <- res_clf
  predomics_res_list$clf <- clf
  
  # save digest of results
  res_clf.dig <- digest(obj = res_clf, penalty = 0.75/100, plot = TRUE)
  predomics_res_list$digest <- res_clf.dig
  
  ### Testing the model in the testing dataset (20%)
  
  # best.model.test <- evaluateModel(mod = res_clf.dig$best$model, X = t(ilevels.df.test), y = ilevels.df.test.class, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE, mode = "test")
  predomics_res_list$test_eval <- best.model.test
  
  #### Family of Best Models (FBM) ####
  
  # get the population of models scrambled by model size
  pop <- modelCollectionToPopulation(res_clf$classifier$models)
  pop.df <- populationToDataFrame(pop)
  predomics_res_list$model_pop <- pop.df
  
  # select the best
  fbm <- selectBestPopulation(pop)
  fbm.df <- populationToDataFrame(fbm)
  predomics_res_list$fbm <- fbm.df
  
  #### Family of Best Models (FBM) - Feature Annotation ####
  
  fa <- makeFeatureAnnot(pop = fbm,
                         X = t(ilevels.df.train),
                         y = ilevels.df.train.class,
                         clf = clf)
  
  predomics_res_list$featureAnnotFBM <- fa
  
  #### Feature importance ####
  
  feat1.import <- mergeMeltImportanceCV(list.results = list(terBeam = res_clf),
                                        filter.cv.prev = 0,
                                        min.kfold.nb = FALSE,
                                        learner.grep.pattern = "*",
                                        #nb.top.features = 50,
                                        feature.selection = rownames(fa$pop.noz),
                                        scaled.importance = TRUE,
                                        make.plot = TRUE,
                                        cv.prevalence = FALSE)
  
  feat2.import <- mergeMeltImportanceCV(list.results = list(terBeam = res_clf),
                                        filter.cv.prev = 0,
                                        min.kfold.nb = FALSE,
                                        learner.grep.pattern = "*",
                                        nb.top.features = 148,
                                        #feature.selection = rownames(fa$pop.noz),
                                        scaled.importance = TRUE,
                                        make.plot = TRUE,
                                        cv.prevalence = FALSE)
  
  predomics_res_list$FI_fmbFeats <- feat1.import
  predomics_res_list$FI_allfeat <- feat2.import
  
  
  #Merge the featImport and pvalues output in single table + save
  if(identical(names(predomics_res_list$Comp_data$X_train), as.character(predomics_res_list$FI_allfeat$summary$feature)))
  {
    # get all species and specify if there are from FBM or not with IsIndSp
    subdf <- predomics_res_list$FI_allfeat$summary[, c('feature', 'value','sign')]
    
    subdf$value[is.nan(subdf$value)] <- 0
    subdf$value = round(subdf$value, digits = 3)
    subdf$sign = ifelse(subdf$sign == -1, "INSHORE", "OFFSHORE")
    fbm.species = predomics_res_list$FI_fmbFeats$summary$feature
    subdf$IsIndsp = as.integer(subdf$feature %in% fbm.species)
    colnames(subdf) <- c("feature","featureImportance", "class","IsIndSp")
    #compute prevalence of each species
    subdf$prevalence <- (colSums(ilevels.df.train > 0)/nrow(ilevels.df.train))* 100
    predomics_res_list$pred_out_fbm <- subdf
    predout.bin[["INSHORE_OFFSHORE"]] <- predomics_res_list
    ##Analyse permanova; get the species to include (from FBM)
    ##Exclude species absents in all subset samples
    ilevels.df.train <- ilevels.df.train[,colSums(ilevels.df.train>0)>0]
    ilevels.df.test <- ilevels.df.test[,colSums(ilevels.df.test>0)>0]
    #Subset the abundance table to these species
    subdf.features.df <- ilevels.df.train[, fbm.species]
    ## Exclude sample with no species
    subdf.features.df <- subdf.features.df[rowSums(subdf.features.df) != 0, ]
    # get class according to rows of subdf.features.df
    subdf.features.df.class <- ilevels.df[rownames(subdf.features.df),]
    subdf.features.df.class <- ifelse(subdf.features.df.class$hab_inoff=='INSHORE', -1,
                                      ifelse(subdf.features.df.class$hab_inoff=='OFFSHORE', 1, NA))
    which(is.na(subdf.features.df.class))
    #Compute bray curtis distances
    subdf.features.df.jaccard <- vegdist(subdf.features.df, method="jaccard", binary=TRUE)
    #Do the permanova with subset species
    subdf.meta <- data.frame(sample=rownames(subdf.features.df), class=subdf.features.df.class)
    set.seed(100)
    subdf.features.df.jaccard.adonis <- adonis2(subdf.features.df.jaccard~class, data=subdf.meta)
    subdf.features.df.jaccard.adonis <- data.frame(subdf.features.df.jaccard.adonis)
    subdf.features.df.jaccard.adonis <- subdf.features.df.jaccard.adonis["class",,drop=FALSE]
    subdf.features.df.jaccard.adonis$comparison <- "INSHORE_OFFSHORE"
    subdf.features.df.jaccard.adonis$data <- "bin"
    subdf.features.df.jaccard.adonis$source <- paste0(algorithm_language,"Species")
    subdf.features.df.jaccard.adonis$features <- length(subdf$feature)
    
    #Compute bray curtis distances; all community
    ilevels.df.jaccard <- vegdist(ilevels.df.train, method="jaccard", binary=TRUE)
    alldf.meta <- data.frame(sample=rownames(ilevels.df.train), class=ilevels.df.train.class)
    set.seed(100)
    ilevels.df.jaccard.adonis <- adonis2(ilevels.df.jaccard~class, data=alldf.meta)
    ilevels.df.jaccard.adonis <- data.frame(ilevels.df.jaccard.adonis)
    ilevels.df.jaccard.adonis <- ilevels.df.jaccard.adonis["class",,drop=FALSE]
    ilevels.df.jaccard.adonis$comparison <- "INSHORE_OFFSHORE"
    ilevels.df.jaccard.adonis$data <- "bin"
    ilevels.df.jaccard.adonis$source <- paste("allSpecies", algorithm_language, sep = '_')
    ilevels.df.jaccard.adonis$features <- ncol(ilevels.df.train)
    #put all together + save
    ilevels.df.jaccard.adonis.all <- rbind(ilevels.df.jaccard.adonis, subdf.features.df.jaccard.adonis)
    adonis_pred.bin[["INSHORE_OFFSHORE"]] <- ilevels.df.jaccard.adonis.all
  }
  
  ##Get the pred_out table
  predout.bin.sub <- lapply(predout.bin, function(x){x[["pred_out_fbm"]]})
  for(i in names(predout.bin.sub))
  {
    predout.bin.sub[[i]][,"comparison"] <- i
  }
  predout.bin.sub <- do.call("rbind", predout.bin.sub)
  predout.bin.sub$source <- algorithm_language
  predout.bin.sub$data <- "bin"
  
  
  #### Save the datasets
  save_path='/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/'
  save(adonis_pred.bin, #adonis results all comparisons; binary data
       adonis_pred.maxn, # adonis results all comparisons; abudnance data
       predout.bin, # indval results all comparisons, abundance data
       predout.maxn,
       predout.bin.sub,
       predout.maxn.sub, file = paste0(save_path, paste(algorithm_language, "Predomics_all_analyses_final.Rda", sep = "_")))
  
}
