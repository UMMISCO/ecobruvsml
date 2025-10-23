# To run this script, use this command: Rscript codePredomics_all_analyses_prev_treshold.R terinter 10

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
percentage_prev <- args[2]

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
  
  # split data into train and test
  set.seed(123)  # For reproducibility
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
  
  # get species with a fix prevalence percentage
  ilevels.df.train <- filterFeaturesByPrevalence(ilevels.df.train,ilevels.df.train.class,perc.prevalence = percentage_prev, by.class = TRUE)
  ilevels.df.train.class <- ifelse(ilevels.df.train$hab==ilevels[[1]], -1,
                                   ifelse(ilevels.df.train$hab==ilevels[[2]], 1, NA))
  
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
    # get only species from FBM
    subdf <- predomics_res_list$FI_fmbFeats$summary[, c('feature', 'value')]
    colnames(subdf) <- c("feature","featureImportance")
    predomics_res_list$pred_out_fbm <- subdf
    predout.maxn[[paste(ilevels, collapse = "_")]] <- predomics_res_list
    ##Analyse permanova; get the species to include (from FBM)
    ##Exclude species absents in all subset samples
    ilevels.df.train <- ilevels.df.train[,colSums(ilevels.df.train>0)>0]
    ilevels.df.test <- ilevels.df.test[,colSums(ilevels.df.test>0)>0]
    #Subset the abundance table to these species
    subdf.features.df <- ilevels.df.train[, subdf$feature]
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
set.seed(123)  # For reproducibility
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

# get species with a fix prevalence percentage
ilevels.df.train <- filterFeaturesByPrevalence(ilevels.df.train,ilevels.df.train.class,perc.prevalence = percentage_prev, by.class = TRUE)
ilevels.df.train.class <- ifelse(ilevels.df.train$hab==ilevels[[1]], -1,
                                 ifelse(ilevels.df.train$hab==ilevels[[2]], 1, NA))

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
  # get only species from FBM
  subdf <- predomics_res_list$FI_fmbFeats$summary[, c('feature', 'value')]
  colnames(subdf) <- c("feature","featureImportance")
  predomics_res_list$pred_out_fbm <- subdf
  predout.maxn[["INSHORE_OFFSHORE"]] <- predomics_res_list
  ##Analyse permanova; get the species to include (from FBM)
  ##Exclude species absents in all subset samples
  ilevels.df.train <- ilevels.df.train[,colSums(ilevels.df.train>0)>0]
  ilevels.df.test <- ilevels.df.test[,colSums(ilevels.df.test>0)>0]
  #Subset the abundance table to these species
  subdf.features.df <- ilevels.df.train[, subdf$feature]
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
  set.seed(123)  # For reproducibility
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
  
  # get species with a fix prevalence percentage
  ilevels.df.train <- filterFeaturesByPrevalence(ilevels.df.train,ilevels.df.train.class,perc.prevalence = percentage_prev, by.class = TRUE)
  ilevels.df.train.class <- ifelse(ilevels.df.train$hab==ilevels[[1]], -1,
                                   ifelse(ilevels.df.train$hab==ilevels[[2]], 1, NA))
  
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
    # get only species from FBM
    subdf <- predomics_res_list$FI_fmbFeats$summary[, c('feature', 'value')]
    colnames(subdf) <- c("feature","featureImportance")
    predomics_res_list$pred_out_fbm <- subdf
    predout.bin[[paste(ilevels, collapse = "_")]] <- predomics_res_list
    ##Analyse permanova; get the species to include (from FBM)
    ##Exclude species absents in all subset samples
    ilevels.df.train <- ilevels.df.train[,colSums(ilevels.df.train>0)>0]
    ilevels.df.test <- ilevels.df.test[,colSums(ilevels.df.test>0)>0]
    #Subset the abundance table to these species
    subdf.features.df <- ilevels.df.train[, subdf$feature]
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
    subdf.features.df.jaccard.adonis$data <- "pres/abs"
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
    ilevels.df.jaccard.adonis$data <- "pres/abs"
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
set.seed(123)  # For reproducibility
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

# get species with a fix prevalence percentage
ilevels.df.train <- filterFeaturesByPrevalence(ilevels.df.train,ilevels.df.train.class,perc.prevalence = percentage_prev, by.class = TRUE)
ilevels.df.train.class <- ifelse(ilevels.df.train$hab==ilevels[[1]], -1,
                                 ifelse(ilevels.df.train$hab==ilevels[[2]], 1, NA))

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
  # get only species from FBM
  subdf <- predomics_res_list$FI_fmbFeats$summary[, c('feature', 'value')]
  colnames(subdf) <- c("feature","featureImportance")
  predomics_res_list$pred_out_fbm <- subdf
  predout.bin[["INSHORE_OFFSHORE"]] <- predomics_res_list
  ##Analyse permanova; get the species to include (from FBM)
  ##Exclude species absents in all subset samples
  ilevels.df.train <- ilevels.df.train[,colSums(ilevels.df.train>0)>0]
  ilevels.df.test <- ilevels.df.test[,colSums(ilevels.df.test>0)>0]
  #Subset the abundance table to these species
  subdf.features.df <- ilevels.df.train[, subdf$feature]
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
  subdf.features.df.jaccard.adonis$data <- "pres/abs"
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
  ilevels.df.jaccard.adonis$data <- "pres/abs"
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
predout.bin.sub$data <- "pres/abs"


#### Save the datasets
save_path='/data/projects/aime/analyses/seamount/metabarcoding/2.analyses_predomics/analysis_results'

# Create the directory if it doesn't exist
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
}

save(adonis_pred.bin, #adonis results all comparisons; binary data
     adonis_pred.maxn, # adonis results all comparisons; abudnance data
     predout.bin, # indval results all comparisons, abundance data
     predout.maxn,
     predout.bin.sub,
     predout.maxn.sub, file = paste0(save_path, paste(algorithm_language,paste0(percentage_prev,"%"),"prev", "Predomics_all_analyses.Rda", sep = "_")))
