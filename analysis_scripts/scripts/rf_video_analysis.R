# RF analysis
# This script is used to identify the set of key species that characterize habitats with efficiency 

if(!require("randomForest")) install.packages("randomForest") ; library(randomForest)
if(!require("caret")) install.packages("caret") ; library(caret)
if (!require("readr")) install.packages("readr"); library(readr)
if (!require("ROCR")) install.packages("ROCR"); library(ROCR)

# load train data from pred/bininter in pres/abs
load(file='/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/Predomics_analysis/terBeam_bininter/Offshore_vs_Inshore/predomicsResults_Offshore_Inshore_pres_abs.Rda')
data_pred_train = list()
data_pred_test = list()
data_pred_train$X = predomics_res_list$comparison_data$X_train
data_pred_train$y = predomics_res_list$comparison_data$y_train
data_pred_test$X = predomics_res_list$comparison_data$X_test
data_pred_test$y = predomics_res_list$comparison_data$y_test

# filter species of train data by indicator species of Indval or ter/bin
indval_species= as.data.frame(read_csv("/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/indicsSpeciesIndval/InValSpecies_video_data_pres_abs.csv", show_col_types = FALSE))
indval_species= indval_species[indval_species$IsIndSpec == 1,]

bin_species = as.data.frame(read_csv("/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/indicSpeciesPredomics/terBeam_bininter/IndicSpecies_predomicsResults_Offshore_Inshore_pres_abs.csv", show_col_types = FALSE))
bin_species= bin_species[bin_species$IsInFBM == 1, ]


ter_species = as.data.frame(read_csv("/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/indicSpeciesPredomics/terBeam_terinter/IndicSpecies_predomicsResults_Offshore_Inshore_pres_abs.csv", show_col_types = FALSE))
ter_species= ter_species[ter_species$IsInFBM == 1, ]

# Rf analysis on indval
set.seed(222)
#Get the species in the indval results on train dataset
indval_species.names <- as.character(indval_species$Species)
#Subset the data_pred_train to indval_species.spnames
data_pred_train_indval <- data_pred_train
data_pred_train_indval$X <- data_pred_train_indval$X[indval_species.names,]
#Add class variable (inshore/offshore) to the X table
data_pred_train_indval_rfinput <- as.data.frame(t(data_pred_train_indval$X))
data_pred_train_indval_rfinput$class <- as.factor(data_pred_train_indval$y)
##need to remove spaces in colnames (species names) before running the rf model (if not error)
colnames(data_pred_train_indval_rfinput) <- gsub(" ",".", colnames(data_pred_train_indval_rfinput))
indval.rf <- randomForest(class~., data = data_pred_train_indval_rfinput, proximity=TRUE) 
# show results of RF on indval species
print(indval.rf)

##Get feature importance
indval.rf_importance <- data.frame(indval.rf$importance)
indval.rf_importance$Species <- gsub("\\."," ", rownames(indval.rf_importance))
##Compare with indicator_value in indval_species results
indval_species_sub <- indval_species[,c("indicator_value","Species")]
indval_species_sub <- merge(indval_species_sub, indval.rf_importance, by="Species")
indval_species_sub$indicator_value <- as.numeric(as.character(indval_species_sub$indicator_value))
plot(indval_species_sub$indicator_value, indval_species_sub$MeanDecreaseGini)

#Add class variable (inshore/offshore) to the X table of test data
data_pred_test_indval_rfinput <- as.data.frame(t(data_pred_test$X))
data_pred_test_indval_rfinput$class <- as.factor(data_pred_test$y)
##need to remove spaces in colnames (species names) before running the rf model (if not error)
colnames(data_pred_test_indval_rfinput) <- gsub(" ",".", colnames(data_pred_test_indval_rfinput))

# make prediction on train data
pred_train <- predict(indval.rf, type= "class")
pred_train

# get confusion matrix for prediction results
indval.pred_train.cm=confusionMatrix(table(pred_train,data_pred_train_indval_rfinput$class))
indval.pred_train.cm

# get acc, auc and f1 metrics of rf on indval species using train data
indval.pred_train.acc <- indval.pred_train.cm$overall['Accuracy']
indval.pred_train.f1 <- indval.pred_train.cm$byClass['F1']

# get the auc of rf model on train data
indval_rf.train.pred=predict(indval.rf, type = "prob")
indval_rf.train.prediction = prediction(indval_rf.train.pred[,2], data_pred_train_indval_rfinput$class)
indval.pred_train.auc= performance(indval_rf.train.prediction, measure = "auc")@y.values[[1]] 

# get metric performances of rf with indval species in training data
rf.metrics = list()
rf.metrics[["indval.rf.train"]] = data.frame(metric= c("auc","accuracy","f1"),
                                       value= c(indval.pred_train.auc, indval.pred_train.acc, indval.pred_train.f1),
                                       dataset= "train",
                                       model="IndVal",
                                       Source= "Presence/Absence")

# make prediction on test data
pred_test <- predict(indval.rf, newdata = data_pred_test_indval_rfinput, type= "class")
pred_test

indval.pred_test.cm=confusionMatrix(table(pred_test,data_pred_test_indval_rfinput$class)) # The prediction to compute the confusion matrix and see the accuracy score 
indval.pred_test.cm

# get acc, auc and f1 metrics of rf on indval species using test data
indval.pred_test.acc <- indval.pred_test.cm$overall['Accuracy']
indval.pred_test.f1 <- indval.pred_test.cm$byClass['F1']

# get the auc of rf model on test data
indval_rf.test.pred=predict(indval.rf, newdata = data_pred_test_indval_rfinput, type = "prob")
indval_rf.test.prediction = prediction(indval_rf.test.pred[,2], data_pred_test_indval_rfinput$class)
indval.pred_test.auc= performance(indval_rf.test.prediction, measure = "auc")@y.values[[1]] 

#get metric performances of rf with indval species in testing data
rf.metrics[["indval.rf.test"]] = data.frame(metric= c("auc","accuracy","f1"),
                                       value= c(indval.pred_test.auc, indval.pred_test.acc, indval.pred_test.f1),
                                       dataset= "test",
                                       model="IndVal",
                                       Source= "Presence/Absence")

# Rf analysis on bininter
set.seed(222)
#Get the species in the bininter results on train dataset
bin_species.names <- as.character(bin_species$Species)
#Subset the data_pred_train to bin_species.spnames
data_pred_train_bin <- data_pred_train
data_pred_train_bin$X <- data_pred_train_bin$X[bin_species.names,]
#Add class variable (inshore/offshore) to the X table
data_pred_train_bin_rfinput <- as.data.frame(t(data_pred_train_bin$X))
data_pred_train_bin_rfinput$class <- as.factor(data_pred_train_bin$y)
##need to remove spaces in colnames (species names) before running the rf model (if not error)
colnames(data_pred_train_bin_rfinput) <- gsub(" ",".", colnames(data_pred_train_bin_rfinput))
bin.rf <- randomForest(class~., data = data_pred_train_bin_rfinput, proximity=TRUE) 
# show results of RF on bininter species
print(bin.rf)

##Get feature importance
bin.rf_importance <- data.frame(bin.rf$importance)
bin.rf_importance$Species <- gsub("\\."," ", rownames(bin.rf_importance))
##Compare with indicator_value in indval_species results
bin_species_sub <- bin_species[,c("featImport","Species")]
bin_species_sub <- merge(bin_species_sub, bin.rf_importance, by="Species")
bin_species_sub$featImport <- as.numeric(as.character(bin_species_sub$featImport))
plot(bin_species_sub$featImport, bin_species_sub$MeanDecreaseGini)

#Add class variable (inshore/offshore) to the X table of test data
data_pred_test_bin_rfinput <- as.data.frame(t(data_pred_test$X))
data_pred_test_bin_rfinput$class <- as.factor(data_pred_test$y)
##need to remove spaces in colnames (species names) before running the rf model (if not error)
colnames(data_pred_test_bin_rfinput) <- gsub(" ",".", colnames(data_pred_test_bin_rfinput))

# make prediction on train data
bin.pred_train <- predict(bin.rf, type= "class")
bin.pred_train

# get confusion matrix for prediction results
bin.pred_train.cm=confusionMatrix(table(bin.pred_train,data_pred_train_bin_rfinput$class))
bin.pred_train.cm

# get acc, auc and f1 metrics of rf on bin species using train data
bin.pred_train.acc <- bin.pred_train.cm$overall['Accuracy']
bin.pred_train.f1 <- bin.pred_train.cm$byClass['F1']

# get the auc of rf model on train data
bin_rf.train.pred=predict(bin.rf, type = "prob")
bin_rf.train.prediction = prediction(bin_rf.train.pred[,2], data_pred_train_bin_rfinput$class)
bin.pred_train.auc= performance(bin_rf.train.prediction, measure = "auc")@y.values[[1]] 

#get metric performances of rf with bininter species in training data
rf.metrics[["bin.rf.train"]] = data.frame(metric= c("auc","accuracy","f1"),
                                            value= c(bin.pred_train.auc, bin.pred_train.acc, bin.pred_train.f1),
                                            dataset= "train",
                                            model="bininter",
                                            Source= "Presence/Absence")

# make prediction on test data
bin.pred_test <- predict(bin.rf, newdata = data_pred_test_bin_rfinput, type= "class")
bin.pred_test

bin.pred_test.cm=confusionMatrix(table(bin.pred_test,data_pred_test_bin_rfinput$class)) # The prediction to compute the confusion matrix and see the accuracy score 
bin.pred_test.cm

# get acc, auc and f1 metrics of rf on indval species using test data
bin.pred_test.acc <- bin.pred_test.cm$overall['Accuracy']
bin.pred_test.f1 <- bin.pred_test.cm$byClass['F1']

# get the auc of rf model on test data
bin_rf.test.pred=predict(bin.rf, newdata = data_pred_test_bin_rfinput, type = "prob")
bin_rf.test.prediction = prediction(bin_rf.test.pred[,2], data_pred_test_bin_rfinput$class)
bin.pred_test.auc= performance(bin_rf.test.prediction, measure = "auc")@y.values[[1]] 

#get metric performances of rf with bininter species in testing data
rf.metrics[["bin.rf.test"]] = data.frame(metric= c("auc","accuracy","f1"),
                                          value= c(bin.pred_test.auc, bin.pred_test.acc, bin.pred_test.f1),
                                          dataset= "test",
                                          model="bininter",
                                          Source= "Presence/Absence")

# Rf analysis on terinter
set.seed(222)
#Get the species in the terinter results on train dataset
ter_species.names <- as.character(ter_species$Species)
#Subset the data_pred_train to ter_species.spnames
data_pred_train_ter <- data_pred_train
data_pred_train_ter$X <- data_pred_train_ter$X[ter_species.names,]
#Add class variable (inshore/offshore) to the X table
data_pred_train_ter_rfinput <- as.data.frame(t(data_pred_train_ter$X))
data_pred_train_ter_rfinput$class <- as.factor(data_pred_train_ter$y)
##need to remove spaces in colnames (species names) before running the rf model (if not error)
colnames(data_pred_train_ter_rfinput) <- gsub(" ",".", colnames(data_pred_train_ter_rfinput))
ter.rf <- randomForest(class~., data = data_pred_train_ter_rfinput, proximity=TRUE) 
# show results of RF on terinter species
print(ter.rf)

##Get feature importance
ter.rf_importance <- data.frame(ter.rf$importance)
ter.rf_importance$Species <- gsub("\\."," ", rownames(ter.rf_importance))
##Compare with indicator_value in indval_species results
ter_species_sub <- ter_species[,c("featImport","Species")]
ter_species_sub <- merge(ter_species_sub, ter.rf_importance, by="Species")
ter_species_sub$featImport <- as.numeric(as.character(ter_species_sub$featImport))
plot(ter_species_sub$featImport, ter_species_sub$MeanDecreaseGini)

#Add class variable (inshore/offshore) to the X table of test data
data_pred_test_ter_rfinput <- as.data.frame(t(data_pred_test$X))
data_pred_test_ter_rfinput$class <- as.factor(data_pred_test$y)
##need to remove spaces in colnames (species names) before running the rf model (if not error)
colnames(data_pred_test_ter_rfinput) <- gsub(" ",".", colnames(data_pred_test_ter_rfinput))

# make prediction on train data
ter.pred_train <- predict(ter.rf, type= "class")
ter.pred_train

# get confusion matrix for prediction results
ter.pred_train.cm=confusionMatrix(table(ter.pred_train,data_pred_train_ter_rfinput$class))
ter.pred_train.cm

# get acc, auc and f1 metrics of rf on terinter species using train data
ter.pred_train.acc <- ter.pred_train.cm$overall['Accuracy']
ter.pred_train.f1 <- ter.pred_train.cm$byClass['F1']

# get the auc of rf model on train data
ter_rf.train.pred=predict(ter.rf, type = "prob")
ter_rf.train.prediction = prediction(ter_rf.train.pred[,2], data_pred_train_ter_rfinput$class)
ter.pred_train.auc= performance(ter_rf.train.prediction, measure = "auc")@y.values[[1]] 

#get metric performances of rf with terinter species in training data
rf.metrics[["ter.rf.train"]] = data.frame(metric= c("auc","accuracy","f1"),
                                         value= c(ter.pred_train.auc, ter.pred_train.acc, ter.pred_train.f1),
                                         dataset= "train",
                                         model="terinter",
                                         Source= "Presence/Absence")

# make prediction on test data
ter.pred_test <- predict(ter.rf, newdata = data_pred_test_ter_rfinput, type= "class")
ter.pred_test

ter.pred_test.cm=confusionMatrix(table(ter.pred_test,data_pred_test_ter_rfinput$class)) # The prediction to compute the confusion matrix and see the accuracy score 
ter.pred_test.cm

# get acc, auc and f1 metrics of rf on indval species using test data
ter.pred_test.acc <- ter.pred_test.cm$overall['Accuracy']
ter.pred_test.f1 <- ter.pred_test.cm$byClass['F1']

# get the auc of rf model on test data
ter_rf.test.pred=predict(ter.rf, newdata = data_pred_test_ter_rfinput, type = "prob")
ter_rf.test.prediction = prediction(ter_rf.test.pred[,2], data_pred_test_ter_rfinput$class)
ter.pred_test.auc= performance(ter_rf.test.prediction, measure = "auc")@y.values[[1]] 

#get metric performances of rf with terinter species in testing data
rf.metrics[["ter.rf.test"]] = data.frame(metric= c("auc","accuracy","f1"),
                                          value= c(ter.pred_test.auc, ter.pred_test.acc, ter.pred_test.f1),
                                          dataset= "test",
                                          model="terinter",
                                          Source= "Presence/Absence")

############################# Random Forest Analysis on abundance data ######################################

# load train data from pred/bininter in abund
load(file='/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/terBeam_bininter/Offshore_vs_Inshore/predomicsResults_Offshore_Inshore.Rda')
data_abund_pred_train = list()
data_abund_pred_test = list()
data_abund_pred_train$X = predomics_res_list$comparison_data$X_train
data_abund_pred_train$y = predomics_res_list$comparison_data$y_train
data_abund_pred_test$X = predomics_res_list$comparison_data$X_test
data_abund_pred_test$y = predomics_res_list$comparison_data$y_test

# filter species of train data by indicator species of Indval or ter/bin
indval_ab_species= as.data.frame(read_csv("/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/InValSpecies_video_data.csv", show_col_types = FALSE))
indval_ab_species= indval_ab_species[indval_ab_species$IsIndSpec == 1,]

bin_ab_species = as.data.frame(read_csv("/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/indicSpeciesPredomics/terBeam_bininter/IndicSpecies_predomicsResults_Offshore_Inshore.csv", show_col_types = FALSE))
bin_ab_species= bin_ab_species[bin_ab_species$IsInFBM == 1, ]

ter_ab_species = as.data.frame(read_csv("/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/indicSpeciesPredomics/terBeam_terinter/IndicSpecies_predomicsResults_Offshore_Inshore.csv", show_col_types = FALSE))
ter_ab_species= ter_ab_species[ter_ab_species$IsInFBM == 1, ]

# Rf analysis on indval
set.seed(222)
#Get the species in the indval results on train dataset
indval_ab_species.names <- as.character(indval_ab_species$Species)
#Subset the data_pred_train to indval_species.spnames
data_pred_train_indval_ab <- data_abund_pred_train
data_pred_train_indval_ab$X <- data_pred_train_indval_ab$X[indval_ab_species.names,]
#Add class variable (inshore/offshore) to the X table
data_pred_train_indval_ab_rfinput <- as.data.frame(t(data_abund_pred_train$X))
data_pred_train_indval_ab_rfinput$class <- as.factor(data_abund_pred_train$y)
##need to remove spaces in colnames (species names) before running the rf model (if not error)
colnames(data_pred_train_indval_ab_rfinput) <- gsub(" ",".", colnames(data_pred_train_indval_ab_rfinput))
indval_ab.rf <- randomForest(class~., data = data_pred_train_indval_ab_rfinput, proximity=TRUE) 
# show results of RF on indval species
print(indval_ab.rf)

##Get feature importance
indval_ab.rf_importance <- data.frame(indval_ab.rf$importance)
indval_ab.rf_importance$Species <- gsub("\\."," ", rownames(indval_ab.rf_importance))
##Compare with indicator_value in indval_species results
indval_ab_species_sub <- indval_ab_species[,c("indicator_value","Species")]
indval_ab_species_sub <- merge(indval_ab_species_sub, indval_ab.rf_importance, by="Species")
indval_ab_species_sub$indicator_value <- as.numeric(as.character(indval_ab_species_sub$indicator_value))
plot(indval_ab_species_sub$indicator_value, indval_ab_species_sub$MeanDecreaseGini)

#Add class variable (inshore/offshore) to the X table of test data
data_pred_test_indval_ab_rfinput <- as.data.frame(t(data_abund_pred_test$X))
data_pred_test_indval_ab_rfinput$class <- as.factor(data_abund_pred_test$y)
##need to remove spaces in colnames (species names) before running the rf model (if not error)
colnames(data_pred_test_indval_ab_rfinput) <- gsub(" ",".", colnames(data_pred_test_indval_ab_rfinput))

# make prediction on train data
pred_train_ab <- predict(indval_ab.rf, type= "class")
pred_train_ab

# get confusion matrix for prediction results
indval_ab.pred_train.cm=confusionMatrix(table(pred_train_ab,data_pred_train_indval_ab_rfinput$class))
indval_ab.pred_train.cm

# get acc, auc and f1 metrics of rf on indval species using train data
indval_ab.pred_train.acc <- indval_ab.pred_train.cm$overall['Accuracy']
indval_ab.pred_train.f1 <- indval_ab.pred_train.cm$byClass['F1']

# get the auc of rf model on train data
indval_ab_rf.train.pred=predict(indval_ab.rf, type = "prob")
indval_ab_rf.train.prediction = prediction(indval_ab_rf.train.pred[,2], data_pred_train_indval_ab_rfinput$class)
indval_ab.pred_train.auc= performance(indval_ab_rf.train.prediction, measure = "auc")@y.values[[1]] 

#get metric performances of rf with indval species in training data in abundance
rf.metrics[["indval_ab.rf.train"]] = data.frame(metric= c("auc","accuracy","f1"),
                                          value= c(indval_ab.pred_train.auc, indval_ab.pred_train.acc, indval_ab.pred_train.f1),
                                          dataset= "train",
                                          model="IndVal",
                                          Source= "Abundance")

# make prediction on test data
pred_test_ab <- predict(indval_ab.rf, newdata = data_pred_test_indval_ab_rfinput, type= "class")
pred_test_ab

indval_ab.pred_test.cm=confusionMatrix(table(pred_test_ab,data_pred_test_indval_ab_rfinput$class)) # The prediction to compute the confusion matrix and see the accuracy score 
indval_ab.pred_test.cm

# get acc, auc and f1 metrics of rf on indval species using test data
indval_ab.pred_test.acc <- indval.pred_test.cm$overall['Accuracy']
indval_ab.pred_test.f1 <- indval.pred_test.cm$byClass['F1']

# get the auc of rf model on test data
indval_ab_rf.test.pred=predict(indval_ab.rf, newdata = data_pred_test_indval_ab_rfinput, type = "prob")
indval_ab_rf.test.prediction = prediction(indval_ab_rf.test.pred[,2], data_pred_test_indval_ab_rfinput$class)
indval_ab.pred_test.auc= performance(indval_ab_rf.test.prediction, measure = "auc")@y.values[[1]] 

#get metric performances of rf with indval species in testing data in abundance
rf.metrics[["indval_ab.rf.test"]] = data.frame(metric= c("auc","accuracy","f1"),
                                                value= c(indval_ab.pred_test.auc, indval_ab.pred_test.acc, indval_ab.pred_test.f1),
                                                dataset= "test",
                                                model="IndVal",
                                                Source= "Abundance")

# Rf analysis on bininter
set.seed(222)
#Get the species in the bininter results on train dataset
bin_ab_species.names <- as.character(bin_ab_species$Species)
#Subset the data_pred_train to bin_species.spnames
data_pred_train_bin_ab <- data_abund_pred_train
data_pred_train_bin_ab$X <- data_pred_train_bin_ab$X[bin_ab_species.names,]
#Add class variable (inshore/offshore) to the X table
data_pred_train_bin_ab_rfinput <- as.data.frame(t(data_abund_pred_train$X))
data_pred_train_bin_ab_rfinput$class <- as.factor(data_abund_pred_train$y)
##need to remove spaces in colnames (species names) before running the rf model (if not error)
colnames(data_pred_train_bin_ab_rfinput) <- gsub(" ",".", colnames(data_pred_train_bin_ab_rfinput))
bin_ab.rf <- randomForest(class~., data = data_pred_train_bin_ab_rfinput, proximity=TRUE) 
# show results of RF on bininter species
print(bin_ab.rf)

##Get feature importance
bin_ab.rf_importance <- data.frame(bin_ab.rf$importance)
bin_ab.rf_importance$Species <- gsub("\\."," ", rownames(bin_ab.rf_importance))
##Compare with indicator_value in bin_species results
bin_ab_species_sub <- bin_ab_species[,c("featImport","Species")]
bin_ab_species_sub <- merge(bin_ab_species_sub, bin_ab.rf_importance, by="Species")
bin_ab_species_sub$featImport <- as.numeric(as.character(bin_ab_species_sub$featImport))
plot(bin_ab_species_sub$featImport, bin_ab_species_sub$MeanDecreaseGini)

#Add class variable (inshore/offshore) to the X table of test data
data_pred_test_bin_ab_rfinput <- as.data.frame(t(data_abund_pred_test$X))
data_pred_test_bin_ab_rfinput$class <- as.factor(data_abund_pred_test$y)
##need to remove spaces in colnames (species names) before running the rf model (if not error)
colnames(data_pred_test_bin_ab_rfinput) <- gsub(" ",".", colnames(data_pred_test_bin_ab_rfinput))

# make prediction on train data
bin.pred_train_ab <- predict(bin_ab.rf, type= "class")
bin.pred_train_ab

# get confusion matrix for prediction results
bin_ab.pred_train.cm=confusionMatrix(table(bin.pred_train_ab,data_pred_train_bin_ab_rfinput$class))
bin_ab.pred_train.cm

# get acc, auc and f1 metrics of rf on indval species using train data
bin_ab.pred_train.acc <- bin_ab.pred_train.cm$overall['Accuracy']
bin_ab.pred_train.f1 <- bin_ab.pred_train.cm$byClass['F1']

# get the auc of rf model on train data
bin_ab_rf.train.pred=predict(bin_ab.rf, type = "prob")
bin_ab_rf.train.prediction = prediction(bin_ab_rf.train.pred[,2], data_pred_train_bin_ab_rfinput$class)
bin_ab.pred_train.auc= performance(bin_ab_rf.train.prediction, measure = "auc")@y.values[[1]] 

#get metric performances of rf with bininter species in training data in abundance
rf.metrics[["bin_ab.rf.train"]] = data.frame(metric= c("auc","accuracy","f1"),
                                               value= c(bin_ab.pred_train.auc, bin_ab.pred_train.acc, bin_ab.pred_train.f1),
                                               dataset= "train",
                                               model="bininter",
                                               Source= "Abundance")

# make prediction on test data
bin.pred_test_ab <- predict(bin_ab.rf, newdata = data_pred_test_bin_ab_rfinput, type= "class")
bin.pred_test_ab

bin_ab.pred_test.cm=confusionMatrix(table(bin.pred_test_ab,data_pred_test_bin_ab_rfinput$class)) # The prediction to compute the confusion matrix and see the accuracy score 
bin_ab.pred_test.cm

# get acc, auc and f1 metrics of rf on bininter species using test data
bin_ab.pred_test.acc <- bin_ab.pred_test.cm$overall['Accuracy']
bin_ab.pred_test.f1 <- bin_ab.pred_test.cm$byClass['F1']

# get the auc of rf model on test data
bin_ab_rf.test.pred=predict(bin_ab.rf, newdata = data_pred_test_bin_ab_rfinput, type = "prob")
bin_ab_rf.test.prediction = prediction(bin_ab_rf.test.pred[,2], data_pred_test_bin_ab_rfinput$class)
bin_ab.pred_test.auc= performance(bin_ab_rf.test.prediction, measure = "auc")@y.values[[1]] 

#get metric performances of rf with bininter species in testing data in abundance
rf.metrics[["bin_ab.rf.test"]] = data.frame(metric= c("auc","accuracy","f1"),
                                             value= c(bin_ab.pred_test.auc, bin_ab.pred_test.acc, bin_ab.pred_test.f1),
                                             dataset= "test",
                                             model="bininter",
                                             Source= "Abundance")

# Rf analysis on terinter
set.seed(222)
#Get the species in the terinter results on train dataset
ter_ab_species.names <- as.character(ter_ab_species$Species)
#Subset the data_pred_train to ter_species.spnames
data_pred_train_ter_ab <- data_abund_pred_train
data_pred_train_ter_ab$X <- data_pred_train_ter_ab$X[ter_ab_species.names,]
#Add class variable (inshore/offshore) to the X table
data_pred_train_ter_ab_rfinput <- as.data.frame(t(data_abund_pred_train$X))
data_pred_train_ter_ab_rfinput$class <- as.factor(data_abund_pred_train$y)
##need to remove spaces in colnames (species names) before running the rf model (if not error)
colnames(data_pred_train_ter_ab_rfinput) <- gsub(" ",".", colnames(data_pred_train_ter_ab_rfinput))
ter_ab.rf <- randomForest(class~., data = data_pred_train_ter_ab_rfinput, proximity=TRUE) 
# show results of RF on indval species
print(ter_ab.rf)

##Get feature importance
ter_ab.rf_importance <- data.frame(ter_ab.rf$importance)
ter_ab.rf_importance$Species <- gsub("\\."," ", rownames(ter_ab.rf_importance))
##Compare with indicator_value in indval_species results
ter_ab_species_sub <- ter_ab_species[,c("featImport","Species")]
ter_ab_species_sub <- merge(ter_ab_species_sub, ter_ab.rf_importance, by="Species")
ter_ab_species_sub$featImport <- as.numeric(as.character(ter_ab_species_sub$featImport))
plot(ter_ab_species_sub$featImport, ter_ab_species_sub$MeanDecreaseGini)

#Add class variable (inshore/offshore) to the X table of test data
data_pred_test_ter_ab_rfinput <- as.data.frame(t(data_abund_pred_test$X))
data_pred_test_ter_ab_rfinput$class <- as.factor(data_abund_pred_test$y)
##need to remove spaces in colnames (species names) before running the rf model (if not error)
colnames(data_pred_test_ter_ab_rfinput) <- gsub(" ",".", colnames(data_pred_test_ter_ab_rfinput))

# make prediction on train data
ter.pred_train_ab <- predict(ter_ab.rf, type= "class")
ter.pred_train_ab

# get confusion matrix for prediction results
ter_ab.pred_train.cm=confusionMatrix(table(ter.pred_train_ab,data_pred_train_ter_ab_rfinput$class))
ter_ab.pred_train.cm

# get acc, auc and f1 metrics of rf on indval species using train data
ter_ab.pred_train.acc <- ter_ab.pred_train.cm$overall['Accuracy']
ter_ab.pred_train.f1 <- ter_ab.pred_train.cm$byClass['F1']

# get the auc of rf model on train data
ter_ab_rf.train.pred=predict(ter_ab.rf, type = "prob")
ter_ab_rf.train.prediction = prediction(ter_ab_rf.train.pred[,2], data_pred_train_ter_ab_rfinput$class)
ter_ab.pred_train.auc= performance(ter_ab_rf.train.prediction, measure = "auc")@y.values[[1]] 

#get metric performances of rf with terinter species in training data in abundance
rf.metrics[["ter_ab.rf.train"]] = data.frame(metric= c("auc","accuracy","f1"),
                                            value= c(ter_ab.pred_train.auc, ter_ab.pred_train.acc, ter_ab.pred_train.f1),
                                            dataset= "train",
                                            model="terinter",
                                            Source= "Abundance")

# make prediction on test data
ter.pred_test_ab <- predict(ter_ab.rf, newdata = data_pred_test_ter_ab_rfinput, type= "class")
ter.pred_test_ab

ter_ab.pred_test.cm=confusionMatrix(table(ter.pred_test_ab,data_pred_test_ter_ab_rfinput$class)) # The prediction to compute the confusion matrix and see the accuracy score 
ter_ab.pred_test.cm

# get acc, auc and f1 metrics of rf on terinter species using test data
ter_ab.pred_test.acc <- ter_ab.pred_test.cm$overall['Accuracy']
ter_ab.pred_test.f1 <- ter_ab.pred_test.cm$byClass['F1']

# get the auc of rf model on test data
ter_ab_rf.test.pred=predict(ter_ab.rf, newdata = data_pred_test_ter_ab_rfinput, type = "prob")
ter_ab_rf.test.prediction = prediction(ter_ab_rf.test.pred[,2], data_pred_test_ter_ab_rfinput$class)
ter_ab.pred_test.auc= performance(ter_ab_rf.test.prediction, measure = "auc")@y.values[[1]]

#get metric performances of rf with terinter species in testing data in abundance
rf.metrics[["ter_ab.rf.test"]] = data.frame(metric= c("auc","accuracy","f1"),
                                             value= c(ter_ab.pred_test.auc, ter_ab.pred_test.acc, ter_ab.pred_test.f1),
                                             dataset= "test",
                                             model="terinter",
                                             Source= "Abundance")

#For each element in rf.metrics, do the rbind of the dataframes of each binary prediction task
rf.metrics.df <- do.call("rbind",rf.metrics)

# remove rownames of metrics df
rownames(rf.metrics.df) <- NULL

# round metrics' values
rf.metrics.df$value <- round(rf.metrics.df$value, digits = 2)

# plot the metrics' value
rf.metrics.plot1 <- ggplot(data=rf.metrics.df, aes(x=metric, y=value, color=metric)) +
                    geom_bar(stat="identity", fill="lightblue", width = 0.6)+
                    geom_text(aes(label=value), vjust=2, size=3, position = position_dodge(0.9), color="black")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
                    labs(color = "Metrics") + 
                    facet_grid(dataset+Source~model)+
                    theme(aspect.ratio = 1)
rf.metrics.plot2 <- ggplot(data=rf.metrics.df, aes(x=metric, y=value, fill=dataset)) +
                geom_bar(stat="identity", position = position_dodge())+
                geom_text(aes(label=value), vjust=2, size=3, position = position_dodge(0.9), color="black")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
                facet_grid(Source~model)+
                theme(aspect.ratio = 1)

# There is no need to perform cross validation on a rf model to get an unbiased estimate of the test set error:
  #reference: http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm#ooberr

# But there is a cross validation function in rf for feature selection

indval.rf.cv <- rfcv(data_pred_test_indval_rfinput[, -ncol(data_pred_test_indval_rfinput)],  data_pred_test_indval_rfinput$class, cv.fold = 10)

plot(indval.rf.cv)

# Plot cross validation versus model producers accuracy
par(mfrow=c(1,2))
  plot(indval.rf.cv$n.var, indval.rf.cv$error.cv, main = "CV producers accuracy")
  plot(indval.rf.cv$n.var, indval.rf.cv$predicted, main = "Model producers accuracy")

# show results of RF on indval species
print(indval.rf)