# How to run this script
## make sure you are in the script directory
# cd /data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/scripts
# Rscript codeIndval_all_analyses_prev.R 0

if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("reshape2")) install.packages("reshape2"); library(reshape2)
if (!require("vegan")) install.packages("vegan"); library(vegan)
if (!require("stats")) install.packages("stats"); library(stats)
if (!require("labdsv")) install.packages("labdsv"); library(labdsv)

# define arguments for the script
args <- commandArgs(trailingOnly = TRUE)
prevalence_rate <- as.numeric(args[1])

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
           X <- sample_data[,-match(c("Station","Site","hab","hab_inoff"), colnames(sample_data))]
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
           sample_data.X <- sample_data[,-match(c("Station","Site","hab","hab_inoff"), colnames(sample_data))]
           
           # Return sample.info, sample.X, sample.class
           return(list(sample.info = sample_data, sample.X = sample_data.X, sample.class = sample_data.class))
         }
  )
  
}

# Function to save predomics results object
save_indval_results <- function(results_indval, X, ilevels){
  
  # get result table from indval object
  subdf <- data.frame(results_indval$indcls, results_indval$pval, results_indval$maxcls)
  colnames(subdf) <- c("featureImportance","pval", "class")
  subdf$feature <- rownames(subdf)
  # replace NA values by 0
  subdf$pval[is.nan(subdf$pval)] <- 0
  # round feature importance value
  subdf$featureImportance = round(subdf$featureImportance, digits = 3)
  # round pvalue
  subdf$pval = round(subdf$pval, digits = 3)
  # assign the class 1 to the first habitat and 2 to the second for each comparison
  subdf$class <- ifelse(subdf$class == 1, ilevels[[1]], ifelse(subdf$class == 2, ilevels[[2]], NA))
  # If the species pval <0.05, the species is an indicator 
  subdf$IsIndSp <- ifelse(subdf$pval<0.05 , 1, 0)
  #compute prevalence of each species
  subdf$prevalence <- (colSums(X> 0)/nrow(X))* 100
  # round pvalue
  subdf$prevalence = round(subdf$prevalence, digits = 3)
  # reorder colnames
  subdf <- subdf[, c("feature","featureImportance","pval","class","IsIndSp","prevalence")]
  results_indval$indval_pval <- subdf
  
  return(results_indval)
  
}
# compute permanova analysis
permanova_analysis <- function(sample.info, sample.X, sample.class, results_indval, hab_type, ilevels, data_type) {
  
  # Initialize a list to store results
  adonis <- list()
  
  if (data_type == "maxN") {
    # Get indics species names
    indics.species <- names(results_indval$indcls[results_indval$pval < 0.05])
    
    # Subset the abundance table to these species
    subdf.features.df <- sample.X[, indics.species]
    
    # Exclude samples with no species
    subdf.features.df <- subdf.features.df[rowSums(subdf.features.df) > 0, ]
    
    # Get class according to rows of subdf.features.df
    subdf.features.df.class <- sample.info[rownames(subdf.features.df),]
    subdf.features.df.class <- ifelse(subdf.features.df.class[[hab_type]] == ilevels[[1]], 1,
                                      ifelse(subdf.features.df.class[[hab_type]] == ilevels[[2]], 2, NA))
    
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
    subdf.features.df.bray.adonis$source <- "indvalSpecies"
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
    ilevels.df.bray.adonis$source <- "allSpecies"
    ilevels.df.bray.adonis$features <- ncol(sample.X)
    
    # Combine results
    adonis_pred <- rbind(ilevels.df.bray.adonis, subdf.features.df.bray.adonis)
    
  } else if (data_type == "pres/abs") {
    # Similar steps for presence/absence data
    # Get indics species names
    indics.species <- names(results_indval$indcls[results_indval$pval < 0.05])
    
    # Subset the abundance table to these species
    subdf.features.df <- sample.X[, indics.species]
    
    # Exclude samples with no species
    subdf.features.df <- subdf.features.df[rowSums(subdf.features.df) > 0, ]
    
    # Get class according to rows of subdf.features.df
    subdf.features.df.class <- sample.info[rownames(subdf.features.df),]
    subdf.features.df.class <- ifelse(subdf.features.df.class[[hab_type]] == ilevels[[1]], 1,
                                      ifelse(subdf.features.df.class[[hab_type]] == ilevels[[2]], 2, NA))
    
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
    subdf.features.df.jaccard.adonis$source <- "indvalSpecies"
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
    ilevels.df.jaccard.adonis$source <- "allSpecies"
    ilevels.df.jaccard.adonis$features <- ncol(sample.X)
    
    # Combine results
    adonis_pred <- rbind(ilevels.df.jaccard.adonis, subdf.features.df.jaccard.adonis)
  } else {
    stop("Invalid data_type specified. Use 'maxN' or 'pres/abs'.")
  }
  
  return(adonis_pred)
}

####################
## Analyses on abundance data (MaxN)
####################
#library(labdsv)
indvalout.maxn <- list()
adonis.maxn <- list()
for(i in 1:length(comp))
{
  # print comparison
  print(comp[[i]])
  ##Get the levels to compare
  ilevels <- comp[[i]]
  ##Get the data limited to the levels to compare
  ilevels.df <- acc[acc$hab %in% ilevels,]
  # prepare data for indval analysis
  data <- prepare_data_for_analysis(ilevels.df, habitat_type = "hab", ilevels[[1]], ilevels[[2]])
  # get elements returned
  ilevels.df= data$sample.info
  ilevels.df.X= data$sample.X
  ilevels.df.class= data$sample.class
  ##Do the indval analyses
  ilevels.df.indval <- labdsv::indval(x = ilevels.df.X, clustering = ilevels.df.class, numitr = 9999)
  ##Merge the indval and pvalues output in single table + save
  if(identical(names(ilevels.df.indval$indcls), names(ilevels.df.indval$pval)))
  {
    # save indval results
    indvalout.maxn[[paste(ilevels, collapse = "_")]] <- save_indval_results(ilevels.df.indval, ilevels.df.X, ilevels)
    ##Analyse permanova; get the species to include (indval pvalue<0.05)
    adonis.maxn[[paste(ilevels, collapse = "_")]] <- permanova_analysis(ilevels.df, ilevels.df.X,ilevels.df.class, ilevels.df.indval, hab_type = "hab", ilevels, data_type = "maxN")
  }
}

#Build a comparison vector for INSHORE_OFFSHORE  
ilevels.df <- acc

# define a list for habitat INSHORE/OFFSHORE
ilevels_hab <- list("INSHORE", "OFFSHORE")

# print comparison
print(ilevels_hab)

# prepare data for indval analysis
data <- prepare_data_for_analysis(ilevels.df, habitat_type = "hab_inoff", ilevels_hab[[1]], ilevels_hab[[2]])
# get elements returned
ilevels.df= data$sample.info
ilevels.df.X= data$sample.X
ilevels.df.class= data$sample.class

##Do the indval analyses
ilevels.df.indval <- labdsv::indval(x = ilevels.df.X, clustering = ilevels.df.class, numitr = 9999)
##Merge the indval and pvalues output in single table + save
if(identical(names(ilevels.df.indval$indcls), names(ilevels.df.indval$pval)))
{
  # save indval results
  indvalout.maxn[[paste(ilevels_hab, collapse = "_")]] <- save_indval_results(ilevels.df.indval, ilevels.df.X, ilevels_hab)
  ##Analyse permanova; get the species to include (indval pvalue<0.05)
  adonis.maxn[[paste(ilevels_hab, collapse = "_")]] <- permanova_analysis(ilevels.df, ilevels.df.X,ilevels.df.class, ilevels.df.indval, hab_type = "hab_inoff", ilevels_hab, data_type = "maxN")
}

##Get the indval_pval table
indvalout.maxn.sub <- lapply(indvalout.maxn, function(x){x[["indval_pval"]]})
for(i in names(indvalout.maxn.sub))
{
  indvalout.maxn.sub[[i]][,"comparison"] <- i
}
indvalout.maxn.sub <- do.call("rbind", indvalout.maxn.sub)
indvalout.maxn.sub$source <- "indval"
indvalout.maxn.sub$data <- "maxN"


####################
## Analyses on presence/absence data
####################

library(labdsv)
indvalout.bin <- list()
adonis.bin <- list()
for(i in 1:length(comp))
{
  # print comparison
  print(comp[[i]])
  ##Get the levels to compare
  ilevels <- comp[[i]]
  ##Get the data limited to the levels to compare
  ilevels.df <- acc[acc$hab %in% ilevels,]
  # prepare data for indval analysis
  data <- prepare_data_for_analysis(ilevels.df, habitat_type = "hab", ilevels[[1]], ilevels[[2]])
  # get elements returned
  ilevels.df= data$sample.info
  ilevels.df.X= data$sample.X
  ilevels.df.class= data$sample.class
  ##Transform on binary
  ilevels.df.X <- as.data.frame(apply(ilevels.df.X, 2, function(x){ifelse(x==0,0,1)}))
  ##Do the indval analyses
  ilevels.df.indval <- labdsv::indval(x = ilevels.df.X, clustering = ilevels.df.class, numitr = 9999)
  ##Merge the indval and pvalues output in single table + save
  if(identical(names(ilevels.df.indval$indcls), names(ilevels.df.indval$pval)))
  {
    # save indval results
    indvalout.bin[[paste(ilevels, collapse = "_")]] <- save_indval_results(ilevels.df.indval, ilevels.df.X, ilevels)
    ##Analyse permanova; get the species to include (indval pvalue<0.05)
    adonis.bin[[paste(ilevels, collapse = "_")]] <- permanova_analysis(ilevels.df, ilevels.df.X,ilevels.df.class, ilevels.df.indval, hab_type = "hab", ilevels, data_type = "pres/abs")
  }
}

#Build a comparison vector for INSHORE_OFFSHORE  
ilevels.df <- acc

# define a list for habitat INSHORE/OFFSHORE
ilevels_hab <- list("INSHORE", "OFFSHORE")

# print comparison
print(ilevels_hab)

# prepare data for indval analysis
data <- prepare_data_for_analysis(ilevels.df, habitat_type = "hab_inoff", ilevels_hab[[1]], ilevels_hab[[2]])
# get elements returned
ilevels.df= data$sample.info
ilevels.df.X= data$sample.X
ilevels.df.class= data$sample.class

##Transform to binary
ilevels.df.X <- as.data.frame(apply(ilevels.df.X, 2, function(x){ifelse(x==0,0,1)}))
##Do the indval analyses
ilevels.df.indval <- labdsv::indval(x = ilevels.df.X, clustering = ilevels.df.class, numitr = 9999)
##Merge the indval and pvalues output in single table + save
if(identical(names(ilevels.df.indval$indcls), names(ilevels.df.indval$pval)))
{
  
  # save indval results
  indvalout.bin[[paste(ilevels_hab, collapse = "_")]] <- save_indval_results(ilevels.df.indval, ilevels.df.X, ilevels_hab)
  ##Analyse permanova; get the species to include (indval pvalue<0.05)
  adonis.bin[[paste(ilevels_hab, collapse = "_")]] <- permanova_analysis(ilevels.df, ilevels.df.X,ilevels.df.class, ilevels.df.indval, hab_type = "hab_inoff", ilevels_hab, data_type = "pres/abs")
  
}

##Get the indval_pval table
indvalout.bin.sub <- lapply(indvalout.bin, function(x){x[["indval_pval"]]})
for(i in names(indvalout.bin.sub))
{
  indvalout.bin.sub[[i]][,"comparison"] <- i
}
indvalout.bin.sub <- do.call("rbind", indvalout.bin.sub)
indvalout.bin.sub$source <- "indval"
indvalout.bin.sub$data <- "pres/abs"


##### Save the datasets
save_path=paste0('/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/indval_output_data/') 
save(adonis.bin, #adonis results all comparisons; binary data
     adonis.maxn, # adonis results all comparisons; abudnance data
     indvalout.bin, # indval results all comparisons, abundance data
     indvalout.maxn, 
     indvalout.bin.sub, 
     indvalout.maxn.sub, file = paste0(save_path, paste("Indval_all_analyses_overall_data", "prev", paste0(prevalence_rate,".Rda"), sep = "_")))