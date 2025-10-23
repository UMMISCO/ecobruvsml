if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("reshape2")) install.packages("reshape2"); library(reshape2)
if (!require("vegan")) install.packages("vegan"); library(vegan)
if (!require("stats")) install.packages("stats"); library(stats)
if (!require("labdsv")) install.packages("labdsv"); library(labdsv)
if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if (!require("ggrepel")) install.packages("ggrepel"); library(ggrepel)

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
#library(labdsv)
indvalout.maxn <- list()
adonis.maxn <- list()
for(i in 1:length(comp))
{
  print(i)
  ##Get the levels to compare
  ilevels <- comp[[i]]
  ##Get the data limited to the levels to compare
  ilevels.df <- acc[acc$hab %in% ilevels,]
  ##Build the class vector to do the indval comparison
  ilevels.df.class <- ilevels.df$hab
  ##Exclude unnecesary columns (keep abundance table only samples x species)
  ilevels.df <- ilevels.df[,-match(c("Station","Site","hab","hab_inoff"), colnames(acc))]
  ##Exclude species absents in all subset samples
  ilevels.df <- ilevels.df[,colSums(ilevels.df>0)>0]
  ##Do the indval analyses
  ilevels.df.indval <- labdsv::indval(x = ilevels.df, clustering = ilevels.df.class, numitr = 9999)
  ##Merge the indval and pvalues output in single table + save
  if(identical(names(ilevels.df.indval$indcls), names(ilevels.df.indval$pval)))
  {
    subdf <- data.frame(ilevels.df.indval$indcls, ilevels.df.indval$pval, ilevels.df.indval$maxcls)
    colnames(subdf) <- c("featureImportance","pval", "class")
    subdf$feature <- rownames(subdf)
    subdf$featureImportance = round(subdf$featureImportance, digits = 3)
    subdf$pval = round(subdf$pval, digits = 3)
    # assign the class 1 to the first habitat and 2 to the second for each comparison
    subdf$class <- ifelse(subdf$class == 1, ilevels[[1]], ifelse(subdf$class == 2, ilevels[[2]], NA))
    # If the species pval <0.05, the species is an indicator 
    subdf$IsIndSp <- ifelse(subdf$pval<0.05 , 1, 0)
    #compute prevalence of each species
    subdf$prevalence <- (colSums(ilevels.df > 0)/nrow(ilevels.df))* 100
    #ilevels.df.indval <- subdf
    # reorder colnames
    subdf <- subdf[, c("feature","featureImportance","pval","class","IsIndSp","prevalence")]
    ilevels.df.indval$indval_pval <- subdf
    indvalout.maxn[[paste(ilevels, collapse = "_")]] <- ilevels.df.indval
    ##Analyse permanova; get the species to include (indval pvalue<0.05)
    subdf.features <- subdf$feature[subdf$pval<0.05]
    #Subset the abundance table to these species
    subdf.features.df <- ilevels.df[,subdf.features]
    #Compute bray curtis distances 
    subdf.features.df.bray <- vegdist(subdf.features.df, method="bray")
    #Do the permanova with subset species
    subdf.meta <- data.frame(sample=rownames(ilevels.df), class=ilevels.df.class)
    set.seed(100)
    subdf.features.df.bray.adonis <- adonis2(subdf.features.df.bray~class, data=subdf.meta)
    subdf.features.df.bray.adonis <- data.frame(subdf.features.df.bray.adonis)
    subdf.features.df.bray.adonis <- subdf.features.df.bray.adonis["class",,drop=FALSE]
    subdf.features.df.bray.adonis$comparison <- paste(ilevels, collapse = "_")
    subdf.features.df.bray.adonis$data <- "maxN"
    subdf.features.df.bray.adonis$source <- "indvalSpecies"
    subdf.features.df.bray.adonis$features <- length(subdf.features)
    
    #Compute bray curtis distances; all community
    ilevels.df.bray <- vegdist(ilevels.df, method="bray")
    set.seed(100)
    ilevels.df.bray.adonis <- adonis2(ilevels.df.bray~class, data=subdf.meta)
    ilevels.df.bray.adonis <- data.frame(ilevels.df.bray.adonis)
    ilevels.df.bray.adonis <- ilevels.df.bray.adonis["class",,drop=FALSE]
    ilevels.df.bray.adonis$comparison <- paste(ilevels, collapse = "_")
    ilevels.df.bray.adonis$data <- "maxN"
    ilevels.df.bray.adonis$source <- "allSpecies"
    ilevels.df.bray.adonis$features <- ncol(ilevels.df)
    
    #put all together + save
    ilevels.df.bray.adonis.all <- rbind(ilevels.df.bray.adonis, subdf.features.df.bray.adonis)
    adonis.maxn[[paste(ilevels, collapse = "_")]] <- ilevels.df.bray.adonis.all
  }
}

##Build the class vector to do the indval comparison
ilevels.df.class <- acc$hab_inoff
##Exclude unnecesary columns (keep abundance table only samples x species)
ilevels.df <- acc[,-match(c("Station","Site","hab","hab_inoff"), colnames(acc))]
##Exclude species absents in all subset samples
ilevels.df <- ilevels.df[,colSums(ilevels.df>0)>0]
##Do the indval analyses
ilevels.df.indval <- labdsv::indval(x = ilevels.df, clustering = ilevels.df.class, numitr = 9999)
##Merge the indval and pvalues output in single table + save
if(identical(names(ilevels.df.indval$indcls), names(ilevels.df.indval$pval)))
{
  print("INSHORE_OFFSHORE")
  subdf <- data.frame(ilevels.df.indval$indcls, ilevels.df.indval$pval, ilevels.df.indval$maxcls)
  colnames(subdf) <- c("featureImportance","pval", "class")
  subdf$feature <- rownames(subdf)
  subdf$featureImportance = round(subdf$featureImportance, digits = 3)
  subdf$pval = round(subdf$pval, digits = 3)
  # assign the class 1 to the first habitat and 2 to the second for each comparison
  subdf$class <- ifelse(subdf$class == 1, "INSHORE", ifelse(subdf$class == 2, "OFFSHORE", NA))
  # If the species pval <0.05, the species is an indicator 
  subdf$IsIndSp <- ifelse(subdf$pval<0.05 , 1, 0)
  #compute prevalence of each species
  subdf$prevalence <- (colSums(ilevels.df > 0)/nrow(ilevels.df))* 100
  #ilevels.df.indval <- subdf
  # reorder colnames
  subdf <- subdf[, c("feature","featureImportance","pval","class","IsIndSp","prevalence")]
  ilevels.df.indval$indval_pval <- subdf
  indvalout.maxn[["INSHORE_OFFSHORE"]] <- ilevels.df.indval
  ##Analyse permanova; get the species to include (indval pvalue<0.05)
  subdf.features <- subdf$feature[subdf$pval<0.05]
  #Subset the abundance table to these species
  subdf.features.df <- ilevels.df[,subdf.features]
  #Compute bray curtis distances
  subdf.features.df.bray <- vegdist(subdf.features.df, method="bray")
  #Do the permanova with subset species
  subdf.meta <- data.frame(sample=rownames(ilevels.df), class=ilevels.df.class)
  set.seed(100)
  subdf.features.df.bray.adonis <- adonis2(subdf.features.df.bray~class, data=subdf.meta)
  subdf.features.df.bray.adonis <- data.frame(subdf.features.df.bray.adonis)
  subdf.features.df.bray.adonis <- subdf.features.df.bray.adonis["class",,drop=FALSE]
  subdf.features.df.bray.adonis$comparison <- "INSHORE_OFFSHORE"
  subdf.features.df.bray.adonis$data <- "maxN"
  subdf.features.df.bray.adonis$source <- "indvalSpecies"
  subdf.features.df.bray.adonis$features <- length(subdf.features)
  
  #Compute bray curtis distances; all community
  ilevels.df.bray <- vegdist(ilevels.df, method="bray")
  set.seed(100)
  ilevels.df.bray.adonis <- adonis2(ilevels.df.bray~class, data=subdf.meta)
  ilevels.df.bray.adonis <- data.frame(ilevels.df.bray.adonis)
  ilevels.df.bray.adonis <- ilevels.df.bray.adonis["class",,drop=FALSE]
  ilevels.df.bray.adonis$comparison <- "INSHORE_OFFSHORE"
  ilevels.df.bray.adonis$data <- "maxN"
  ilevels.df.bray.adonis$source <- "allSpecies"
  ilevels.df.bray.adonis$features <- ncol(ilevels.df)
  
  #put all together + save
  ilevels.df.bray.adonis.all <- rbind(ilevels.df.bray.adonis, subdf.features.df.bray.adonis)
  adonis.maxn[['INSHORE_OFFSHORE']] <- ilevels.df.bray.adonis.all
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
  print(i)
  ##Get the levels to compare
  ilevels <- comp[[i]]
  ##Get the data limited to the levels to compare
  ilevels.df <- acc[acc$hab %in% ilevels,]
  ##Build the class vector to do the indval comparison
  ilevels.df.class <- ilevels.df$hab
  ##Exclude unnecesary columns (keep abundance table only samples x species)
  ilevels.df <- ilevels.df[,-match(c("Station","Site","hab","hab_inoff"), colnames(acc))]
  ##Exclude species absents in all subset samples
  ilevels.df <- ilevels.df[,colSums(ilevels.df>0)>0]
  ##Transform on binary
  ilevels.df <- as.data.frame(apply(ilevels.df, 2, function(x){ifelse(x==0,0,1)}))
  ##Do the indval analyses
  ilevels.df.indval <- labdsv::indval(x = ilevels.df, clustering = ilevels.df.class, numitr = 9999)
  ##Merge the indval and pvalues output in single table + save
  if(identical(names(ilevels.df.indval$indcls), names(ilevels.df.indval$pval)))
  {
    subdf <- data.frame(ilevels.df.indval$indcls, ilevels.df.indval$pval, ilevels.df.indval$maxcls)
    colnames(subdf) <- c("featureImportance","pval", "class")
    subdf$feature <- rownames(subdf)
    subdf$featureImportance = round(subdf$featureImportance, digits = 3)
    subdf$pval = round(subdf$pval, digits = 3)
    # assign the class 1 to the first habitat and 2 to the second for each comparison
    subdf$class <- ifelse(subdf$class == 1, ilevels[[1]], ifelse(subdf$class == 2, ilevels[[2]], NA))
    # If the species pval <0.05, the species is an indicator 
    subdf$IsIndSp <- ifelse(subdf$pval<0.05 , 1, 0)
    #compute prevalence of each species
    subdf$prevalence <- (colSums(ilevels.df > 0)/nrow(ilevels.df))* 100
    #ilevels.df.indval <- subdf
    # reorder colnames
    subdf <- subdf[, c("feature","featureImportance","pval","class","IsIndSp","prevalence")]
    ilevels.df.indval$indval_pval <- subdf
    indvalout.bin[[paste(ilevels, collapse = "_")]] <- ilevels.df.indval
    ##Analyse permanova; get the species to include (indval pvalue<0.05)
    subdf.features <- subdf$feature[subdf$pval<0.05]
    #Subset the abundance table to these species
    subdf.features.df <- ilevels.df[,subdf.features]
    ## Exclude sample with no species
    subdf.features.df <- subdf.features.df[rowSums(subdf.features.df) != 0, ]
    # get class according to rows of subdf.features.df
    subdf.features.df.class <- acc[acc$hab %in% ilevels,]
    subdf.features.df.class <- subdf.features.df.class[rownames(subdf.features.df),]
    subdf.features.df.class <- ifelse(subdf.features.df.class$hab==ilevels[[1]], -1,
                                      ifelse(subdf.features.df.class$hab==ilevels[[2]], 1, NA))
    #Compute jaccard distances
    subdf.features.df.jaccard <- vegdist(subdf.features.df, method="jaccard", binary=TRUE)
    #Do the permanova
    subdf.meta <- data.frame(sample=rownames(subdf.features.df), class=subdf.features.df.class)
    set.seed(100)
    subdf.features.df.jaccard.adonis <- adonis2(subdf.features.df.jaccard~class, data=subdf.meta)
    subdf.features.df.jaccard.adonis <- data.frame(subdf.features.df.jaccard.adonis)
    subdf.features.df.jaccard.adonis <- subdf.features.df.jaccard.adonis["class",,drop=FALSE]
    subdf.features.df.jaccard.adonis$comparison <- paste(ilevels, collapse = "_")
    subdf.features.df.jaccard.adonis$data <- "bin"
    subdf.features.df.jaccard.adonis$source <- "indvalSpecies"
    subdf.features.df.jaccard.adonis$features <- length(subdf.features)
    
    #Compute jaccard distances; all community
    ilevels.df.jaccard <- vegdist(ilevels.df, method="jaccard", binary = TRUE)
    subdf.meta <- data.frame(sample=rownames(ilevels.df), class=ilevels.df.class)
    set.seed(100)
    ilevels.df.jaccard.adonis <- adonis2(ilevels.df.jaccard~class, data=subdf.meta)
    ilevels.df.jaccard.adonis <- data.frame(ilevels.df.jaccard.adonis)
    ilevels.df.jaccard.adonis <- ilevels.df.jaccard.adonis["class",,drop=FALSE]
    ilevels.df.jaccard.adonis$comparison <- paste(ilevels, collapse = "_")
    ilevels.df.jaccard.adonis$data <- "bin"
    ilevels.df.jaccard.adonis$source <- "allSpecies"
    ilevels.df.jaccard.adonis$features <- ncol(ilevels.df)
    
    #put all together + save
    ilevels.df.jaccard.adonis.all <- rbind(ilevels.df.jaccard.adonis, subdf.features.df.jaccard.adonis)
    adonis.bin[[paste(ilevels, collapse = "_")]] <- ilevels.df.jaccard.adonis.all
  }
}

##Build the class vector to do the indval comparison
ilevels.df.class <- acc$hab_inoff
##Exclude unnecesary columns (keep abundance table only samples x species)
ilevels.df <- acc[,-match(c("Station","Site","hab","hab_inoff"), colnames(acc))]
##Exclude species absents in all subset samples
ilevels.df <- ilevels.df[,colSums(ilevels.df>0)>0]
##Transform to binary
ilevels.df <- as.data.frame(apply(ilevels.df, 2, function(x){ifelse(x==0,0,1)}))
##Do the indval analyses
ilevels.df.indval <- labdsv::indval(x = ilevels.df, clustering = ilevels.df.class, numitr = 9999)
##Merge the indval and pvalues output in single table + save
if(identical(names(ilevels.df.indval$indcls), names(ilevels.df.indval$pval)))
{
  print("INSHORE_OFFSHORE")
  subdf <- data.frame(ilevels.df.indval$indcls, ilevels.df.indval$pval, ilevels.df.indval$maxcls)
  colnames(subdf) <- c("featureImportance","pval", "class")
  subdf$feature <- rownames(subdf)
  subdf$featureImportance = round(subdf$featureImportance, digits = 3)
  subdf$pval = round(subdf$pval, digits = 3)
  # assign the class 1 to the first habitat and 2 to the second for each comparison
  subdf$class <- ifelse(subdf$class == 1, "INSHORE", ifelse(subdf$class == 2, "OFFSHORE", NA))
  # If the species pval <0.05, the species is an indicator 
  subdf$IsIndSp <- ifelse(subdf$pval<0.05 , 1, 0)
  #compute prevalence of each species
  subdf$prevalence <- (colSums(ilevels.df > 0)/nrow(ilevels.df))* 100
  # reorder colnames
  subdf <- subdf[, c("feature","featureImportance","pval","class","IsIndSp","prevalence")]
  ilevels.df.indval$indval_pval <- subdf
  indvalout.bin[["INSHORE_OFFSHORE"]] <- ilevels.df.indval
  ##Analyse permanova; get the species to include (indval pvalue<0.05)
  subdf.features <- subdf$feature[subdf$pval<0.05]
  #Subset the abundance table to these species
  subdf.features.df <- ilevels.df[,subdf.features]
  #Compute jaccard distances
  subdf.features.df.jaccard <- vegdist(subdf.features.df, method="jaccard", binary=TRUE)
  #Do the permanova
  subdf.meta <- data.frame(sample=rownames(ilevels.df), class=ilevels.df.class)
  set.seed(100)
  subdf.features.df.jaccard.adonis <- adonis2(subdf.features.df.jaccard~class, data=subdf.meta)
  subdf.features.df.jaccard.adonis <- data.frame(subdf.features.df.jaccard.adonis)
  subdf.features.df.jaccard.adonis <- subdf.features.df.jaccard.adonis["class",,drop=FALSE]
  subdf.features.df.jaccard.adonis$comparison <- "INSHORE_OFFSHORE"
  subdf.features.df.jaccard.adonis$data <- "bin"
  subdf.features.df.jaccard.adonis$source <- "indvalSpecies"
  subdf.features.df.jaccard.adonis$features <- length(subdf.features)
  
  #Compute jaccard distances; all community
  ilevels.df.jaccard <- vegdist(ilevels.df, method="jaccard", binary = TRUE)
  subdf.meta <- data.frame(sample=rownames(ilevels.df), class=ilevels.df.class)
  set.seed(100)
  ilevels.df.jaccard.adonis <- adonis2(ilevels.df.jaccard~class, data=subdf.meta)
  ilevels.df.jaccard.adonis <- data.frame(ilevels.df.jaccard.adonis)
  ilevels.df.jaccard.adonis <- ilevels.df.jaccard.adonis["class",,drop=FALSE]
  ilevels.df.jaccard.adonis$comparison <- "INSHORE_OFFSHORE"
  ilevels.df.jaccard.adonis$data <- "bin"
  ilevels.df.jaccard.adonis$source <- "allSpecies"
  ilevels.df.jaccard.adonis$features <- ncol(ilevels.df)
  
  #put all together + save
  ilevels.df.jaccard.adonis.all <- rbind(ilevels.df.jaccard.adonis, subdf.features.df.jaccard.adonis)
  adonis.bin[["INSHORE_OFFSHORE"]] <- ilevels.df.jaccard.adonis.all
}

##Get the indval_pval table
indvalout.bin.sub <- lapply(indvalout.bin, function(x){x[["indval_pval"]]})
for(i in names(indvalout.bin.sub))
{
  indvalout.bin.sub[[i]][,"comparison"] <- i
}
indvalout.bin.sub <- do.call("rbind", indvalout.bin.sub)
indvalout.bin.sub$source <- "indval"
indvalout.bin.sub$data <- "bin"


##### Save the datasets
save_path='/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/'
save(adonis.bin, #adonis results all comparisons; binary data
     adonis.maxn, # adonis results all comparisons; abudnance data
     indvalout.bin, # indval results all comparisons, abundance data
     indvalout.maxn, 
     indvalout.bin.sub, 
     indvalout.maxn.sub, file = paste0(save_path, "Indval_all_analyses_final.Rda"))