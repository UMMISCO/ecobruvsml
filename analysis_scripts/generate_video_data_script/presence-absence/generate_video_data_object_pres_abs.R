if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("reshape2")) install.packages("reshape2"); library(reshape2)
if (!require("vegan")) install.packages("vegan"); library(vegan)
if (!require("stats")) install.packages("stats"); library(stats)
if (!require("labdsv")) install.packages("labdsv"); library(labdsv)
if (!require("indicspecies")) install.packages("indicspecies"); library(indicspecies)
if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if (!require("ggrepel")) install.packages("ggrepel"); library(ggrepel)

# setwd("D:/Florian/Travail/Stage_Fonds_Meubles_2016/Papier_NC_Soft_bottoms/Analyses_Redone/04_Random_Forest")
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

sm = list()
video_data=read.delim("/data/projects/aime/data/metabarcoding/AllDataexportefixed.txt", header = TRUE, sep = "\t", quote="\"", dec=".",fill = TRUE)
head(video_data)

# Changing names in data
# mydata$G_SP=ess$SP
video_data$G_SP=gsub("_", " ",video_data$G_SP)
video_data$Transect=paste(substr(video_data$Transect,3,7))
video_data$Site=factor(video_data$Site,levels=c("ABAR","ALAG","ABAI","MBAR","MLAG","MBAI"))
video_data$Site <- do.call(rbind, lapply(video_data$Site, gsub, pattern = "BAI", replacement = "BAY"))[,1] 
video_data$Station <- do.call(rbind, lapply(video_data$Station, gsub, pattern = "BAI", replacement = "BAY"))[,1] 
head(video_data)

sm$db_long = video_data

#summarize data by station and species
matablepivot=video_data %>%
  group_by(Site,G_SP)
matablepivot=as.data.frame(matablepivot)
matablepivot=matablepivot[c("Site","Station","G_SP","MaxN")]
head(matablepivot)

# Species community matrix
abund = dcast(matablepivot, Station + Site ~ G_SP, value.var = "MaxN")
abund[is.na(abund)] = 0    # fill with zeroes for summarise below
rownames(abund)=abund$Station
# Recreate habitat because lazy
pres_abs <- abund
pres_abs[, -c(1,2)] = decostand(pres_abs[, -c(1,2)], method = "pa")
sm$X = t(pres_abs[, -c(1,2)])
pres_abs$hab <- substring(pres_abs$Site, 2)

# redefine habitat in two classes 'offshore' and 'inshore
pres_abs$Habitat <- ifelse(pres_abs$hab == "BAR", "Offshore","Inshore")

pres_abs$Zone <- ifelse(pres_abs$hab == "BAR", "Barrier", ifelse(pres_abs$hab == "BAY", "Bay", "Lagoon"))

pres_abs$Spygen = pres_abs$Station

sm$sample_info = pres_abs

# save data object
save(sm, file = "/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/video_data_object_pres_abs.rda")

