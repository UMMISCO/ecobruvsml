# =============================================================================
# Script name: save_generalization_results_with_twinspan.R
# Author: Estephe Kana, Eugeni Belda, Edi Prifti
# Date created: 2025-01-01
# Purpose: Aggregate and save generalization performance results across all methods (IndVal, Predomics, TWINSPAN, RF, SVM) for the 80/20 split analyses.
# Inputs:  analyses/analysis_outputs/all_analysis_results_with_twinspan_80_20.Rda
# Outputs: analyses/analysis_outputs/generalization_results_with_twinspan.Rda
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script:
# cd /home/estephe/Documents/Papers/2026-aime_key_fish_bruvs/Bruvs-analysis-code/analyses/scripts
# Rscript save_generalization_results_with_twinspan.R
# -----------------------------------------------------------------------------

library(ggplot2) # for plots
library(gridExtra) # for multiple plots
library(patchwork) # to display multiples graphs in the same figure
library(dplyr)
library(ggVennDiagram)
library(reshape2)
library(ggpubr)
library(stringr)

# load result list of 80_20 analysis (Indval/terinter/bininter)
# Detect project root automatically (works for any user/machine after cloning)
if (!require("here", quietly = TRUE)) install.packages("here")
base_dir <- here::here()

load(file.path(base_dir, "analyses", "analysis_outputs", "all_analysis_results_with_twinspan_80_20.Rda"))

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
# get permanova results and analysis results by Indval, predomics_ter and predomics_bin
adonis.list <- list()
all_species.list <- list()

for (i in names(alldatalist.video_80_20)) {
  for (j in c("prev_0", "prev_10")) {
    methods <- intersect(c("Indval", "bininter", "terinter", "Twinspan"), names(alldatalist.video_80_20[[i]][[j]]))
    for (k in methods) {
      ### get object for adonis results in all comparison from each algorithm
      if(k=="Indval"){
        
        # bind results for all comparisons in abundance and pres/abs
        adonis.list[[i]][[j]][[k]][['presAbs']] <- bind_rows(alldatalist.video_80_20[[i]][[j]][[k]][["adonis.bin"]])
        adonis.list[[i]][[j]][[k]][['maxN']] <- bind_rows(alldatalist.video_80_20[[i]][[j]][[k]][["adonis.maxn"]])
        
        # get the list of species used for training
        #get colnames to remove pval and match with predomics table
        colnames.indvalout.bin <- colnames(alldatalist.video_80_20[[i]][[j]][[k]][["indvalout.bin.sub"]])
        colnames.indvalout.maxn <- colnames(alldatalist.video_80_20[[i]][[j]][[k]][["indvalout.maxn.sub"]])
        
        all_species.list[[i]][[j]][[k]][['presAbs']] <- alldatalist.video_80_20[[i]][[j]][[k]][["indvalout.bin.sub"]][ , !(colnames.indvalout.bin %in% "pval")]
        all_species.list[[i]][[j]][[k]][['maxN']] <- alldatalist.video_80_20[[i]][[j]][[k]][["indvalout.maxn.sub"]][ , !(colnames.indvalout.maxn %in% "pval")]
        
        
      }
      if(k %in% c("bininter","terinter")){
        
        # bind results for all comparisons in abundance and pres/abs
        adonis.list[[i]][[j]][[k]][['presAbs']] <- bind_rows(alldatalist.video_80_20[[i]][[j]][[k]][["adonis_pred.bin"]])
        adonis.list[[i]][[j]][[k]][['maxN']] <- bind_rows(alldatalist.video_80_20[[i]][[j]][[k]][["adonis_pred.maxn"]])
        
        # get the list of species used for training
        all_species.list[[i]][[j]][[k]][['presAbs']] <- alldatalist.video_80_20[[i]][[j]][[k]][["predout.bin.sub"]]
        all_species.list[[i]][[j]][[k]][['maxN']] <- alldatalist.video_80_20[[i]][[j]][[k]][["predout.maxn.sub"]]
      }
      if(k=="Twinspan"){
        adonis.list[[i]][[j]][[k]][['presAbs']] <- bind_rows(alldatalist.video_80_20[[i]][[j]][[k]][["adonis.bin"]])
        adonis.list[[i]][[j]][[k]][['presAbs']]$source <- gsub("allSpecies","allSpecies_twinspan", adonis.list[[i]][[j]][[k]][['presAbs']]$source)
        adonis.list[[i]][[j]][[k]][['maxN']] <- bind_rows(alldatalist.video_80_20[[i]][[j]][[k]][["adonis.maxn"]])
        adonis.list[[i]][[j]][[k]][['maxN']]$source <- gsub("allSpecies","allSpecies_twinspan", adonis.list[[i]][[j]][[k]][['maxN']]$source)
        all_species.list[[i]][[j]][[k]][['presAbs']] <- alldatalist.video_80_20[[i]][[j]][[k]][["twinspan.bin.sub"]][ , !(colnames(alldatalist.video_80_20[[i]][[j]][[k]][["twinspan.bin.sub"]]) %in% "pval")]
        all_species.list[[i]][[j]][[k]][['maxN']] <- alldatalist.video_80_20[[i]][[j]][[k]][["twinspan.maxn.sub"]][ , !(colnames(alldatalist.video_80_20[[i]][[j]][[k]][["twinspan.maxn.sub"]]) %in% "pval")]
      }
      
      # bind data by rows
      all_species.list[[i]][[j]][[k]] <- bind_rows(all_species.list[[i]][[j]][[k]])
      adonis.list[[i]][[j]][[k]] <- bind_rows(adonis.list[[i]][[j]][[k]])
      
      # add seed and prev_rate to adonis table
      adonis.list[[i]][[j]][[k]][['seed']] <- i
      adonis.list[[i]][[j]][[k]][['prev_rate']] <- j
      
      all_species.list[[i]][[j]][[k]][['prev_rate']] <- j
      
    }
  }
}

# Combine adonis list into a dataframe
adonis.df <- flatten_and_bind(adonis.list)

# Rename the "Pr..F." column of adonis.df
colnames(adonis.df)[colnames(adonis.df) == "Pr..F."] <- "pvalue"

# Replace all occurrences of the value 'allSpecies_terinter' or 'allSpecies_bininter' with 'allSpecies' in adonis.df
adonis.df$source[adonis.df$source %in% c('allSpecies_terinter', 'allSpecies_bininter', 'allSpecies_twinspan')] <- 'allSpecies'

# compute the ratio between r2 of individual sites comparison and r2 of the community
adonis.df <- adonis.df %>%
  group_by(comparison, data) %>%
  mutate(
    ratio = round(R2 / R2[source == "allSpecies"], digits = 3)
  ) %>%
  ungroup()

# reorder seed from 01, 02, to 10
adonis.df$seed= ifelse(adonis.df$seed == "seed_10", "seed_10", gsub("seed_", "seed_0", adonis.df$seed))

## Combine indics_species list into a dataframe
all_species.df <- flatten_and_bind(all_species.list)

# rename seed column from 1 to seed_01, 10 to seed_10
all_species.df$seed= ifelse(all_species.df$seed == "10", "seed_10", paste0("seed_0",all_species.df$seed))

# Remove row names of all_species.list
rownames(all_species.df) <- NULL

# get only indicator species
indics_species.df= all_species.df[all_species.df$IsIndSp==1,]

# Plot permanova results
adonis.df$colour <- ifelse(adonis.df$pvalue<=0.001,"blue","red")


# save permanova results for analyses in 80_20
root_path <- file.path(base_dir, "analyses", "analysis_outputs")
# write.csv(adonis.df, file = paste(root_path, "permava_80_20_all_results.csv", sep="/"), row.names = FALSE)

# save analysis results for all species or only indics species with Indval, Predomics_bininter, Predomics_terinter
write.csv(all_species.df, file = paste(root_path, "all_species_analysis_results_80_20_with_twinspan.csv", sep="/"), row.names = FALSE)
write.csv(indics_species.df, file = paste(root_path, "indics_species_analysis_results_80_20_with_twinspan.csv", sep="/"), row.names = FALSE)