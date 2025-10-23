library(ggplot2) # for plots
library(gridExtra) # for multiple plots
library(patchwork) # to display multiples graphs in the same figure
library(dplyr)
library(ggVennDiagram)
library(reshape2)
library(ggpubr)
library(stringr)

root_path="/data/projects/aime/analyses/bruvs/FigureS5/results_replicas_analysis"

# Get a list of all directories within the root_path recursively
all_dirs <- list.dirs(path = root_path, full.names = TRUE, recursive = TRUE)

# Filter directories that contain "output_data" in their names
output_data_dirs <- grep("output_data", all_dirs, value = TRUE)

files <- list.files(output_data_dirs, recursive = TRUE)

# alldatalist.video_split_replicas <- list()
# 
# for(f in files)
# {
#   # Extract method and prevalence_rate from result file name
#   matches <- str_match(f, "^([^_]+).*_prev_(\\d+)\\.Rda$")
#   
#   if (!is.na(matches[1, 1])) {  
#     method <- matches[1, 2]  # First part (before first "_")
#     prevalence_rate <- paste0("prev_", matches[1, 3])  # Extract prev_X
#     
#     # Get the path of each result
#     file_dir_path <- paste(root_path, paste(method, "output_data", sep="_"), sep="/")
#     
#     # Load results data
#     fobj <- load(paste(file_dir_path, f, sep = "/"))
#     
#     # Get result object
#     fobj <- mget(fobj)
#     
#     # Save result object
#     alldatalist.video_split_replicas[[prevalence_rate]][[method]] <- fobj
#   }
# }

# Initialize a named list to store all result lists
alldatalist.video_split_replicas <- list()

for(f in files)
{
  # Extract the comparison and prevalence_rate from the filename
  matches <- str_match(f, "^([^_]+).*_replicas_([^_]+_vs_[^_]+)_prev_(\\d+)\\.Rda$")
  
  if (!is.na(matches[1, 1])) { 
    method <- matches[1, 2]
    replicas_comp <- matches[1, 3]  # e.g., "ABORE_vs_MBERE"
    prevalence_rate <- paste0("prev_", matches[1, 4])  # e.g., "prev_0"
    
    # Ensure the sublist exists
    if (is.null(alldatalist.video_split_replicas[[replicas_comp]])) {
      alldatalist.video_split_replicas[[replicas_comp]] <- list()
    }
    
    # Get the path of each result
    file_dir_path <- paste(root_path, paste(method, "output_data", sep="_"), sep="/")
    
    # Load data
    fobj_name <- load(paste(file_dir_path, f, sep = "/"))
    fobj <- mget(fobj_name)
    
    # Save to the dynamically-named sublist
    alldatalist.video_split_replicas[[replicas_comp]][[prevalence_rate]][[method]] <- fobj
  }
}

# # save Indval/Predomics analysis results data
# save(alldatalist.video_split_replicas, file = paste(root_path,"all_analysis_results_split_replicas.Rda", sep = "/"))
# # load result list of split by replicas analysis (Indval/terinter/bininter)
# load("/data/projects/aime/analyses/bruvs/results_replicas_analysis/all_analysis_results_split_replicas.Rda")

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

for (i in names(alldatalist.video_split_replicas)) {
  for (j in c("prev_0","prev_10")) {
    for (k in c("Indval","bininter","terinter")) {
      
      ### get object for adonis results in all comparison from each algorithm
      if(k=="Indval"){
        
        # bind results for all comparisons in abundance and pres/abs
        adonis.list[[i]][[j]][[k]][['presAbs']] <- bind_rows(alldatalist.video_split_replicas[[i]][[j]][[k]][["adonis.bin"]])
        adonis.list[[i]][[j]][[k]][['maxN']] <- bind_rows(alldatalist.video_split_replicas[[i]][[j]][[k]][["adonis.maxn"]])
        
        # get the list of species used for training
        #get colnames to remove pval and match with predomics table
        colnames.indvalout.bin <- colnames(alldatalist.video_split_replicas[[i]][[j]][[k]][["indvalout.bin.sub"]])
        colnames.indvalout.maxn <- colnames(alldatalist.video_split_replicas[[i]][[j]][[k]][["indvalout.maxn.sub"]])
        
        all_species.list[[i]][[j]][[k]][['presAbs']] <- alldatalist.video_split_replicas[[i]][[j]][[k]][["indvalout.bin.sub"]][ , !(colnames.indvalout.bin %in% "pval")]
        all_species.list[[i]][[j]][[k]][['maxN']] <- alldatalist.video_split_replicas[[i]][[j]][[k]][["indvalout.maxn.sub"]][ , !(colnames.indvalout.maxn %in% "pval")]
        
      }
      if(k %in% c("bininter","terinter")){
        
        # bind results for all comparisons in abundance and pres/abs
        adonis.list[[i]][[j]][[k]][['presAbs']] <- bind_rows(alldatalist.video_split_replicas[[i]][[j]][[k]][["adonis_pred.bin"]])
        adonis.list[[i]][[j]][[k]][['maxN']] <- bind_rows(alldatalist.video_split_replicas[[i]][[j]][[k]][["adonis_pred.maxn"]])
        
        # get the list of species used for training
        all_species.list[[i]][[j]][[k]][['presAbs']] <- alldatalist.video_split_replicas[[i]][[j]][[k]][["predout.bin.sub"]]
        all_species.list[[i]][[j]][[k]][['maxN']] <- alldatalist.video_split_replicas[[i]][[j]][[k]][["predout.maxn.sub"]]
      }
      
      # bind data by rows
      all_species.list[[i]][[j]][[k]] <- bind_rows(all_species.list[[i]][[j]][[k]])
      adonis.list[[i]][[j]][[k]] <- bind_rows(adonis.list[[i]][[j]][[k]])
      
      # add replica_comp and prev_rate to adonis table
      adonis.list[[i]][[j]][[k]][['prev_rate']] <- j
      adonis.list[[i]][[j]][[k]][['replica_comp']] <- i
      
      all_species.list[[i]][[j]][[k]][['prev_rate']] <- j
      all_species.list[[i]][[j]][[k]][['replica_comp']] <- i
      
    }
    
  }
}

# Combine adonis list into a dataframe
adonis.df <- flatten_and_bind(adonis.list)

# Rename the "Pr..F." column of adonis.df
colnames(adonis.df)[colnames(adonis.df) == "Pr..F."] <- "pvalue"

# Replace all occurrences of the value 'allSpecies_terinter' or 'allSpecies_bininter' with 'allSpecies' in adonis.df
adonis.df$source[adonis.df$source %in% c('allSpecies_terinter', 'allSpecies_bininter')] <- 'allSpecies'

# compute the ratio between r2 of individual sites comparison and r2 of the community
adonis.df <- adonis.df %>%
  group_by(comparison, data) %>%
  mutate(
    ratio = round(R2 / R2[source == "allSpecies"], digits = 3)
  ) %>%
  ungroup()

# # reorder seed from 01, 02, to 10
# adonis.df$seed= ifelse(adonis.df$seed == "seed_10", "seed_10", gsub("seed_", "seed_0", adonis.df$seed))

## Combine indics_species list into a dataframe
all_species.df <- flatten_and_bind(all_species.list)

# # rename seed column from 1 to seed_01, 10 to seed_10
# all_species.df$seed= ifelse(all_species.df$seed == "10", "seed_10", paste0("seed_0",all_species.df$seed))

# Remove row names of all_species.list
rownames(all_species.df) <- NULL

# get only indicator species
indics_species.df= all_species.df[all_species.df$IsIndSp==1,]

# Plot permanova results
adonis.df$colour <- ifelse(adonis.df$pvalue<=0.001,"blue","red")


# save permanova results for analyses split by replicas
write.csv(adonis.df, file = paste(root_path, "permava_split_by_replicas_all_results.csv", sep="/"), row.names = FALSE)

# save analysis results for all species or only indics species with Indval, Predomics_bininter, Predomics_terinter
write.csv(all_species.df, file = paste(root_path, "all_species_analysis_split_by_replicas.csv", sep="/"), row.names = FALSE)
write.csv(indics_species.df, file = paste(root_path, "indics_species_analysis_split_by_replicas.csv", sep="/"), row.names = FALSE)
