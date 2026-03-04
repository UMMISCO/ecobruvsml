# =============================================================================
# Script name: format_results_split_by_replicas_with_twinspan.R
# Author: Estephe Kana, Eugeni Belda, Edi Prifti
# Date created: 2025-01-01
# Purpose: Aggregate and format performance metrics from replica-split analyses
#          including TWINSPAN for Figure S5.
# Inputs:  analyses/analysis_outputs/replicas_analysis_results/ (multiple .Rda)
#          analyses/analysis_outputs/rf_analysis_by_replicas/ (multiple .Rda)
# Outputs: analyses/analysis_outputs/replicas_analysis_results/
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script:
# cd /home/estephe/Documents/Papers/2026-aime_key_fish_bruvs/Bruvs-analysis-code/analyses/Figures/FigureS5
# Rscript format_results_split_by_replicas_with_twinspan.R
# -----------------------------------------------------------------------------


library(ggplot2)
library(gridExtra)
library(patchwork)
library(dplyr)
library(ggVennDiagram)
library(reshape2)
library(ggpubr)
library(stringr)

# Detect project root automatically (works for any user/machine after cloning)
if (!require("here", quietly = TRUE)) install.packages("here")
base_dir <- here::here()
root_path <- file.path(base_dir, "analyses", "Figures", "FigureS5", "results_replicas_analysis")

# ── 1. Discover output files ──────────────────────────────────────────────────
all_dirs         <- list.dirs(path = root_path, full.names = TRUE, recursive = TRUE)
output_data_dirs <- grep("output_data", all_dirs, value = TRUE)
files            <- list.files(output_data_dirs, recursive = TRUE)

# ── 2. Load all result files into a nested list ───────────────────────────────
# Structure: alldatalist[[replicas_comp]][[prev_rate]][[method]]
alldatalist.video_split_replicas <- list()

for (f in files) {
  # Match pattern: <Method>_all_*_replicas_<R1_vs_R2>_prev_<N>.Rda
  matches <- str_match(f, "^([^_]+).*_replicas_([^_]+_vs_[^_]+)_prev_(\\d+)\\.Rda$")
  
  if (!is.na(matches[1, 1])) {
    method          <- matches[1, 2]                      # e.g. "Indval", "bininter", "terinter", "Twinspan"
    replicas_comp   <- matches[1, 3]                      # e.g. "ABORE_vs_MBERE"
    prevalence_rate <- paste0("prev_", matches[1, 4])     # e.g. "prev_0", "prev_10"
    
    if (is.null(alldatalist.video_split_replicas[[replicas_comp]]))
      alldatalist.video_split_replicas[[replicas_comp]] <- list()
    
    file_dir_path <- paste(root_path, paste(method, "output_data", sep = "_"), sep = "/")
    fobj_name     <- load(paste(file_dir_path, f, sep = "/"))
    fobj          <- mget(fobj_name)
    
    alldatalist.video_split_replicas[[replicas_comp]][[prevalence_rate]][[method]] <- fobj
  }
}

# save Indval/Predomics analysis results data
save(alldatalist.video_split_replicas, file = paste(root_path,"all_analysis_results_split_replicas_with_twinspan.Rda", sep = "/"))

# ── 3. Helper ─────────────────────────────────────────────────────────────────
flatten_and_bind <- function(nested_list) {
  if (all(sapply(nested_list, is.data.frame))) {
    do.call(rbind, nested_list)
  } else {
    do.call(rbind, lapply(nested_list, flatten_and_bind))
  }
}

# ── 4. Extract adonis + species tables for all methods ────────────────────────
# Methods and the object-name conventions they use:
#
#  Indval    : adonis.bin / adonis.maxn        | indvalout.bin.sub / indvalout.maxn.sub
#  bininter  : adonis_pred.bin / .maxn         | predout.bin.sub   / predout.maxn.sub
#  terinter  : adonis_pred.bin / .maxn         | predout.bin.sub   / predout.maxn.sub
#  Twinspan  : adonis.bin / adonis.maxn        | twinspan.bin.sub  / twinspan.maxn.sub

adonis.list      <- list()
all_species.list <- list()

for (i in names(alldatalist.video_split_replicas)) {
  for (j in c("prev_0", "prev_10")) {
    for (k in c("Indval", "bininter", "terinter", "Twinspan")) {
      
      # Skip silently if this method/prev combination was not loaded
      if (is.null(alldatalist.video_split_replicas[[i]][[j]][[k]])) next
      
      obj <- alldatalist.video_split_replicas[[i]][[j]][[k]]
      
      # ---- Indval ----
      if (k == "Indval") {
        adonis.list[[i]][[j]][[k]][["presAbs"]] <- bind_rows(obj[["adonis.bin"]])
        adonis.list[[i]][[j]][[k]][["maxN"]]    <- bind_rows(obj[["adonis.maxn"]])
        
        cols.bin  <- colnames(obj[["indvalout.bin.sub"]])
        cols.maxn <- colnames(obj[["indvalout.maxn.sub"]])
        all_species.list[[i]][[j]][[k]][["presAbs"]] <- obj[["indvalout.bin.sub"]][ , !(cols.bin  %in% "pval")]
        all_species.list[[i]][[j]][[k]][["maxN"]]    <- obj[["indvalout.maxn.sub"]][ , !(cols.maxn %in% "pval")]
      }
      
      # ---- Predomics (bininter / terinter) ----
      if (k %in% c("bininter", "terinter")) {
        adonis.list[[i]][[j]][[k]][["presAbs"]] <- bind_rows(obj[["adonis_pred.bin"]])
        adonis.list[[i]][[j]][[k]][["maxN"]]    <- bind_rows(obj[["adonis_pred.maxn"]])
        
        all_species.list[[i]][[j]][[k]][["presAbs"]] <- obj[["predout.bin.sub"]]
        all_species.list[[i]][[j]][[k]][["maxN"]]    <- obj[["predout.maxn.sub"]]
      }
      
      # ---- Twinspan ----
      if (k == "Twinspan") {
        adonis.list[[i]][[j]][[k]][["presAbs"]] <- bind_rows(obj[["adonis.bin"]])
        adonis.list[[i]][[j]][[k]][["presAbs"]]$source    <- gsub("allSpecies", "allSpecies_twinspan", adonis.list[[i]][[j]][[k]][["presAbs"]]$source)
        adonis.list[[i]][[j]][[k]][["maxN"]]    <- bind_rows(obj[["adonis.maxn"]])
        adonis.list[[i]][[j]][[k]][["maxN"]]$source    <- gsub("allSpecies", "allSpecies_twinspan", adonis.list[[i]][[j]][[k]][["maxN"]]$source)
        
        all_species.list[[i]][[j]][[k]][["presAbs"]] <- obj[["twinspan.bin.sub"]]
        all_species.list[[i]][[j]][[k]][["maxN"]]    <- obj[["twinspan.maxn.sub"]]
      }
      
      # ---- Bind data types and annotate ----
      all_species.list[[i]][[j]][[k]] <- bind_rows(all_species.list[[i]][[j]][[k]])
      adonis.list[[i]][[j]][[k]]      <- bind_rows(adonis.list[[i]][[j]][[k]])
      
      adonis.list[[i]][[j]][[k]][["prev_rate"]]    <- j
      adonis.list[[i]][[j]][[k]][["replica_comp"]] <- i
      
      all_species.list[[i]][[j]][[k]][["prev_rate"]]    <- j
      all_species.list[[i]][[j]][[k]][["replica_comp"]] <- i
    }
  }
}

# ── 5. Flatten to data frames ─────────────────────────────────────────────────
adonis.df <- flatten_and_bind(adonis.list)
colnames(adonis.df)[colnames(adonis.df) == "Pr..F."] <- "pvalue"

# remove NAs
adonis.df <- na.omit(adonis.df)

# Harmonise "allSpecies" source label across methods
adonis.df$source[adonis.df$source %in% c("allSpecies_terinter", "allSpecies_bininter", "allSpecies_twinspan")] <- "allSpecies"

# Compute R² ratio relative to full-community model
adonis.df <- adonis.df %>%
  group_by(comparison, data, replica_comp, prev_rate) %>%
  mutate(
    ratio = round(R2 / unique(R2[source == "allSpecies"]), digits = 3)
  ) %>%
  ungroup()

all_species.df <- flatten_and_bind(all_species.list)
rownames(all_species.df) <- NULL

# Indicator species only
indics_species.df <- all_species.df[all_species.df$IsIndSp == 1, ]

# Colour flag for significance
adonis.df$colour <- ifelse(adonis.df$pvalue <= 0.001, "blue", "red")

# ── 6. Save outputs ───────────────────────────────────────────────────────────
write.csv(adonis.df,
          file = paste(root_path, "permanova_split_by_replicas_all_results_with_twinspan.csv", sep = "/"),
          row.names = FALSE)

write.csv(all_species.df,
          file = paste(root_path, "all_species_analysis_split_by_replicas_with_twinspan.csv", sep = "/"),
          row.names = FALSE)

write.csv(indics_species.df,
          file = paste(root_path, "indics_species_analysis_split_by_replicas_with_twinspan.csv", sep = "/"),
          row.names = FALSE)

message("Done. Results written to: ", root_path)
