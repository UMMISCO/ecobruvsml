# =============================================================================
# Script name: codeTwinspan_all_analyses_prev.R
# Author: Estephe Kana, Eugeni Belda, Edi Prifti
# Date created: 2025-01-01
# Purpose: Perform Two-Way Indicator Species Analysis (TWINSPAN) with prevalence filtering for all pairwise habitat comparisons and INSHORE vs OFFSHORE.
# Inputs:  data/AllDataexportefixed.txt; args: prevalence_rate (%)
# Outputs: analyses/analysis_outputs/twinspan_output_data/
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script:
# cd /home/estephe/Documents/Papers/2026-aime_key_fish_bruvs/Bruvs-analysis-code/analyses/scripts
# Rscript codeTwinspan_all_analyses_prev.R 10
# -----------------------------------------------------------------------------


if (!require("dplyr"))      install.packages("dplyr");      library(dplyr)
if (!require("reshape2"))   install.packages("reshape2");   library(reshape2)
if (!require("vegan"))      install.packages("vegan");      library(vegan)
if (!require("stringr"))    install.packages("stringr");    library(stringr)
if (!require("twinspanR"))  install.packages("twinspanR"); library(twinspanR)

# ─────────────────────────────────────────────────────────────────────────────
# PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
prevalence_rate <- as.numeric(args[1])
# prevalence_rate <- 10   # % – species must occur in at least this % of samples

# Detect project root automatically (works for any user/machine after cloning)
if (!require("here", quietly = TRUE)) install.packages("here")
base_dir <- here::here()

dataset <- file.path(base_dir, "data", "AllDataexportefixed.txt")

# ─────────────────────────────────────────────────────────────────────────────
# DATA LOADING & CLEANING
# ─────────────────────────────────────────────────────────────────────────────
mydata <- read.delim(dataset, header = TRUE, sep = "\t",
                     quote = "\"", dec = ".", fill = TRUE)

mydata$G_SP     <- gsub("_", " ", mydata$G_SP)
mydata$Transect <- paste(substr(mydata$Transect, 3, 7))
mydata$Site     <- factor(mydata$Site,
                          levels = c("ABAR","ALAG","ABAI","MBAR","MLAG","MBAI"))
# mydata$Site     <- do.call(rbind, lapply(mydata$Site,    gsub,
#                                          pattern = "BAI", replacement = "BAY"))[,1]
mydata$Site <- gsub("BAI", "BAY", mydata$Site)
# mydata$Station  <- do.call(rbind, lapply(mydata$Station, gsub,
#                                          pattern = "BAI", replacement = "BAY"))[,1]
mydata$Station <- gsub("BAI", "BAY", mydata$Station)

# Summarise to station × species
matablepivot <- mydata[, c("Site","Station","G_SP","MaxN")]

# Build species community matrix (Station × Species)
abund <- reshape2::dcast(matablepivot, Station + Site ~ G_SP,
                         value.var = "MaxN")
abund[is.na(abund)] <- 0
rownames(abund) <- abund$Station

# Derive habitat from Site label
abund$hab <- substring(abund$Site, 2)
acc <- abund

# Inshore / offshore grouping
acc$hab_inoff <- ifelse(acc$hab %in% c("BAY","LAG"), "INSHORE",
                        ifelse(acc$hab == "BAR",      "OFFSHORE", NA))

table(acc$hab, acc$hab_inoff)

# All pairwise habitat combinations
comp <- combn(x = unique(acc$hab), m = 2, simplify = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

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


## Function to parse and get class of all the species
parse_species_class <- function(twi_lines) {
  
  # Ensure input is a plain character vector (tw$twi[[1]] can be a list)
  twi_lines <- as.character(unlist(twi_lines))
  
  # ── Locate "ORDER OF SAMPLES" as the start of the tabulation block ─────────
  tab_start <- which(str_detect(twi_lines, "ORDER OF SAMPLES"))
  if (length(tab_start) == 0)
    stop("Could not find 'ORDER OF SAMPLES' in TWINSPAN output.")
  
  # ── Find the classification bar: a line of only 0s and 1s (+ whitespace) ───
  tail_lines   <- twi_lines[(tab_start + 1):length(twi_lines)]
  bar_rel_idx  <- which(str_detect(tail_lines, "^\\s+[01]+\\s*$"))
  if (length(bar_rel_idx) == 0)
    stop("Could not find the species classification bar (0/1 line) in TWINSPAN output.")
  bar_abs_idx  <- bar_rel_idx[1] + tab_start   # absolute index in twi_lines
  
  # ── Extract species lines between tab_start and the bar ───────────────────
  # A valid species line: leading whitespace, integer, space, "Abcd efgh",
  # abundance pattern, trailing whitespace, then a single 0 or 1.
  candidate_lines <- twi_lines[(tab_start + 1):(bar_abs_idx - 1)]
  species_lines   <- candidate_lines[
    str_detect(candidate_lines,
               "^\\s+\\d+\\s+[A-Z][a-z]+ [a-z]+\\s+\\S+\\s+[01]\\s*$")]
  
  if (length(species_lines) == 0)
    stop("No species lines matched in the tabulation block.")
  
  # ── Parse each line ────────────────────────────────────────────────────────
  parsed <- lapply(species_lines, function(line) {
    # Trailing class: last non-space token (0 or 1)
    sp_class  <- str_trim(str_extract(line, "[01]\\s*$"))
    # Species abbreviation: "Xxxx xxxx" immediately after the leading integer
    sp_abbrev <- str_match(line,
                           "^\\s*\\d+\\s+([A-Z][a-z]+ [a-z]+)")[1, 2]
    data.frame(abbrev.name    = sp_abbrev,
               twinspan_class = paste0("*", sp_class),
               stringsAsFactors = FALSE)
  })
  
  do.call(rbind, parsed)
}

## TWINSPAN analysis ──────────────────────────────────────────────────────────
run_twinspan_analysis <- function(sample.X, sample.info, ilevels,
                                  hab_type, data_type = "maxN") {
  
  # ── Optionally binarise ────────────────────────────────────────────────────
  if (data_type == "pres/abs") {
    X_input <- as.data.frame(apply(sample.X, 2,
                                   function(x) ifelse(x == 0, 0, 1)))
    dist <- "jaccard"
  } else {
    X_input <- sample.X
    dist    <- "bray"
  }
  
  # ── Run TWINSPAN ───────────────────────────────────────────────────────────
  tw <- twinspan(X_input, modif = TRUE, clusters = 2, diss = dist)
  
  # ── Parse TWINSPAN indicator species (abbreviated names + group) ───────────
  twi_lines     <- as.character(unlist(tw$twi[[1]]))
  indicator_idx <- which(trimws(twi_lines) == "INDICATORS, together with their SIGN")
  if (length(indicator_idx) == 0)
    stop("Could not find 'INDICATORS, together with their SIGN' in TWINSPAN output.")
  indicator_line <- twi_lines[indicator_idx + 1]
  
  indicators_raw <- str_trim(unlist(str_split(indicator_line, "\\s{2,}")))
  indicators_raw <- indicators_raw[indicators_raw != ""]
  
  ind_parsed <- data.frame(species             = character(),
                           pseudospecies_level = character(),
                           group               = character(),
                           stringsAsFactors    = FALSE)
  for (ind in indicators_raw) {
    if (!grepl("\\(", ind)) next
    m <- str_match(ind, "(.+?)(\\d+)\\(([+-])\\)")
    if (!is.na(m[1, 1])) {
      ind_parsed <- rbind(ind_parsed,
                          data.frame(
                            species             = gsub(" ", "", str_trim(m[1, 2])),
                            pseudospecies_level = m[1, 3],
                            group               = ifelse(m[1, 4] == "+", "*1", "*0"),
                            stringsAsFactors    = FALSE))
    }
  }
  # Full names for the indicator set (used only for IsIndSp flagging)
  ind_fullnames <- character(0)
  if (nrow(ind_parsed) > 0) {
    ind_merged    <- merge(ind_parsed, tw$spnames,
                           by.x = "species", by.y = "abbrev.name", all.x = TRUE)
    ind_fullnames <- na.omit(ind_merged$full.name)
  }
  
  # ── Cluster membership per plot ────────────────────────────────────────────
  tw_classif <- tw$classif
  tw_group   <- merge(tw_classif,
                      sample.info[, c("Station", "Site", hab_type)],
                      by.x = "plot.no", by.y = 0, all.x = TRUE)
  
  # ── Habitat label per TWINSPAN class (majority vote) ──────────────────────
  class_hab <- tw_group %>%
    group_by(class, .data[[hab_type]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(class) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(class, all_of(hab_type)) %>%
    as.data.frame()
  
  class_purity <- tw_group %>%
    group_by(class) %>%
    summarise(n_hab = n_distinct(.data[[hab_type]]), .groups = "drop")
  mixed <- class_purity$class[class_purity$n_hab > 1]
  if (length(mixed) > 0)
    warning("TWINSPAN class(es) ", paste(mixed, collapse = ", "),
            " contain samples from multiple habitats — ",
            "habitat label assigned by majority vote.")
  
  # ── All-species TWINSPAN class from ordered tabulation ────────────────────
  all_sp_class <- parse_species_class(twi_lines)        # abbrev.name + twinspan_class
  all_sp_class$abbrev.name <-  gsub(" ", "", all_sp_class$abbrev.name)
  all_sp_class <- merge(all_sp_class, tw$spnames,
                        by = "abbrev.name", all.x = TRUE)
  all_sp_class <- merge(all_sp_class,
                        class_hab,
                        by.x = "twinspan_class", by.y = "class",
                        all.x = TRUE)
  
  # ── Chi-square: species values × habitat group ────────────────────────────
  # For every species, cross-tabulate its per-sample values (MaxN counts or
  # 0/1) against the habitat label of the sample, then extract X² and p-value.
  hab_vec <- sample.info[rownames(X_input), hab_type]
  
  chisq_stats <- lapply(colnames(X_input), function(sp) {
    tbl <- table(X_input[[sp]], hab_vec)
    ct  <- suppressWarnings(chisq.test(tbl))
    data.frame(full.name         = sp,
               featureImportance = round(unname(ct$statistic), 4),
               pval              = round(ct$p.value,           6),
               stringsAsFactors  = FALSE)
  })
  chisq_df <- do.call(rbind, chisq_stats)
  
  # ── Prevalence for all species ─────────────────────────────────────────────
  prev_vec <- round((colSums(X_input > 0) / nrow(X_input)) * 100, 3)
  
  # ── Assemble final per-species output table ────────────────────────────────
  indicator_df <- all_sp_class %>%
    rename(class = all_of(hab_type)) %>%
    left_join(chisq_df, by = "full.name") %>%
    mutate(
      IsIndSp    = as.integer(full.name %in% ind_fullnames),
      prevalence = prev_vec[match(full.name, colnames(X_input))],
      comparison = paste(ilevels, collapse = "_"),
      source     = "twinspan",
      data       = data_type
    ) %>%
    select(
      feature           = full.name,
      featureImportance,
      pval,
      class,
      IsIndSp,
      prevalence,
      comparison,
      source,
      data
    ) %>%
    as.data.frame()
  
  # ── PERMANOVA ─────────────────────────────────────────────────────────────
  adonis_results <- NULL
  ind_species    <- indicator_df$feature[indicator_df$IsIndSp == 1]
  ind_species    <- ind_species[!is.na(ind_species) & ind_species %in% colnames(X_input)]
  binary_flag    <- data_type == "pres/abs"
  
  if (length(ind_species) >= 2) {
    sub_X    <- X_input[, ind_species, drop = FALSE]
    sub_X    <- sub_X[rowSums(sub_X) > 0, , drop = FALSE]
    sub_meta <- data.frame(sample = rownames(sub_X),
                           class  = sample.info[rownames(sub_X), hab_type])
    
    dist_sub <- vegdist(sub_X, method = dist, binary = binary_flag)
    set.seed(100)
    adonis_sub <- adonis2(dist_sub ~ class, data = sub_meta)
    adonis_sub <- data.frame(adonis_sub)["Model", , drop = FALSE]
    adonis_sub$comparison <- paste(ilevels, collapse = "_")
    adonis_sub$data       <- data_type
    adonis_sub$source     <- "twinspanSpecies"
    adonis_sub$features   <- length(ind_species)
    
    full_meta <- data.frame(sample = rownames(X_input),
                            class  = sample.info[rownames(X_input), hab_type])
    dist_all  <- vegdist(X_input, method = dist, binary = binary_flag)
    set.seed(100)
    adonis_all <- adonis2(dist_all ~ class, data = full_meta)
    adonis_all <- data.frame(adonis_all)["Model", , drop = FALSE]
    adonis_all$comparison <- paste(ilevels, collapse = "_")
    adonis_all$data       <- data_type
    adonis_all$source     <- "allSpecies"
    adonis_all$features   <- ncol(X_input)
    
    adonis_results <- rbind(adonis_all, adonis_sub)
  } else {
    message("  Not enough indicator species for PERMANOVA.")
  }
  
  list(twinspan_obj = tw,
       classif      = tw_group,
       indicator_df = indicator_df,
       adonis       = adonis_results)
}


# ─────────────────────────────────────────────────────────────────────────────
# MAIN LOOP – pairwise habitat comparisons
# ─────────────────────────────────────────────────────────────────────────────
twinspan.maxn  <- list()
adonis.maxn    <- list()
twinspan.bin   <- list()
adonis.bin     <- list()

for (i in seq_along(comp)) {
  
  ilevels   <- comp[[i]]
  key       <- paste(ilevels, collapse = "_")
  message("\n── Comparison: ", key, " ──────────────────────────────────────────")
  
  ilevels_df <- acc[acc$hab %in% ilevels, ]
  data_prep  <- prepare_data_for_analysis(ilevels_df, "hab", ilevels[[1]], ilevels[[2]])
  
  # ── MaxN ──────────────────────────────────────────────────────────────────
  message("  Running TWINSPAN on MaxN data …")
  res_maxn <- run_twinspan_analysis(
    sample.X    = data_prep$sample.X,
    sample.info = data_prep$sample.info,
    ilevels     = ilevels,
    hab_type    = "hab",
    data_type   = "maxN")
  
  twinspan.maxn[[key]] <- res_maxn
  if (!is.null(res_maxn$adonis)) adonis.maxn[[key]] <- res_maxn$adonis
  
  # message("  Indicator species (MaxN):")
  # print(res_maxn$indicator_df[, c("feature","group","pseudospecies_level",
  #                                 "prevalence","habitat","comparison")])
  
  # ── Presence / Absence ───────────────────────────────────────────────────
  message("  Running TWINSPAN on presence/absence data …")
  res_bin <- run_twinspan_analysis(
    sample.X    = data_prep$sample.X,
    sample.info = data_prep$sample.info,
    ilevels     = ilevels,
    hab_type    = "hab",
    data_type   = "pres/abs")
  
  twinspan.bin[[key]] <- res_bin
  if (!is.null(res_bin$adonis)) adonis.bin[[key]] <- res_bin$adonis
  
  # message("  Indicator species (pres/abs):")
  # print(res_bin$indicator_df[, c("full.name","group","pseudospecies_level",
  #                                "prevalence","habitat","comparison")])
}

# ─────────────────────────────────────────────────────────────────────────────
# INSHORE vs OFFSHORE comparison
# ─────────────────────────────────────────────────────────────────────────────
message("\n── Comparison: INSHORE_OFFSHORE ─────────────────────────────────────────")

ilevels_hab <- list("INSHORE","OFFSHORE")
key_io      <- paste(ilevels_hab, collapse = "_")

data_io <- prepare_data_for_analysis(acc, "hab_inoff",
                                     ilevels_hab[[1]], ilevels_hab[[2]])

# MaxN
message("  Running TWINSPAN on MaxN data …")
res_io_maxn <- run_twinspan_analysis(
  sample.X    = data_io$sample.X,
  sample.info = data_io$sample.info,
  ilevels     = ilevels_hab,
  hab_type    = "hab_inoff",
  data_type   = "maxN")

twinspan.maxn[[key_io]] <- res_io_maxn
if (!is.null(res_io_maxn$adonis)) adonis.maxn[[key_io]] <- res_io_maxn$adonis

# Presence / Absence
message("  Running TWINSPAN on presence/absence data …")
res_io_bin <- run_twinspan_analysis(
  sample.X    = data_io$sample.X,
  sample.info = data_io$sample.info,
  ilevels     = ilevels_hab,
  hab_type    = "hab_inoff",
  data_type   = "pres/abs")

twinspan.bin[[key_io]] <- res_io_bin
if (!is.null(res_io_bin$adonis)) adonis.bin[[key_io]] <- res_io_bin$adonis

# ─────────────────────────────────────────────────────────────────────────────
# FLATTEN INDICATOR SPECIES TABLES
# ─────────────────────────────────────────────────────────────────────────────
flatten_indicators <- function(twinspan_list, data_type_label) {
  sub <- lapply(twinspan_list, function(x) x[["indicator_df"]])
  do.call(rbind, sub) |>
    mutate(source    = "twinspan",
           data      = data_type_label)
}

twinspan.maxn.sub <- flatten_indicators(twinspan.maxn, "maxN")
twinspan.bin.sub  <- flatten_indicators(twinspan.bin,  "pres/abs")

# ─────────────────────────────────────────────────────────────────────────────
# FLATTEN PERMANOVA TABLES
# ─────────────────────────────────────────────────────────────────────────────
adonis.maxn.sub <- do.call(rbind, adonis.maxn)
adonis.bin.sub  <- do.call(rbind, adonis.bin)

# ─────────────────────────────────────────────────────────────────────────────
# SAVE OUTPUTS
# ─────────────────────────────────────────────────────────────────────────────
save_path <- file.path(base_dir, "analyses", "analysis_outputs", "twinspan_output_data")

save(
  twinspan.maxn,       # full TWINSPAN objects – MaxN, all comparisons
  twinspan.bin,        # full TWINSPAN objects – pres/abs, all comparisons
  twinspan.maxn.sub,   # flattened indicator-species table – MaxN
  twinspan.bin.sub,    # flattened indicator-species table – pres/abs
  adonis.maxn,         # PERMANOVA results – MaxN, all comparisons
  adonis.bin,          # PERMANOVA results – pres/abs, all comparisons
  adonis.maxn.sub,     # flattened PERMANOVA table – MaxN
  adonis.bin.sub,      # flattened PERMANOVA table – pres/abs
  file = file.path(save_path,
                   paste("Twinspan_all_analyses_overall_data",
                         "prev", paste0(prevalence_rate, ".Rda"),
                         sep = "_"))
)

message("\nAll done. Results saved to: ", save_path)