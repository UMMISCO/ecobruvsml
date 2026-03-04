# =============================================================================
# Script name: codeTwinspan_all_analyses_split_by_replicas.R
# Author: Estephe Kana, Eugeni Belda, Edi Prifti
# Date created: 2025-02-01
# Purpose: Run TWINSPAN analysis on BRUVS data split by replicates across all
#          habitat comparisons. Takes command-line arguments for prevalence
#          threshold and replica IDs.
# Inputs:  data/AllDataexportefixed.txt
# Outputs: Figures/FigureS5/results_replicas_analysis/twinspan_output_data/
#          (one .Rda file per replica comparison and prevalence threshold)
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script:
# cd /home/estephe/Documents/Papers/2026-aime_key_fish_bruvs/Bruvs-analysis-code/analyses/Figures/FigureS5
# Rscript codeTwinspan_all_analyses_split_by_replicas.R 10 MBERE ABORE
# -----------------------------------------------------------------------------

if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("reshape2")) install.packages("reshape2"); library(reshape2)
if (!require("vegan")) install.packages("vegan"); library(vegan)
if (!require("stringr")) install.packages("stringr"); library(stringr)
if (!require("twinspanR")) install.packages("twinspanR"); library(twinspanR)

# Detect project root automatically (works for any user/machine after cloning)
if (!require("here", quietly = TRUE)) install.packages("here")
base_dir <- here::here()

args <- commandArgs(trailingOnly = TRUE)
prevalence_rate <- as.numeric(args[1])
replica1 <- args[2]  # train models on replica1
replica2 <- args[3]  # test models on replica2

dataset <- file.path(base_dir, "data", "AllDataexportefixed.txt")
mydata <- read.delim(dataset, header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE)

mydata$G_SP <- gsub("_", " ", mydata$G_SP)
mydata$Transect <- paste(substr(mydata$Transect, 3, 7))
mydata$Site <- factor(mydata$Site, levels = c("ABAR", "ALAG", "ABAI", "MBAR", "MLAG", "MBAI"))
mydata$Site <- gsub("BAI", "BAY", mydata$Site)
mydata$Station <- gsub("BAI", "BAY", mydata$Station)

matablepivot <- mydata[, c("Site", "Station", "Transect", "G_SP", "MaxN")]
abund <- reshape2::dcast(matablepivot, Station + Transect + Site ~ G_SP, value.var = "MaxN")
abund[is.na(abund)] <- 0
rownames(abund) <- abund$Station
abund$hab <- substring(abund$Site, 2)
acc <- abund
acc$hab_inoff <- ifelse(acc$hab %in% c("BAY", "LAG"), "INSHORE", ifelse(acc$hab == "BAR", "OFFSHORE", NA))

comp <- combn(x = unique(acc$hab), m = 2, simplify = FALSE)
ilevels_hab <- list("INSHORE", "OFFSHORE")

# Function to get sample data by prevalence
get_sample_by_prevalence <- function(abundance_matrix, prevalence_rate) {
  prevalence_rate <- prevalence_rate / 100
  if (prevalence_rate < 0 || prevalence_rate > 1) stop("Prevalence rate must be between 0 and 100.")
  num_samples <- nrow(abundance_matrix)
  min_occurrences <- ceiling(prevalence_rate * num_samples)
  valid_features <- colSums(abundance_matrix > 0) >= min_occurrences
  sampled_matrix <- abundance_matrix[, valid_features, drop = FALSE]
  sampled_matrix <- sampled_matrix[rowSums(sampled_matrix > 0) > 0, , drop = FALSE]
  sampled_matrix
}

# Function to prepare data for analysis (replica-based split)
prepare_data_for_analysis <- function(sample_data, habitat_type, hab1, hab2) {
  if (hab1 == hab2) stop("It seems that the habitats you entered are the same.")
  sample_data <- sample_data[!is.na(sample_data[[habitat_type]]) &
                               sample_data[[habitat_type]] %in% c(hab1, hab2), ]
  sample_data <- sample_data[, colSums(sample_data > 0) > 0]
  X <- sample_data[, -match(c("Station", "Transect", "Site", "hab", "hab_inoff"), colnames(sample_data))]
  sample_X <- get_sample_by_prevalence(X, prevalence_rate)
  unprevalent_species <- setdiff(colnames(X), colnames(sample_X))
  sample_data <- sample_data[, !(colnames(sample_data) %in% unprevalent_species)]
  sample_data <- sample_data[rownames(sample_data) %in% rownames(sample_X), ]
  sample_data <- sample_data[rowSums(sample_data > 0) > 0, ]
  sample_data.class <- ifelse(sample_data[[habitat_type]] == hab1, -1,
                              ifelse(sample_data[[habitat_type]] == hab2, 1, NA))
  sample_data.X <- sample_data[, -match(c("Station", "Transect", "Site", "hab", "hab_inoff"), colnames(sample_data))]
  list(sample.info = sample_data,
       sample.X    = sample_data.X,
       sample.class = sample_data.class)
}

# Function to split data by replica transect labels
split_train_test_by_replica <- function(sample.info, sample.X, sample.class,
                                        replica1, replica2, hab_type = "hab") {
  train_idx <- sample.info$Transect == replica1
  test_idx  <- sample.info$Transect == replica2
  
  if (sum(train_idx) == 0) stop("No training samples found for replica: ", replica1)
  if (sum(test_idx)  == 0) stop("No test samples found for replica: ", replica2)
  
  train_classes <- unique(sample.info[[hab_type]][train_idx])
  test_classes  <- unique(sample.info[[hab_type]][test_idx])
  if (length(train_classes) < 2)
    stop("Training set (", replica1, ") has only one habitat class.")
  if (length(test_classes) < 2)
    stop("Test set (", replica2, ") has only one habitat class.")
  
  list(
    X_train    = sample.X[ train_idx, , drop = FALSE],
    y_train    = sample.class[ train_idx],
    info_train = sample.info[ train_idx, , drop = FALSE],
    X_test     = sample.X[ test_idx, , drop = FALSE],
    y_test     = sample.class[ test_idx],
    info_test  = sample.info[ test_idx, , drop = FALSE]
  )
}

parse_species_class <- function(twi_lines) {
  twi_lines <- as.character(unlist(twi_lines))
  tab_start <- which(str_detect(twi_lines, "ORDER OF SAMPLES"))
  if (length(tab_start) == 0) stop("Could not find 'ORDER OF SAMPLES' in TWINSPAN output.")
  
  tail_lines  <- twi_lines[(tab_start + 1):length(twi_lines)]
  bar_rel_idx <- which(str_detect(tail_lines, "^\\s+[01]+\\s*$"))
  if (length(bar_rel_idx) == 0) stop("Could not find the species classification bar.")
  
  bar_abs_idx     <- bar_rel_idx[1] + tab_start
  candidate_lines <- twi_lines[(tab_start + 1):(bar_abs_idx - 1)]
  
  species_lines <- candidate_lines[
    str_detect(candidate_lines,
               "^\\s+\\d+\\s+[A-Z][A-Za-z.]+(?:\\s+[a-z.]+)+\\s+\\S+\\s+[01]\\s*$")
  ]
  if (length(species_lines) == 0) stop("No species lines matched in the tabulation block.")
  
  parsed <- lapply(species_lines, function(line) {
    sp_class  <- str_trim(str_extract(line, "[01]\\s*$"))
    sp_abbrev <- str_match(line,
                           "^\\s*\\d+\\s+((?:[A-Z][A-Za-z.]+)(?:\\s+[a-z.]+)+)\\s+\\S+\\s+[01]\\s*$")[1, 2]
    data.frame(abbrev.name = sp_abbrev, twinspan_class = paste0("*", sp_class),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, parsed)
}

run_twinspan_analysis <- function(Comp_data, ilevels, hab_type, data_type = "maxN") {
  X_train    <- Comp_data$X_train
  y_train    <- Comp_data$y_train
  sample.info <- Comp_data$info_train
  
  if (data_type == "pres/abs") {
    X_input     <- as.data.frame(apply(X_train, 2, function(x) ifelse(x == 0, 0, 1)))
    dist        <- "jaccard"
    binary_flag <- TRUE
  } else {
    X_input     <- X_train
    dist        <- "bray"
    binary_flag <- FALSE
  }
  
  full_names  <- colnames(X_input)
  split_names <- strsplit(full_names, " ")
  genus       <- sapply(split_names, `[`, 1)
  species     <- sapply(split_names, `[`, 2)
  abbrev_keys <- make.unique(paste0(substr(genus, 1, 4), substr(species, 1, 4)), sep = ".")
  name_map    <- setNames(full_names, abbrev_keys)
  colnames(X_input) <- abbrev_keys
  
  tw       <- twinspan(X_input, modif = TRUE, clusters = 2, diss = dist)
  twi_lines <- as.character(unlist(tw$twi[[1]]))
  
  indicator_idx <- which(trimws(twi_lines) == "INDICATORS, together with their SIGN")
  if (length(indicator_idx) == 0) stop("Could not find TWINSPAN indicators block.")
  
  indicator_line  <- twi_lines[indicator_idx + 1]
  indicators_raw  <- str_trim(unlist(str_split(indicator_line, "\\s{2,}")))
  indicators_raw  <- indicators_raw[indicators_raw != ""]
  
  ind_parsed <- data.frame(species = character(), stringsAsFactors = FALSE)
  for (ind in indicators_raw) {
    if (!grepl("\\(", ind)) next
    m <- str_match(ind, "(.+?)(\\d+)\\(([+-])\\)")
    if (!is.na(m[1, 1]))
      ind_parsed <- rbind(ind_parsed,
                          data.frame(species = str_trim(m[1, 2]), stringsAsFactors = FALSE))
  }
  ind_abbrevnames <- gsub(" ", "", unique(ind_parsed$species))
  
  tw_classif <- tw$classif
  tw_group   <- merge(tw_classif, sample.info[, c("Station", hab_type)],
                      by.x = "plot.no", by.y = "Station", all.x = TRUE)
  tw_group$class <- as.character(tw_group$class)
  
  class_hab <- tw_group %>%
    group_by(class, .data[[hab_type]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(class) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(class, all_of(hab_type)) %>%
    as.data.frame()
  
  all_sp_class <- parse_species_class(twi_lines)
  all_sp_class$abbrev.name <- gsub(" ", "", all_sp_class$abbrev.name)
  all_sp_class$full.name   <- unname(name_map[all_sp_class$abbrev.name])
  all_sp_class <- merge(all_sp_class, class_hab, by.x = "twinspan_class", by.y = "class", all.x = TRUE)
  
  ind_fullnames <- unname(name_map[ind_abbrevnames])
  ind_fullnames <- ind_fullnames[!is.na(ind_fullnames)]
  
  # Restore full names before chi-sq and distances
  colnames(X_input) <- unname(name_map[colnames(X_input)])
  
  hab_vec <- sample.info[rownames(X_input), hab_type]
  
  chisq_stats <- lapply(colnames(X_input), function(sp) {
    tbl <- table(X_input[[sp]], hab_vec)
    ct  <- suppressWarnings(chisq.test(tbl))
    data.frame(full.name = sp,
               featureImportance = round(unname(ct$statistic), 4),
               pval = round(ct$p.value, 4),
               stringsAsFactors = FALSE)
  })
  chisq_df <- do.call(rbind, chisq_stats)
  
  prev_vec <- round((colSums(X_input > 0) / nrow(X_input)) * 100, 3)
  
  unmatched_prev <- all_sp_class$full.name[
    !is.na(all_sp_class$full.name) & !(all_sp_class$full.name %in% names(prev_vec))
  ]
  if (length(unmatched_prev) > 0)
    warning("prevalence NA for species not found in X_input columns: ",
            paste(unmatched_prev, collapse = ", "))
  
  indicator_df <- all_sp_class %>%
    rename(class = all_of(hab_type)) %>%
    left_join(chisq_df, by = "full.name") %>%
    mutate(
      IsIndSp    = as.integer(full.name %in% ind_fullnames),
      prevalence = prev_vec[match(full.name, names(prev_vec))],
      comparison = paste(ilevels, collapse = "_"),
      source     = "twinspan",
      data       = data_type
    ) %>%
    select(feature = full.name, featureImportance, class,
           IsIndSp, prevalence, comparison, source, data) %>%
    as.data.frame()
  
  adonis_results <- NULL
  ind_species    <- indicator_df$feature[indicator_df$IsIndSp == 1]
  ind_species    <- ind_species[!is.na(ind_species) & ind_species %in% colnames(X_input)]
  
  if (length(ind_species) == 1) {
    warning("Only 1 indicator species found for ", paste(ilevels, collapse = "_"),
            "; skipping adonis on indicator subset.")
  }
  if (length(ind_species) >= 2) {
    sub_X    <- X_input[, ind_species, drop = FALSE]
    sub_X    <- sub_X[rowSums(sub_X) > 0, , drop = FALSE]
    if (nrow(sub_X) < 2) {
      warning("sub_X has fewer than 2 rows after removing all-zero rows; skipping adonis for indicator species.")
    } else {
      sub_meta <- data.frame(sample = rownames(sub_X),
                             class  = sample.info[rownames(sub_X), hab_type])
      if (length(unique(sub_meta$class)) < 2) {
        warning("sub_X has only one class level after filtering; skipping adonis for indicator species.")
      } else {
        dist_sub    <- vegdist(sub_X, method = dist, binary = binary_flag)
        adonis_seed <- 100
        set.seed(adonis_seed)
        adonis_sub  <- adonis2(dist_sub ~ class, data = sub_meta)
        adonis_sub  <- data.frame(adonis_sub)["Model", , drop = FALSE]
        adonis_sub$comparison <- paste(ilevels, collapse = "_")
        adonis_sub$data       <- data_type
        adonis_sub$source     <- "twinspanSpecies"
        adonis_sub$features   <- length(ind_species)
        
        full_meta <- data.frame(sample = rownames(X_input),
                                class  = sample.info[rownames(X_input), hab_type])
        dist_all  <- vegdist(X_input, method = dist, binary = binary_flag)
        set.seed(adonis_seed)
        adonis_all <- adonis2(dist_all ~ class, data = full_meta)
        adonis_all <- data.frame(adonis_all)["Model", , drop = FALSE]
        adonis_all$comparison <- paste(ilevels, collapse = "_")
        adonis_all$data       <- data_type
        adonis_all$source     <- "allSpecies"
        adonis_all$features   <- ncol(X_input)
        
        adonis_results <- rbind(adonis_all, adonis_sub)
      }
    }
  }
  
  list(Comp_data    = Comp_data,
       twinspan_obj = tw,
       indicator.df = indicator_df,
       adonis       = adonis_results)
}

flatten_indicators <- function(twinspan_list, replica1, replica2) {
  sub <- lapply(twinspan_list, function(x) x[["indicator.df"]])
  sub <- Filter(Negate(is.null), sub)
  if (length(sub) == 0) return(NULL)
  sub <- lapply(sub, function(df) { rownames(df) <- NULL; df })
  out <- do.call(rbind, sub)
  # out$replica_train <- replica1
  # out$replica_test  <- replica2
  out
}

# ── Output directory ──────────────────────────────────────────────────────────
save_path_twinspan <- file.path(base_dir, "analyses", "Figures", "FigureS5", "results_replicas_analysis", "twinspan_output_data", "")

if (!dir.exists(save_path_twinspan)) {
  dir.create(save_path_twinspan, recursive = TRUE)
}

file_path_twinspan <- paste0(
  save_path_twinspan,
  "Twinspan_all_comparison_replicas_",
  paste(replica1, "vs", replica2, sep = "_"),
  "_prev_", prevalence_rate, ".Rda"
)

if (file.exists(file_path_twinspan)) {
  message("File already exists, skipping analysis: ", file_path_twinspan)
} else {
  
  twinspan.maxn <- list(); adonis.maxn <- list()
  twinspan.bin  <- list(); adonis.bin  <- list()
  
  # ── Pairwise habitat comparisons (BAY, LAG, BAR) ──────────────────────────
  for (i in seq_along(comp)) {
    ilevels    <- comp[[i]]
    key        <- paste(ilevels, collapse = "_")
    print(ilevels)
    
    ilevels_df <- acc[acc$hab %in% ilevels, ]
    
    # --- train on replica1 ---
    ilevels_df.train <- ilevels_df[ilevels_df$Transect == replica1, ]
    data.train <- prepare_data_for_analysis(ilevels_df.train, "hab", ilevels[[1]], ilevels[[2]])
    
    # --- test on replica2 ---
    ilevels_df.test  <- ilevels_df[ilevels_df$Transect == replica2, ]
    data.test  <- prepare_data_for_analysis(ilevels_df.test,  "hab", ilevels[[1]], ilevels[[2]])
    
    Comp_data <- list(
      X_train    = data.train$sample.X,
      y_train    = data.train$sample.class,
      info_train = data.train$sample.info,
      X_test     = data.test$sample.X,
      y_test     = data.test$sample.class,
      info_test  = data.test$sample.info
    )
    
    twinspan.maxn[[key]] <- run_twinspan_analysis(Comp_data, ilevels, hab_type = "hab", data_type = "maxN")
    if (!is.null(twinspan.maxn[[key]]$adonis)) adonis.maxn[[key]] <- twinspan.maxn[[key]]$adonis
    
    twinspan.bin[[key]]  <- run_twinspan_analysis(Comp_data, ilevels, hab_type = "hab", data_type = "pres/abs")
    if (!is.null(twinspan.bin[[key]]$adonis))  adonis.bin[[key]]  <- twinspan.bin[[key]]$adonis
  }
  
  # ── INSHORE vs OFFSHORE ───────────────────────────────────────────────────
  print(ilevels_hab)
  key_io <- paste(ilevels_hab, collapse = "_")
  
  # --- train on replica1 ---
  io_train   <- acc[acc$Transect == replica1, ]
  data_io_train <- prepare_data_for_analysis(io_train, "hab_inoff", ilevels_hab[[1]], ilevels_hab[[2]])
  
  # --- test on replica2 ---
  io_test    <- acc[acc$Transect == replica2, ]
  data_io_test  <- prepare_data_for_analysis(io_test,  "hab_inoff", ilevels_hab[[1]], ilevels_hab[[2]])
  
  Comp_data_io <- list(
    X_train    = data_io_train$sample.X,
    y_train    = data_io_train$sample.class,
    info_train = data_io_train$sample.info,
    X_test     = data_io_test$sample.X,
    y_test     = data_io_test$sample.class,
    info_test  = data_io_test$sample.info
  )
  
  twinspan.maxn[[key_io]] <- run_twinspan_analysis(Comp_data_io, ilevels_hab,
                                                   hab_type = "hab_inoff", data_type = "maxN")
  if (!is.null(twinspan.maxn[[key_io]]$adonis)) adonis.maxn[[key_io]] <- twinspan.maxn[[key_io]]$adonis
  
  twinspan.bin[[key_io]]  <- run_twinspan_analysis(Comp_data_io, ilevels_hab,
                                                   hab_type = "hab_inoff", data_type = "pres/abs")
  if (!is.null(twinspan.bin[[key_io]]$adonis))  adonis.bin[[key_io]]  <- twinspan.bin[[key_io]]$adonis
  
  # ── Flatten indicator tables ──────────────────────────────────────────────
  twinspan.maxn.sub <- flatten_indicators(twinspan.maxn, replica1, replica2)
  twinspan.bin.sub  <- flatten_indicators(twinspan.bin,  replica1, replica2)
  
  # ── Save ──────────────────────────────────────────────────────────────────
  save(
    adonis.bin, 
    adonis.maxn,
    twinspan.bin, 
    twinspan.maxn,
    twinspan.bin.sub, 
    twinspan.maxn.sub,
    file = file_path_twinspan
  )
  message("File saved: ", file_path_twinspan)
}