# =============================================================================
# Script name: codeTwinspan_all_analyses_split_data.R
# Author: Estephe Kana, Eugeni Belda, Edi Prifti
# Date created: 2025-01-01
# Purpose: Perform TWINSPAN with 80/20 stratified train/test split across multiple seeds.
# Inputs:  data/AllDataexportefixed.txt; args: prevalence_rate (%), maxseed
# Outputs: analyses/analysis_outputs/twinspan_output_data/analysis_80_20/
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script:
# cd /home/estephe/Documents/Papers/2026-aime_key_fish_bruvs/Bruvs-analysis-code/analyses/scripts
# Rscript codeTwinspan_all_analyses_split_data.R 10 10
# -----------------------------------------------------------------------------


if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("reshape2")) install.packages("reshape2"); library(reshape2)
if (!require("vegan")) install.packages("vegan"); library(vegan)
if (!require("stringr")) install.packages("stringr"); library(stringr)
if (!require("twinspanR")) install.packages("twinspanR"); library(twinspanR)
if (!require("caTools")) install.packages("caTools"); library(caTools)

args <- commandArgs(trailingOnly = TRUE)
prevalence_rate <- as.numeric(args[1])
maxseed <- as.numeric(args[2])

# Detect project root automatically (works for any user/machine after cloning)
if (!require("here", quietly = TRUE)) install.packages("here")
base_dir <- here::here()

dataset <- file.path(base_dir, "data", "AllDataexportefixed.txt")
mydata <- read.delim(dataset, header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE)

mydata$G_SP <- gsub("_", " ", mydata$G_SP)
mydata$Transect <- paste(substr(mydata$Transect, 3, 7))
mydata$Site <- factor(mydata$Site, levels = c("ABAR", "ALAG", "ABAI", "MBAR", "MLAG", "MBAI"))
# mydata$Site <- do.call(rbind, lapply(mydata$Site, gsub, pattern = "BAI", replacement = "BAY"))[, 1]
mydata$Site <- gsub("BAI", "BAY", mydata$Site)
# mydata$Station <- do.call(rbind, lapply(mydata$Station, gsub, pattern = "BAI", replacement = "BAY"))[, 1]
mydata$Station <- gsub("BAI", "BAY", mydata$Station)

matablepivot <- mydata[, c("Site", "Station", "G_SP", "MaxN")]
abund <- reshape2::dcast(matablepivot, Station + Site ~ G_SP, value.var = "MaxN")
abund[is.na(abund)] <- 0
rownames(abund) <- unique(abund$Station)
abund$hab <- substring(abund$Site, 2)
acc <- abund
acc$hab_inoff <- ifelse(acc$hab %in% c("BAY", "LAG"), "INSHORE", ifelse(acc$hab == "BAR", "OFFSHORE", NA))

comp <- combn(x = unique(acc$hab), m = 2, simplify = FALSE)
ilevels_hab <- list("INSHORE", "OFFSHORE")

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

prepare_data_for_analysis <- function(sample_data, habitat_type, hab1, hab2) {
  if (hab1 == hab2) stop("It seems that the habitats you entered are the same.")
  sample_data <- sample_data[!is.na(sample_data[[habitat_type]]) &
                               sample_data[[habitat_type]] %in% c(hab1, hab2), ]
  sample_data <- sample_data[, colSums(sample_data > 0) > 0]
  X <- sample_data[, -match(c("Station", "Site", "hab", "hab_inoff"), colnames(sample_data))]
  sample_X <- get_sample_by_prevalence(X, prevalence_rate)
  unprevalent_species <- setdiff(colnames(X), colnames(sample_X))
  sample_data <- sample_data[, !(colnames(sample_data) %in% unprevalent_species)]
  sample_data <- sample_data[rownames(sample_data) %in% rownames(sample_X), ]
  sample_data <- sample_data[rowSums(sample_data > 0) > 0, ]
  sample_data.class <- ifelse(sample_data[[habitat_type]] == hab1, -1, ifelse(sample_data[[habitat_type]] == hab2, 1, NA))
  sample_data.X <- sample_data[, -match(c("Station", "Site", "hab", "hab_inoff"), colnames(sample_data))]
  list(sample.info = sample_data, 
       sample.X = sample_data.X, 
       sample.class = sample_data.class)
}

split_train_test <- function(sample.info, sample.X, sample.class, seed,
                             hab_type = "hab", split_ratio = 0.8) {
  class_counts <- table(sample.info[[hab_type]])
  if (any(class_counts < 2)) {
    stop("Not enough samples per class to perform a stratified train/test split. ",
         "Counts: ", paste(names(class_counts), class_counts, sep = "=", collapse = ", "))
  }
  set.seed(seed)
  split_ind <- caTools::sample.split(Y = sample.info[[hab_type]], SplitRatio = split_ratio)

  # Verify both classes are represented in train and test
  train_classes <- unique(sample.info[[hab_type]][ split_ind])
  test_classes  <- unique(sample.info[[hab_type]][!split_ind])
  if (length(train_classes) < 2 || length(test_classes) < 2) {
    stop("Split resulted in only one class in train or test set. ",
         "Consider using a different seed or adjusting split_ratio.")
  }

  list(
    X_train    = sample.X[ split_ind,  , drop = FALSE],
    y_train    = sample.class[ split_ind],
    info_train = sample.info[ split_ind, , drop = FALSE],
    X_test     = sample.X[!split_ind,  , drop = FALSE],
    y_test     = sample.class[!split_ind],
    info_test  = sample.info[!split_ind, , drop = FALSE]
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
  
  # Relaxed: species name = one or more capitalised/lowercase words
  species_lines <- candidate_lines[
    str_detect(candidate_lines,
               "^\\s+\\d+\\s+[A-Z][A-Za-z.]+(?:\\s+[a-z.]+)+\\s+\\S+\\s+[01]\\s*$")
  ]
  if (length(species_lines) == 0) stop("No species lines matched in the tabulation block.")
  
  parsed <- lapply(species_lines, function(line) {
    sp_class  <- str_trim(str_extract(line, "[01]\\s*$"))
    # Capture everything between the leading number and the trailing abbreviation+class
    sp_abbrev <- str_match(line,
                           "^\\s*\\d+\\s+((?:[A-Z][A-Za-z.]+)(?:\\s+[a-z.]+)+)\\s+\\S+\\s+[01]\\s*$")[1, 2]
    data.frame(abbrev.name = sp_abbrev, twinspan_class = paste0("*", sp_class),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, parsed)
}

run_twinspan_analysis <- function(Comp_data, ilevels, hab_type, data_type = "maxN") {
  X_train <- Comp_data$X_train
  y_train <- Comp_data$y_train
  sample.info <- Comp_data$info_train

  if (data_type == "pres/abs") {
    X_input <- as.data.frame(apply(X_train, 2, function(x) ifelse(x == 0, 0, 1)))
    dist <- "jaccard"
    binary_flag <- TRUE
  } else {
    X_input <- X_train
    dist <- "bray"
    binary_flag <- FALSE
  }
  
  # full_names  <- colnames(X_input)
  # abbrev_keys <- make.unique(substr(gsub("[^A-Za-z0-9]", "", full_names), 1, 8), sep = ".")
  # name_map    <- setNames(full_names, abbrev_keys)   # abbrev -> full
  # colnames(X_input) <- abbrev_keys
  
  full_names <- colnames(X_input)
  
  # Split genus and species
  split_names <- strsplit(full_names, " ")
  
  genus   <- sapply(split_names, `[`, 1)
  species <- sapply(split_names, `[`, 2)
  
  # Create abbreviation: first 4 letters genus + first 4 letters species
  abbrev_keys <- paste0(
    substr(genus, 1, 4),
    substr(species, 1, 4)
  )
  
  # Ensure uniqueness
  abbrev_keys <- make.unique(abbrev_keys, sep = ".")
  
  # Mapping table
  name_map <- setNames(full_names, abbrev_keys)
  
  # Replace colnames
  colnames(X_input) <- abbrev_keys
  
  tw <- twinspan(X_input, modif = TRUE, clusters = 2, diss = dist)
  twi_lines <- as.character(unlist(tw$twi[[1]]))
  
  indicator_idx <- which(trimws(twi_lines) == "INDICATORS, together with their SIGN")
  if (length(indicator_idx) == 0) stop("Could not find TWINSPAN indicators block.")
  
  indicator_line <- twi_lines[indicator_idx + 1]
  indicators_raw <- str_trim(unlist(str_split(indicator_line, "\\s{2,}")))
  indicators_raw <- indicators_raw[indicators_raw != ""]

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
  tw_group <- merge(tw_classif, sample.info[, c("Station", hab_type)], by.x = "plot.no", by.y = "Station", all.x = TRUE)
  tw_group$class <- as.character(tw_group$class)
  
  class_hab <- tw_group %>%
    group_by(class, .data[[hab_type]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(class) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(class, all_of(hab_type)) %>%
    as.data.frame()

  # all_sp_class <- parse_species_class(twi_lines)
  # all_sp_class$abbrev.name <- gsub(" ", "", all_sp_class$abbrev.name)
  # all_sp_class <- merge(all_sp_class, tw$spnames, by = "abbrev.name", all.x = TRUE)
  # all_sp_class <- merge(all_sp_class, class_hab, by.x = "twinspan_class", by.y = "class", all.x = TRUE)
  
  all_sp_class <- parse_species_class(twi_lines)
  all_sp_class$abbrev.name <- gsub(" ", "", all_sp_class$abbrev.name)
  # Restore full names from our map rather than relying on tw$spnames
  all_sp_class$full.name <- unname(name_map[all_sp_class$abbrev.name])
  # (skip the tw$spnames merge since we resolve names directly)
  all_sp_class <- merge(all_sp_class, class_hab, by.x = "twinspan_class", by.y = "class", all.x = TRUE)
  
  # ind_fullnames <- all_sp_class[all_sp_class$abbrev.name %in% ind_abbrevnames, "full.name"]
  ind_fullnames <- unname(name_map[ind_abbrevnames])
  ind_fullnames <- ind_fullnames[!is.na(ind_fullnames)]

  hab_vec <- sample.info[rownames(X_input), hab_type]
  # chisq_stats <- lapply(colnames(X_input), function(sp) {
  #   tbl <- table(X_input[[sp]], hab_vec)
  #   ct <- suppressWarnings(chisq.test(tbl))
  #   data.frame(full.name = sp, 
  #              featureImportance = round(unname(ct$statistic), 4), 
  #              pval = round(ct$p.value, 4), 
  #              stringsAsFactors = FALSE)
  # })
  # chisq_df <- do.call(rbind, chisq_stats)
  
  # Restore full species names on X_input before distance calculations
  colnames(X_input) <- unname(name_map[colnames(X_input)])
  
  chisq_stats <- lapply(colnames(X_input), function(sp) {
    tbl <- table(X_input[[sp]], hab_vec)
    ct <- suppressWarnings(chisq.test(tbl))
    data.frame(full.name = sp,   
               featureImportance = round(unname(ct$statistic), 4), 
               pval = round(ct$p.value, 4), 
               stringsAsFactors = FALSE)
  })
  chisq_df <- do.call(rbind, chisq_stats)
  
  prev_vec <- round((colSums(X_input > 0) / nrow(X_input)) * 100, 3)
  # names(prev_vec) <- unname(name_map[names(prev_vec)])
  
  # prev_vec <- round((colSums(X_input > 0) / nrow(X_input)) * 100, 3)
  
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
      IsIndSp = as.integer(full.name %in% ind_fullnames),
      prevalence = prev_vec[match(full.name, names(prev_vec))],
      comparison = paste(ilevels, collapse = "_"),
      source = "twinspan",
      data = data_type
    ) %>%
    select(feature = full.name, featureImportance, pval, class, 
           IsIndSp, prevalence, comparison, source, data) %>%
    as.data.frame()
  
  adonis_results <- NULL
  ind_species <- indicator_df$feature[indicator_df$IsIndSp == 1]
  ind_species <- ind_species[!is.na(ind_species) & ind_species %in% colnames(X_input)]
  if (length(ind_species) == 1) {
    warning("Only 1 indicator species found for ", paste(ilevels, collapse="_"), 
            "; skipping adonis on indicator subset.")
  }
  if (length(ind_species) >= 2) {
    sub_X <- X_input[, ind_species, drop = FALSE]
    sub_X <- sub_X[rowSums(sub_X) > 0, , drop = FALSE]
    if (nrow(sub_X) < 2) {
      warning("sub_X has fewer than 2 rows after removing all-zero rows; skipping adonis for indicator species.")
    }
    sub_meta  <- data.frame(sample = rownames(sub_X),
                            class  = sample.info[rownames(sub_X), hab_type])
    if (length(unique(sub_meta$class)) < 2) {
      warning("sub_X has only one class level after filtering all-zero rows; skipping adonis for indicator species.")
      adonis_results <- NULL
    }
    dist_sub  <- vegdist(sub_X, method = dist, binary = binary_flag)
    
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
    
    if (length(unique(full_meta$class)) < 2) {
      warning("X_input has only one class level; skipping adonis for all species.")
    } 
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

  list(Comp_data = Comp_data,
       twinspan_obj = tw, 
       indicator.df = indicator_df,
       adonis = adonis_results)
}

flatten_indicators <- function(twinspan_list, seed) {
  sub <- lapply(twinspan_list, function(x) x[["indicator.df"]])
  sub <- Filter(Negate(is.null), sub)
  if (length(sub) == 0) return(NULL)
  sub <- lapply(sub, function(df) {
    rownames(df) <- NULL
    df
  })
  out <- do.call(rbind, sub)
  out$seed <- seed
  out
}

save_path_twinspan <- file.path(base_dir, "analyses", "analysis_outputs", "twinspan_output_data", "analysis_80_20")

for (seed in 1:maxseed) {
  file_path_twinspan <- paste0(save_path_twinspan, "/Twinspan_all_analyses_split_80_20_data_prev_",
                               prevalence_rate, "_seed_", seed, ".Rda")
  if (file.exists(file_path_twinspan)) {
    message("File already exists, skipping analysis: ", file_path_twinspan)
    next
  }
  
  twinspan.maxn <- list(); adonis.maxn <- list()
  twinspan.bin  <- list(); adonis.bin  <- list()
  
  for (i in seq_along(comp)) {
    ilevels    <- comp[[i]]
    key        <- paste(ilevels, collapse = "_")
    ilevels_df <- acc[acc$hab %in% ilevels, ]
    
    data_prep  <- prepare_data_for_analysis(ilevels_df, "hab", ilevels[[1]], ilevels[[2]])
    
    split_data <- split_train_test(data_prep$sample.info, data_prep$sample.X,
                                   data_prep$sample.class, seed = seed, hab_type = "hab")
    Comp_data <- list(
      X_train    = split_data$X_train,  y_train    = split_data$y_train,
      X_test     = split_data$X_test,   y_test     = split_data$y_test,
      info_train = split_data$info_train, info_test = split_data$info_test
    )
    
    twinspan.maxn[[key]] <- run_twinspan_analysis(Comp_data, ilevels, hab_type = "hab", data_type = "maxN")
    if (!is.null(twinspan.maxn[[key]]$adonis)) adonis.maxn[[key]] <- twinspan.maxn[[key]]$adonis
    
    twinspan.bin[[key]]  <- run_twinspan_analysis(Comp_data, ilevels, hab_type = "hab", data_type = "pres/abs")
    if (!is.null(twinspan.bin[[key]]$adonis))  adonis.bin[[key]]  <- twinspan.bin[[key]]$adonis
  }
  
  data_io  <- prepare_data_for_analysis(acc, "hab_inoff", ilevels_hab[[1]], ilevels_hab[[2]])
  
  split_io <- split_train_test(data_io$sample.info, data_io$sample.X,
                               data_io$sample.class, seed = seed, hab_type = "hab_inoff")
  
  Comp_data_io <- list(
    X_train    = split_io$X_train,    y_train    = split_io$y_train,
    X_test     = split_io$X_test,     y_test     = split_io$y_test,
    info_train = split_io$info_train, info_test  = split_io$info_test
  )
  
  key_io <- paste(ilevels_hab, collapse = "_")
  
  twinspan.maxn[[key_io]] <- run_twinspan_analysis(Comp_data_io, ilevels_hab,
                                                   hab_type = "hab_inoff", data_type = "maxN")
  if (!is.null(twinspan.maxn[[key_io]]$adonis)) adonis.maxn[[key_io]] <- twinspan.maxn[[key_io]]$adonis
  
  twinspan.bin[[key_io]]  <- run_twinspan_analysis(Comp_data_io, ilevels_hab,
                                                   hab_type = "hab_inoff", data_type = "pres/abs")
  if (!is.null(twinspan.bin[[key_io]]$adonis))  adonis.bin[[key_io]]  <- twinspan.bin[[key_io]]$adonis
  
  twinspan.maxn.sub <- flatten_indicators(twinspan.maxn, seed)
  twinspan.bin.sub  <- flatten_indicators(twinspan.bin,  seed)
  
  save(
    adonis.bin, adonis.maxn,
    twinspan.bin, twinspan.maxn,
    twinspan.bin.sub, twinspan.maxn.sub,
    file = file_path_twinspan
  )
  message("File saved: ", file_path_twinspan)
  
}
