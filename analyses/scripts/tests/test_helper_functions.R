# =============================================================================
# Script name: test_helper_functions.R
# Author: Estephe Kana, Eugeni Belda, Edi Prifti
# Date created: 2026-03-01
# Purpose: Unit tests for shared helper functions used across BRUVS analysis
#          scripts: get_sample_by_prevalence() and prepare_data_for_analysis().
# Inputs:  (none - uses synthetic test data)
# Outputs: Test results printed to console
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the tests:
# cd /home/estephe/Documents/Papers/2026-aime_key_fish_bruvs/Bruvs-analysis-code/analyses/scripts
# Rscript tests/test_helper_functions.R
# Or with testthat: testthat::test_file("tests/test_helper_functions.R")
# -----------------------------------------------------------------------------

if (!require("testthat")) install.packages("testthat"); library(testthat)

# =============================================================================
# Helper functions under test (copied from the analysis scripts)
# =============================================================================

get_sample_by_prevalence <- function(abundance_matrix, prevalence_rate) {
  message(paste("You set the prevalence rate to:", paste0(prevalence_rate, "%")))
  prevalence_rate <- prevalence_rate / 100
  if (prevalence_rate == 0) message("Using the overall data without any filtering.")
  if (prevalence_rate == 1) message("Using species that appear in all samples.")
  if (prevalence_rate < 0 || prevalence_rate > 1) stop("Prevalence rate must be between 0 and 1.")
  num_samples     <- nrow(abundance_matrix)
  min_occurrences <- ceiling(prevalence_rate * num_samples)
  valid_features  <- colSums(abundance_matrix > 0) >= min_occurrences
  sampled_matrix  <- abundance_matrix[, valid_features, drop = FALSE]
  sampled_matrix  <- sampled_matrix[rowSums(sampled_matrix > 0) > 0, ]
  return(sampled_matrix)
}

prepare_data_for_analysis <- function(sample_data, habitat_type, hab1, hab2,
                                      prevalence_rate = 0) {
  if (hab1 == hab2) stop("It seems that the habitats you entered are the same.")
  sample_data <- sample_data[, colSums(sample_data > 0) > 0]
  X           <- sample_data[, -match(c("Station", "Site", "hab", "hab_inoff"),
                                       colnames(sample_data))]
  sample_X    <- get_sample_by_prevalence(X, prevalence_rate)
  unprevalent <- setdiff(colnames(X), colnames(sample_X))
  sample_data <- sample_data[, !(colnames(sample_data) %in% unprevalent)]
  sample_data <- sample_data[rownames(sample_data) %in% rownames(sample_X), ]
  sample_data <- sample_data[rowSums(sample_data > 0) > 0, ]
  sample_data.class <- ifelse(sample_data[[habitat_type]] == hab1, -1,
                               ifelse(sample_data[[habitat_type]] == hab2, 1, NA))
  sample_data.X <- sample_data[, -match(c("Station", "Site", "hab", "hab_inoff"),
                                         colnames(sample_data))]
  return(list(sample.info  = sample_data,
              sample.X     = sample_data.X,
              sample.class = sample_data.class))
}

# =============================================================================
# Shared synthetic data
# =============================================================================

make_abund <- function() {
  # 9 samples x 5 species; 3 habitats (3 samples each)
  set.seed(42)
  mat <- matrix(
    c(5, 3, 0, 0, 1,   # BAR
      4, 2, 0, 0, 0,
      6, 1, 0, 0, 2,
      0, 0, 7, 5, 0,   # BAY
      0, 0, 6, 4, 0,
      0, 0, 8, 3, 1,
      0, 1, 0, 6, 4,   # LAG
      0, 0, 0, 5, 3,
      0, 2, 0, 7, 2),
    nrow = 9, ncol = 5, byrow = TRUE,
    dimnames = list(
      paste0("S", 1:9),
      c("SpA", "SpB", "SpC", "SpD", "SpE")
    )
  )
  df <- as.data.frame(mat)
  df$Station  <- paste0("S", 1:9)
  df$Site     <- c(rep("ABAR", 3), rep("ABAY", 3), rep("ALAG", 3))
  df$hab      <- c(rep("BAR", 3), rep("BAY", 3), rep("LAG", 3))
  df$hab_inoff <- ifelse(df$hab == "BAR", "OFFSHORE", "INSHORE")
  df
}

# =============================================================================
# Tests: get_sample_by_prevalence
# =============================================================================

test_that("prevalence 0% keeps all species and samples with any non-zero value", {
  mat <- make_abund()[, c("SpA","SpB","SpC","SpD","SpE")]
  result <- get_sample_by_prevalence(mat, 0)
  # All 5 species present; all 9 samples have at least one non-zero
  expect_equal(ncol(result), 5)
  expect_equal(nrow(result), 9)
})

test_that("prevalence 100% keeps only species present in ALL samples", {
  mat <- make_abund()[, c("SpA","SpB","SpC","SpD","SpE")]
  result <- get_sample_by_prevalence(mat, 100)
  # SpD is 0 in row 1-3 (BAR), SpA is 0 in row 4-9, etc. No species is
  # present in all 9 samples -> result should have 0 columns
  expect_equal(ncol(result), 0)
})

test_that("prevalence filters species below threshold correctly", {
  mat <- make_abund()[, c("SpA","SpB","SpC","SpD","SpE")]
  # SpA: present in 3/9 samples (rows 1-3) = 33%
  # SpC: present in 3/9 samples (rows 4-6) = 33%
  # SpD: present in 6/9 samples (rows 4-9) = 67%
  result <- get_sample_by_prevalence(mat, 50)  # ≥ ceiling(0.5*9)=5 samples
  expect_true("SpD" %in% colnames(result))
  expect_false("SpA" %in% colnames(result))
  expect_false("SpC" %in% colnames(result))
})

test_that("prevalence removes samples that become empty after species filtering", {
  # Create matrix where filtering removes all species for some samples
  mat <- data.frame(
    SpRare = c(1, 0, 0, 0, 0),   # only in 1/5 samples (20%)
    SpComm = c(1, 1, 1, 1, 0),   # in 4/5 samples (80%)
    row.names = paste0("S", 1:5)
  )
  result <- get_sample_by_prevalence(mat, 50)  # ceil(0.5*5)=3; SpComm passes
  # S5 has SpRare=0 and SpComm=0 -> removed
  expect_false("S5" %in% rownames(result))
  expect_true("SpComm" %in% colnames(result))
})

test_that("invalid prevalence rate throws error", {
  mat <- make_abund()[, c("SpA","SpB")]
  expect_error(get_sample_by_prevalence(mat, -10), "Prevalence rate must be between 0 and 1")
  expect_error(get_sample_by_prevalence(mat, 110), "Prevalence rate must be between 0 and 1")
})

test_that("result is a data.frame with numeric values only", {
  mat <- make_abund()[, c("SpA","SpB","SpC","SpD","SpE")]
  result <- get_sample_by_prevalence(mat, 10)
  expect_s3_class(result, "data.frame")
  expect_true(all(sapply(result, is.numeric)))
})

# =============================================================================
# Tests: prepare_data_for_analysis
# =============================================================================

test_that("returns list with sample.info, sample.X, sample.class", {
  df  <- make_abund()
  sub <- df[df$hab %in% c("BAR", "BAY"), ]
  out <- prepare_data_for_analysis(sub, "hab", "BAR", "BAY", prevalence_rate = 0)
  expect_named(out, c("sample.info", "sample.X", "sample.class"))
})

test_that("sample.class is -1 for hab1 and 1 for hab2", {
  df  <- make_abund()
  sub <- df[df$hab %in% c("BAR", "BAY"), ]
  out <- prepare_data_for_analysis(sub, "hab", "BAR", "BAY", prevalence_rate = 0)
  expect_true(all(out$sample.class[out$sample.info$hab == "BAR"] == -1))
  expect_true(all(out$sample.class[out$sample.info$hab == "BAY"] ==  1))
})

test_that("sample.X contains no metadata columns", {
  df  <- make_abund()
  sub <- df[df$hab %in% c("BAR", "LAG"), ]
  out <- prepare_data_for_analysis(sub, "hab", "BAR", "LAG", prevalence_rate = 0)
  meta_cols <- c("Station", "Site", "hab", "hab_inoff")
  expect_false(any(meta_cols %in% colnames(out$sample.X)))
})

test_that("error when hab1 equals hab2", {
  df  <- make_abund()
  sub <- df[df$hab == "BAR", ]
  expect_error(
    prepare_data_for_analysis(sub, "hab", "BAR", "BAR"),
    "same"
  )
})

test_that("prevalence filtering works inside prepare_data_for_analysis", {
  df  <- make_abund()
  sub <- df[df$hab %in% c("BAY", "LAG"), ]
  # SpA is absent in all BAY+LAG rows -> already filtered by colSums > 0
  out0  <- prepare_data_for_analysis(sub, "hab", "BAY", "LAG", prevalence_rate = 0)
  out50 <- prepare_data_for_analysis(sub, "hab", "BAY", "LAG", prevalence_rate = 50)
  # With 50% prevalence, fewer species retained
  expect_true(ncol(out50$sample.X) <= ncol(out0$sample.X))
})

test_that("INSHORE_OFFSHORE comparison with hab_inoff column works", {
  df  <- make_abund()
  out <- prepare_data_for_analysis(df, "hab_inoff", "INSHORE", "OFFSHORE",
                                   prevalence_rate = 0)
  expect_true(all(out$sample.class[out$sample.info$hab_inoff == "INSHORE"]  == -1))
  expect_true(all(out$sample.class[out$sample.info$hab_inoff == "OFFSHORE"] ==  1))
})

test_that("nrow(sample.X) == length(sample.class)", {
  df  <- make_abund()
  sub <- df[df$hab %in% c("BAR", "BAY"), ]
  out <- prepare_data_for_analysis(sub, "hab", "BAR", "BAY", prevalence_rate = 0)
  expect_equal(nrow(out$sample.X), length(out$sample.class))
})

# =============================================================================
# Run all tests
# =============================================================================

cat("\n=== Running unit tests for BRUVS helper functions ===\n\n")
# testthat collects results from all test_that() blocks above automatically.
# Run this file directly with:  Rscript tests/test_helper_functions.R
# Or run interactively with:    testthat::test_file("tests/test_helper_functions.R")
