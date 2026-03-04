# =============================================================================
# Script name: generate_BRUVS_data_object.R
# Author: Estephe Kana, Eugeni Belda, Edi Prifti
# Date created: 2026-03-01
# Purpose: Load the raw BRUVS species abundance Excel file and produce two
#          processed data objects (abundance and presence/absence) used by
#          all downstream analysis and figure scripts.
# Inputs:  data/BRUVS_Sandy_Area_Noumea_2016_Vigliola_Zenodo.xlsx
# Outputs: data/video_data_object.rda
#          data/video_data_object_pres_abs.rda
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script:
# cd /path/to/Bruvs-analysis-code/analyses/scripts
# Rscript generate_BRUVS_data_object.R
# -----------------------------------------------------------------------------

# Detect project root automatically (works for any user/machine after cloning)
if (!require("here", quietly = TRUE)) install.packages("here")
base_dir <- here::here()

# -----------------------------
# Load or install required packages
# -----------------------------

packages <- c(
  "dplyr", "reshape2", "readxl", "vegan"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# -----------------------------
# Read data from Excel file
# -----------------------------

input_file <- file.path(base_dir, "data", "BRUVS_Sandy_Area_Noumea_2016_Vigliola_Zenodo.xlsx")

video_data <- readxl::read_excel(input_file, sheet = 1)

head(video_data)

# -----------------------------------------
# Store data in structured object
# -----------------------------------------

sm <- list()
sm_presAbs <- list()

abund <- video_data[, !(colnames(video_data) %in% c("Station", "Date", "Longitude", "Latitude", "Site", "Habitat", "Zone"))]

rownames(abund) <- video_data$Station

# get the community matrix for abundance data
sm$X <- as.matrix(abund)

# get the community matrix for presence/absence data
presAbs <- decostand(abund, method =  "pa")
sm_presAbs$X <- as.matrix(presAbs)


# get sample info for abundance data
rownames(video_data) <- video_data$Station
sm$sample_info <- video_data

# get sample info for abundance data
sm_presAbs$sample_info <- cbind(video_data[, (colnames(video_data) %in% c("Station", "Date", "Longitude", "Latitude", "Site", "Habitat", "Zone"))], presAbs)

# -----------------------------
# Save data objects
# -----------------------------
output_dir <- file.path(base_dir, "data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

abund_output_file   <- file.path(output_dir, "video_data_object.rda")
presAbs_output_file <- file.path(output_dir, "video_data_object_pres_abs.rda")

save(sm, file = abund_output_file)
sm <- sm_presAbs
save(sm, file = presAbs_output_file)
