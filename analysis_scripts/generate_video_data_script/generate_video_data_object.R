# ============================================================
# Script: Prepare video-based metabarcoding data
# Authors: Estephe Kana, Eugeni Belda & Edi Prifti
# Date: Sys.Date()
# Description:
#   - Loads and preprocesses species abundance data
#   - Creates and saves processed data object
# Outputs: video_data_object.rda, video_data_object_pres_abs.rda
# ============================================================

# -------------------------------
# How to run this script
# -------------------------------
# 1. Open a terminal.
# 2. Navigate to the script directory:
#     cd /ecobruvsml/analysis_scripts/generate_video_data_script
# 3. Run the script using Rscript:
#     Rscript generate_video_data_object.R
# -------------------------------

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
# ------------------------------

input_file <- "/data/BRUVS_Sandy_Area_Noumea_2016_Vigliola_Zenodo.xlsx"

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
abund_output_file <- "/data/video_data_object.rda"
presAbs_output_file <- "/data/video_data_object_pres_abs.rda"

save(sm, file = abund_output_file)
sm <- sm_presAbs
save(sm, file = presAbs_output_file)
