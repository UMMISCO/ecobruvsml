# =============================================================================
# Script name: convert_predomics_input_to_xlsx.R
# Author: Estephe Kana, Eugeni Belda, Edi Prifti
# Date created: 2024-06-01
# Purpose: Convert BRUVS data (used in Predomics analysis) from R object format to Excel (.xlsx) for sharing.
# Inputs:  analyses/analysis_outputs/video_data_object.rda
# Outputs: data/environmental_data_lagoon_2016_60_samples_predomics.xlsx
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script:
# cd /home/estephe/Documents/Papers/2026-aime_key_fish_bruvs/Bruvs-analysis-code/analyses/scripts
# Rscript convert_predomics_input_to_xlsx.R
# -----------------------------------------------------------------------------

# Convert the BRUVS data used in Predomics in a xlsx file

install.packages("openxlsx")
library(openxlsx)

# load BRUVS data used in predomics
# Detect project root automatically (works for any user/machine after cloning)
if (!require("here", quietly = TRUE)) install.packages("here")
base_dir <- here::here()

load(file=file.path(base_dir, "data", "video_data_object.rda"))

# Create and populate workbook
wb <- createWorkbook()
for (i in seq_along(sm)) {
  sheetname <- names(sm)[i]
  if (is.null(sheetname) || sheetname == "") {
    sheetname <- paste0("Sheet", i)
  }
  addWorksheet(wb, sheetname)
  writeData(wb, sheetname, sm[[i]])
}
saveWorkbook(wb, file.path(base_dir, "data", "environmental_data_lagoon_2016_60_samples_predomics.xlsx"), overwrite = TRUE)
