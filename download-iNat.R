################################################################################
# Downloading iNaturalist data

# Erin Zylstra
# 2024-04-30
################################################################################

library(rgbif)
library(dplyr)

rm(list = ls())

# Going to try downloading iNaturalist data through the rgbif package, which 
# seems to be better supported and updated than the rinat package. 

# Load csv with priority species
spp <- read.csv("data/priority-species.csv") %>%
  mutate(sci_names = paste0(genus, " ", species))

# Get list of scientific names
name_list <- c(spp$sci_names)
# Check that we can find matches in GBIF for each and get GBIF taxonkeys
gbif_info <- name_backbone_checklist(name_list) %>% data.frame()
if (any(gbif_info$matchType != "EXACT")) {
  prob_indices <- which(gbif_info$matchType != "EXACT")
  warning(paste0("The following scientific names do not have an exact match in GBIF:\n",
          paste(spp$sci_names[prob_indices], collapse = " | ")))
}
spp$taxonkey <- gbif_info$usageKey

occ_download(
  pred_default(),
  pred_in("speciesKey", spp$taxonkey[1:2]),
  pred("datasetKey", "50c9509d-22c7-4a22-a47d-8c48425ef4a7"),
  format = "SIMPLE_CSV"
)

# This prints some stuff, including a Download key.
dk <- "0009318-240425142415019"

# SIMPLE_CSV doesn't seem to be resulting in a comma delimited file.
### PICK UP HERE

# Check the status of download:
occ_download_wait(dk)
# Retrieve the data:
occ_download_get(dk, path = "data/iNat")
# Unzip and load data
unzip(zipfile = paste0("data/iNat/", dk, ".zip"),
      exdir = "data/iNat/")
file.rename(from = paste0("data/iNat/", dk, ".csv"), 
            to = "data/iNat/twospp_test.csv")

# Load data
inat <- read.csv("data/iNat/twospp_test.csv")

      


