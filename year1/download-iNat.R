################################################################################
# Downloading iNaturalist data

# Erin Zylstra
# 2024-07-09
################################################################################

library(rgbif)
library(dplyr)
library(stringr)

# Going to try downloading iNaturalist data through the rgbif package, which 
# seems to be better supported and updated than the rinat package. 

# Load csv with priority species
spp <- read.csv("data/priority-species.csv") %>%
  mutate(sci_names = paste0(genus, " ", species))

# Remove species from list if they don't have species_id in NPN database
spp <- spp %>%
  filter(!is.na(species_id))

# Identify years of interest (using same period we used for NPN data)
yrs <- 2013:2023

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

# To access research grade observations in iNaturalist use:
# datasetKey = "50c9509d-22c7-4a22-a47d-8c48425ef4a7"

occ_download(
  pred_default(),
  pred_in("speciesKey", spp$taxonkey),
  pred("country", "US"),
  pred_in("year", yrs),
  pred("datasetKey", "50c9509d-22c7-4a22-a47d-8c48425ef4a7"),
  format = "SIMPLE_CSV"
)

# This prints some stuff, including a Download key.
dk <- "0018473-240626123714530"

# Check the status of download:
occ_download_wait(dk)
# Retrieve the data:
occ_download_get(dk, path = "data/iNat")

# Unzip and rename datafile
nspp <- nrow(spp)
inat_filename <- paste0("data/iNat/inaturalist-", nspp, "spp.csv")
unzip(zipfile = paste0("data/iNat/", dk, ".zip"),
      exdir = "data/iNat/")
file.rename(from = paste0("data/iNat/", dk, ".csv"), 
            to = inat_filename)

# Load data (note that csv is tab delimited)
inat <- read.csv(inat_filename, sep = "\t", header = TRUE, 
                 na.strings = c(NA, "")) %>%
  select(gbifID, species, stateProvince, occurrenceStatus, decimalLatitude, 
         decimalLongitude, coordinateUncertaintyInMeters,day, month, 
         year, taxonKey, basisOfRecord, collectionCode) %>%
  rename(state = stateProvince,
         lat = decimalLatitude,
         lon = decimalLongitude,
         coord_uncert = coordinateUncertaintyInMeters)

# If all basisOfRecord == HUMAN_OBSERVATION, all collectionCode == Observations,
# and all occurrenceStatus == PRESENT, then we can delete these fields
if (all(inat$basisOfRecord == "HUMAN_OBSERVATION")) {
  inat <- select(inat, -basisOfRecord)
} else {
  warning("Not all basisOfRecord == HUMAN_OBSERVATION")
}
if (all(inat$collectionCode == "Observations")) {
  inat <- select(inat, -collectionCode)
} else {
  warning("Not all collectionCode == Observations")
}
if (all(inat$occurrenceStatus == "PRESENT")) {
  inat <- select(inat, -occurrenceStatus)
} else {
  warning("Not all occurrenceStatus == PRESENCE")
}

# Add in common name
inat <- inat %>%
  left_join(select(spp, c(sci_names, common_name)), 
            by = c("species" = "sci_names"))

# Leave code below commented out until we're sure we want to create new iNat
# dataset

# Save csv and put in a zip file. Remove csv and original zip file (it's huge)
  # write.csv(inat, inat_filename, row.names = FALSE)
  # zip_filename <- str_replace(inat_filename, ".csv", ".zip")
  # zip(zip_filename, files = inat_filename)

# Remove csv
  # file.remove(inat_filename)

# How many records for each species in the 4 target states?
states4 <- c("Louisiana", "New Mexico", "Oklahoma", "Texas")
inat <- inat %>%
  mutate(states4 = ifelse(state %in% states4, 1, 0))
# Counts in CASC states, combined
count(filter(inat, states4 == 1), common_name)
# Counts for each CASC state
table(inat$common_name[inat$states4 == 1], inat$state[inat$states4 == 1])
