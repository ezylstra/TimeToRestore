################################################################################
# Downloading iNaturalist data

# Erin Zylstra
# 2024-05-01
################################################################################

library(rgbif)
library(dplyr)

# Going to try downloading iNaturalist data through the rgbif package, which 
# seems to be better supported and updated than the rinat package. 

# Load csv with priority species
spp <- read.csv("data/priority-species.csv") %>%
  mutate(sci_names = paste0(genus, " ", species))

# Remove species from list if they don't have species_id in NPN database
spp <- spp %>%
  filter(!is.na(species_id))

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
  pred("datasetKey", "50c9509d-22c7-4a22-a47d-8c48425ef4a7"),
  format = "SIMPLE_CSV"
)

# This prints some stuff, including a Download key.
dk <- "0010217-240425142415019"

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

# Remove the original csv file and write this smaller one to file
file.remove(inat_filename)
write.csv(inat, inat_filename, row.names = FALSE)

# Note that we haven't filtered by year. There are <=10 records/year for all
# species combined prior to 2000.

# How many records for each species in the 4 target states?
states4 <- c("Louisiana", "New Mexico", "Oklahoma", "Texas")
inat <- inat %>%
  mutate(states4 = ifelse(state %in% states4, 1, 0))
# Counts in CASC states, combined
count(filter(inat, states4 == 1), species)
# Counts for each CASC state
table(inat$species[inat$states4 == 1], inat$state[inat$states4 == 1])
