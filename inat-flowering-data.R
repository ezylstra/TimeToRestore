library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
# library(shiny)
# library(bslib)
# library(sf)
library(ggplot2)
# library(leaflet)
# library(RColorBrewer)
# library(sortable)
# library(grid)
# library(cowplot)
library(terra)
library(tidyterra)

# Load/format observation data ------------------------------------------------#

df <- read.csv("data/status-intensity-flowers-May2026.csv")
df <- df %>%
  mutate(obsdate = ymd(obsdate)) %>%
  mutate(region = case_when(
    region == "Dallas" ~ "Dallas/Ft. Worth",
    region == "Austin" ~ "Austin/San Antonio",
    .default = region
  ))
# For now limit to 4 state area and 2025-2026
df <- df %>%
  filter(state4 == 1) %>%
  filter(yr >= 2025)

# If there's more than one observation of a plant on a given date, select one 
# with positive status and (higher) intensity value
df <- df %>%
  arrange(id, obsdate, php_id, desc(status), desc(midpoint)) %>%
  distinct(id, obsdate, php_id, .keep_all = TRUE)

# Put all data for a plant, date in the same row (wide form)
dfw <- df %>%
  select(-c(observation_id, php)) %>%
  pivot_wider(names_from = c(php_id),
              values_from = c(status, intensity_value, midpoint)) %>%
  data.frame()

# If flower = 0, then open flower must be 0. Remove observations with open
# flower = 1 and change any open flower = NA to 0.
dfw <- dfw %>%
  filter(!(!is.na(status_500) & !is.na(status_501) & status_500 == 0 & status_501 == 1)) %>%
  mutate(status_501 = ifelse(!is.na(status_500) & status_500 == 0, 0, status_501))

# If open flower = 1, then flower must be 1. Change any that are NA.
dfw <- dfw %>%
  mutate(status_500 = ifelse(!is.na(status_501) & status_501 == 1, 1, status_500))

# Finally, change any midpoint values to 0 when status is 0
dfw <- dfw %>%
  mutate(midpoint_500 = ifelse(!is.na(status_500) & status_500 == 0, 
                               0, midpoint_500)) %>%
  mutate(midpoint_501 = ifelse(!is.na(status_501) & status_501 == 0,
                               0, midpoint_501))

# Create table with location options
loc_options <- distinct(dfw, state, region) %>%
  mutate(region = ifelse(is.na(region), "Statewide", region)) %>%
  rbind(data.frame(state = "SC states (4)", 
                   region = "Entire SC region"))

# Desired minimum sample size for calculating proportions
min_sample_size <- 5

# Load iNat (via Phenobase) data ----------------------------------------------#
# iNat observations with flower annotation (don't have separate open flower
# annotation)

phenofolder <- "C:/Users/erin/Desktop/phenobase"
phenozips <- list.files(phenofolder, full.names = TRUE)

for (i in seq_along(phenozips)) {
  filename <- unz(phenozips[i], "data.csv")
  pheno1 <- read.csv(filename)
  if (i == 1) {
    inatfull <- pheno1
  } else {
    inatfull <- bind_rows(inatfull, pheno1)
  }
}
rm(pheno1)

# Simplify data
inat <- inatfull %>%
  select(scientificName, year, latitude, longitude, eventDate,
         coordinateUncertaintyInMeters) %>%
  mutate(date = ymd(eventDate)) %>%
  rename(coordU = coordinateUncertaintyInMeters,
         spp = scientificName) %>%
  select(-eventDate)

# Look at what points have tons of geographic uncertainty
# inat %>%
#   mutate(uncertain = ifelse(is.na(coordU) | coordU > 10000, 1, 0)) %>%
#   ggplot() +
#   geom_point(aes(x = longitude, y = latitude, color = factor(uncertain)))
# No obvious patterns

# What does it mean if the coordU field is NA?
# For now will delete these, along with records that have very high values
# (location uncertainty > 10 km)
inat <- inat %>%
  filter(!is.na(latitude)) %>%
  filter(!is.na(coordU) & coordU <= 10000) %>%
  select(-coordU)

# Load states shapefile (data/states)
# Assign states
# Filter to 4-state region
# Label wk, wk2, and month




