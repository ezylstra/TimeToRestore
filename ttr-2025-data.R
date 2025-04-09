################################################################################
# Summarize and evaluate quality of 2025 TTR data
# Erin Zylstra
# 2025-04-08
################################################################################

require(rnpn)
require(dplyr)
require(lubridate)
require(stringr)
require(tidyr)
library(terra)
library(tidyterra)

rm(list = ls())

# Logical indicating whether to re-download data
update_data <- TRUE

# Create custom rounding function ---------------------------------------------#

# This is the same as plyr::round_any(), but loading the plyr package often 
# causes conflicts with function names when dplyr is loaded

round_any <- function(x, accuracy) {
  round(x / accuracy) * accuracy
}

# Load shapefiles with state boundaries ---------------------------------------#

# states <- vect("data/states/cb_2017_us_state_500k.shp")
# states <- subset(states, 
#                  !states$STUSPS %in% c("HI", "AK", "VI", "MP", 
#                                        "GU", "PR", "AS"))

# Identify priority plant species ---------------------------------------------#

# List of priority species in ttr-priorityspecies-20250402, which is based on: 
# TimeToRestore_AllPrioritySpecies_2024 google sheet

# Load list of priority species
spp_list <- read.csv("data/ttr-priorityspecies-20250402.csv", 
                     na.strings = c(NA, ""))
spp_list <- spp_list %>%
  mutate(across(c(common_name, scientific_name), str_trim)) %>%
  # Edit spelling of one genus (Vernonia) and one species (cespitosa)
  # Update species based on ITIS for mistflower (Conoclinium), 
  # basket-flower (Centaurea), and Lantana
  mutate(scientific_name = case_when(
    scientific_name == "Oenothera caespitosa" ~ "Oenothera cespitosa",
    scientific_name == "Veronia gigantea" ~ "Vernonia gigantea",
    scientific_name == "Conoclinium greggii" ~ "Conoclinium dissectum",
    scientific_name == "Centaurea americana" ~ "Plectocephalus americanus",
    scientific_name == "Lantana urticoides" ~ "Lantana horrida", 
    .default = scientific_name
  )) %>%
  rename(ttr_common_name = common_name)

# Load information about NN species
nn_spp <- npn_species() %>% data.frame()
nn_spp <- nn_spp %>%
  filter(kingdom == "Plantae") %>%
  select(species_id, common_name, genus, species, functional_type) %>%
  mutate(scientific_name = paste(genus, species))

# Find NN info based on scientific name of priority species (inconsistent 
# capitalization in priority list and a duplicate in NN database [Canada goldenrod])
spp_list <- spp_list %>%
  left_join(nn_spp, by = "scientific_name")
# Check that all priority species have a match in NN database
  # filter(spp_list, is.na(species_id))
# Are there any duplicates? Yes
  count(spp_list, species_id) %>% filter(n > 1)
  filter(spp_list, species_id == 90) 
  # In priority list, Sambucus nigra listed as both black and common elderberry
  filter(spp_list, species_id == 182) 
  # In priority list, Passiflora incarnata listed as both purple passionflower and Maypop

# Remove entries for TTR priority species whose common name doesn't match NPN
spp_list <- spp_list %>%
  filter(!ttr_common_name %in% c("Maypop", "Common elderberry"))
# 53 species

# Identify phenophases of interest --------------------------------------------#

# Focus on 4 phenophases: flowers, open flowers, fruit, ripe fruits
# 500 = Flowers or flower buds
# 501 = Open flowers
# 516 = Fruits
# 390 = Ripe fruits

# First, check that all species have these 4 phenophases
phenophases_byspp <- npn_phenophases_by_species(
  species_ids = c(spp_list$species_id),
  date = "2025-01-01"
) %>% data.frame()
phenophases_byspp %>%
  group_by(species_id, species_name) %>%
  summarize(p500 = ifelse(500 %in% phenophase_id, 1, 0),
            p501 = ifelse(501 %in% phenophase_id, 1, 0),
            p516 = ifelse(516 %in% phenophase_id, 1, 0),
            p390 = ifelse(390 %in% phenophase_id, 1, 0),
            .groups = "keep") %>%
  rowwise() %>%
  filter(sum(c_across(p500:p390)) < 4)
# All species use these 4 phenophases

# What phenophases do the milkweeds use?
milkweed_phps <- phenophases_byspp %>%
  filter(str_detect(species_name, "milkweed")) %>%
  select(species_id, species_name, pheno_class_id, phenophase_id, phenophase_name) %>%
  arrange(species_name, pheno_class_id, phenophase_id)
# Since leaves may be important for monarch eggs and catepillars, we may also 
# want to include:
# 488 = Leaves (for milkweeds only)

# Download and format (or load existing) NPN data for priority plant species --#
  
# Observations in 4 states (LA, NM, OK, TX)
# Focus on 2025 data, but also download 2024 data (for fruits and/or comparison)

phenophases <- c(500, 501, 516, 390)
states4 <- c("LA", "NM", "OK", "TX")
  # Minor note: we could be missing observations in the four states if we use
  # the states argument in the download function because sometimes the state 
  # field is missing or incorrect. Won't worry about that here, but would be 
  # better in the long run to check state values based on lat/lons

data_filename <- "data/ttr-data-20242025.csv"

if (!file.exists(data_filename) | update_data == TRUE) {
  
  # Download flowering, fruiting data for 2024-2025
  status_dl <- npn_download_status_data(
    request_source = "erinz",
    years = 2024:2025,
    species_ids = spp_list$species_id,
    phenophase_ids= phenophases,
    states = states4,
    additional_fields = c("observedby_person_id",
                          "partner_group",
                          "site_name", 
                          "species_functional_type"))
  status_dl <- data.frame(status_dl)
  
  # Download leafing data for milkweeds in 2024-2025
  milkweeds <- spp_list %>%
    filter(str_detect(common_name, "milkweed")) %>%
    pull(species_id)
  status_mwleaf_dl <- npn_download_status_data(
    request_source = "erinz",
    years = 2024:2025,
    species_ids = milkweeds,
    phenophase_ids= 488,
    states = states4,
    additional_fields = c("observedby_person_id",
                          "partner_group",
                          "site_name", 
                          "species_functional_type"))
  status_mwleaf_dl <- data.frame(status_mwleaf_dl)
  
  # Combine everything and format
  status_df <- rbind(status_dl, status_mwleaf_dl) %>%
    mutate(obsdate = ymd(observation_date),
           yr = year(obsdate),
           php = case_when(
             phenophase_id == 500 ~ "flower",
             phenophase_id == 501 ~ "open flower",
             phenophase_id == 516 ~ "fruit",
             phenophase_id == 390 ~ "ripe fruit",
             phenophase_id == 488 ~ "leaves")) %>%
    select(-c(update_datetime, elevation_in_meters, genus, species, kingdom,
              phenophase_description, abundance_value, observation_date)) %>%
    rename(person_id = observedby_person_id,
           func_type = species_functional_type,
           lat = latitude,
           lon = longitude)
  
  # Write to file
  write.csv(status_df, data_filename, row.names = FALSE)
  # Remove objects
  rm(status_df, status_dl, status_mwleaf_dl)
  
}

status <- read.csv(data_filename)

# Append information about intensity categories/values ------------------------#

# Download information about intensity categories
ic <- npn_abundance_categories() %>% data.frame()
ic <- ic %>%
  rename(intensity_category_id = category_id, 
         intensity_value_id = value_id,
         intensity_name = category_name,
         intensity_value = value_name) %>%
  select(-c(category_description, value_description))

# Extract just those categories that appear in status data and format:
ic_subset <- ic %>%
  filter(intensity_category_id %in% unique(status$intensity_category_id)) %>%
  mutate(value1 = NA,
         value2 = NA,
         intensity_type = case_when(
           str_detect(intensity_value, "%") ~ "percent",
           str_detect(intensity_value, "[0-9]") ~ "number",
           .default = "qualitative"
         ))
val12 <- which(colnames(ic_subset) %in% c("value1", "value2"))
for (i in 1:nrow(ic_subset)) {
  if (str_detect(ic_subset$intensity_value[i], " to ")) {
    ic_subset[i, val12] <- str_split_fixed(ic_subset$intensity_value[i], " to ", 2)
    ic_subset[i, val12] <- as.numeric(str_remove(ic_subset[i, val12], ","))
  } else if (str_detect(ic_subset$intensity_value[i], "-")) {
    ic_subset[i, val12] <- str_split_fixed(ic_subset$intensity_value[i], "-", 2)
    ic_subset[i, val12[2]] <- str_remove(ic_subset[i, val12[2]], "%")
  } else if (str_detect(ic_subset$intensity_value[i], "% or more")) {
    ic_subset[i, val12] <- str_remove(ic_subset$intensity_value[i], "% or more")
  } else if (str_detect(ic_subset$intensity_value[i], "Less than ")) {
    ic_subset[i, val12[1]] <- 0
    ic_subset[i, val12[2]] <- str_remove(ic_subset$intensity_value[i], "Less than ")
    ic_subset[i, val12[2]] <- str_remove(ic_subset[i, val12[2]], "%")
  } else if (str_detect(ic_subset$intensity_value[i], "More than ")) {
    ic_subset[i, val12] <- str_remove(ic_subset$intensity_value[i], "More than ")
    ic_subset[i, val12[1]] <- as.numeric(str_remove(ic_subset[i, val12[1]], ",")) + 1
    ic_subset[i, val12[2]] <- as.numeric(str_remove(ic_subset[i, val12[2]], ",")) + 1
  }
}
ic_subset <- ic_subset %>%
  mutate_at(c("value1", "value2"), as.numeric)

# Assigning a middle-ish value for each range (keeping it to nice numbers like 
# 5, 50, 500, and 5000)
ic_subset <- ic_subset %>%
  mutate(mag = nchar(value1) - 1) %>%
  mutate(value = case_when(
    value1 == value2 ~ round(value1),
    intensity_type == "number" & value1 == 0 ~ 1,
    intensity_type == "number" & value1 != 0 ~ 
      round_any(rowMeans(across(value1:value2)), 5 * (10 ^ mag)),
    intensity_type == "percent" ~ round(rowMeans(across(value1:value2))),
    .default = NA
  )) %>%
  select(-c(mag, value1, value2))

ic_append <- ic_subset %>%
  select(intensity_category_id, intensity_name, intensity_value, value) %>%
  rename(intensity_cat = intensity_value, 
         intensity = value)

status <- status %>%
  left_join(ic_append, 
            by = c("intensity_category_id", 
                   "intensity_value" = "intensity_cat")) %>%
  select(-intensity_category_id)
  
# Duplicate observations ------------------------------------------------------###########################

# Any duplicate observations by the same person (same plant, date, phenophase
# status and intensity value)?
dups <- status %>%
  group_by(person_id, individual_id, obsdate, yr, phenophase_id, 
           phenophase_status, intensity) %>%
  summarize(nobs = n(), .groups = "keep") %>%
  data.frame() %>%
  count(nobs)
dups
n_dups <- sum(dups$n * (dups$nobs - 1))
n_dups / nrow(status) * 100 
# 0.24% (n = 39) of all observations in 2024-2025 are duplicates

dups25 <- status %>%
  filter(yr == 2025) %>%
  group_by(person_id, individual_id, obsdate, yr, phenophase_id, 
           phenophase_status, intensity) %>%
  summarize(nobs = n(), .groups = "keep") %>%
  data.frame() %>%
  count(nobs)
dups25
n_dups25 <- sum(dups25$n * (dups25$nobs - 1))
n_dups25 / nrow(status[status$yr == 2025,]) * 100 
# 0.68% (n = 31) observations in 2025 are duplicates

# Removing duplicates (will cause problems later in script)
status <- status %>% distinct(across(-observation_id))

# Any duplicate observations by different people (same plant, date, phenophase
# status and intensity value)?
dups2 <- status %>%
  group_by(individual_id, obsdate, phenophase_id, phenophase_status,
           intensity) %>%
  summarize(nobs = n(), .groups = "keep") %>%
  data.frame() %>%
  count(nobs)
dups2
n_dups2 <- sum(dups2$n * (dups2$nobs - 1))
n_dups2 / nrow(status) * 100
# 2.7% of 2024-2025 observations (could come back to this and look for 
# conflicting reports, but this isn't a priority right now)

# Problems related to intensity values ----------------------------------------###########################

# Are there instances where there is an intensity value but phenophase status is
# not yes?
status <- status %>%
  mutate(intensity_NA = ifelse(is.na(intensity), 1, 0),
         intensity_NotYes = ifelse(intensity_NA == 0 & phenophase_status != 1, 
                                   1, 0))
filter(status, intensity_NotYes == 1)
  # Just one instance in 2024 (American beautyberry; TX)

# Are people reporting the number of flowers instead of inflorescences?
# Highest two intensity values for each intensity category worth investigating
ic_subset <- ic_subset %>%
  group_by(intensity_category_id) %>%
  mutate(intensity_high = ifelse(value %in% tail(sort(value), 2), 1, 0)) %>%
  ungroup() %>%
  mutate(intensity_high = ifelse(!str_detect(intensity_name, 
                                             "Flowers and flower buds"), 
                                 NA, intensity_high)) %>%
  data.frame()
status <- status %>%
  left_join(select(ic_subset, intensity_name, intensity_value, intensity_high),
            by = c("intensity_name", "intensity_value"))

# Table with number and proportion of observations in top 2 highest intensity
# categories by species and year
intensity_high <- status %>%
  group_by(common_name, php) %>%
  summarize(n2024 = sum(yr == 2024),
            n2025 = sum(yr == 2025),
            nhigh_2024 = sum(yr == 2024 & intensity_high == 1 & !is.na(intensity_high)),
            nhigh_2025 = sum(yr == 2025 & intensity_high == 1 & !is.na(intensity_high)),
            .groups = "keep") %>% 
  filter(nhigh_2024 + nhigh_2025 != 0) %>%
  mutate(prophigh_2024 = round(nhigh_2024 / n2024, 2), 
         prophigh_2025 = round(nhigh_2025 / n2025, 2)) %>%
  data.frame()
intensity_high
# Table with the number of observations in each high intensity category by
# species and year
intensitycat_high <- status %>%
  filter(intensity_high == 1) %>%
  group_by(common_name, php, intensity_value) %>%
  summarize(nhigh_2024 = sum(yr == 2024),
            nhigh_2025 = sum(yr == 2025),
            nplants_2024 = n_distinct(individual_id[yr == 2024]),
            nplants_2025 = n_distinct(individual_id[yr == 2025]),
            .groups = "keep") %>%
  data.frame()
intensitycat_high

# Phenophase status inconsistencies -------------------------------------------###########################

# To look at this, can't have more than one observation of a plant per person
# per day. We've already removed duplicates, but now need to resolve instances 
# where somebody made multiple observations of the same plant on the same date 
# that differed in some way.

  # For now, will keep record with more advanced phenophase or higher 
  # intensity value. Will do this by sorting observations in descending
  # order and keeping only the first
  inddateobsp <- status %>%
    group_by(common_name, individual_id, obsdate, person_id, php) %>%
    summarize(n_obs = n(),
              .groups = "keep") %>%
    data.frame()
  inddateobsp$obsnum <- 1:nrow(inddateobsp)
  
  status <- status %>%
    arrange(person_id, individual_id, obsdate, php, 
            desc(phenophase_status), desc(intensity)) %>%
    left_join(select(inddateobsp, -c(n_obs, common_name)), 
              by = c("person_id", "individual_id", 
                     "obsdate", "php")) %>%
    # Create "dups" column, where dups > 1 indicates that the observation can be
    # removed since there's another observation that same day with more advanced
    # phenology or higher intensity/abundance.
    mutate(dups = sequence(rle(as.character(obsnum))$lengths))
  
  # Remove extra observations and unnecessary columns
  status <- status %>%
    filter(dups == 1) %>%
    select(-c(obsnum, dups)) %>%
    arrange(common_name, obsdate, person_id, php)

# To identify inconsistent status values, will need to put flower/fruit data 
# into wide form (all data for a plant visit in the same row). Removing
# unknown status observations first (<0.5% of fruit/flower observations).
statusw <- status %>%
  filter(php != "leaves") %>%
  filter(phenophase_status != -1) %>%
  select(person_id, partner_group, site_id, state, common_name, individual_id,
         yr, obsdate, php, phenophase_status, intensity) %>%
  rename(status = phenophase_status) %>%
  pivot_wider(names_from = php,
              names_glue = "{php}_{.value}",
              values_from = c(status, intensity)) %>%
  # Add month, day to identify observations in water year
  mutate(month = month(obsdate),
         day = day(obsdate)) %>%
  data.frame()

# Identify phenophase status inconsistencies
# NOTE: changing NAs to 999 in order to make this code simpler
statusw <- statusw %>%
  mutate(across(contains("status"), ~replace_na(., 999))) %>%
  # Problem: flower = 0, open = NA or 1
  mutate(flower0_openNot0 = ifelse(flower_status == 0 & open.flower_status != 0, 
                                   1, 0)) %>%
  # Problem: flower = NA, open = 1
  mutate(flowerNA_open1 = ifelse(flower_status == 999 & open.flower_status == 1,
                                 1, 0)) %>%
  # Problem: fruit = 0, ripe = NA or 1
  mutate(fruit0_ripeNot0 = ifelse(fruit_status == 0 & ripe.fruit_status != 0, 
                                   1, 0)) %>%
  # Problem: fruit = NA, ripe = 1
  mutate(fruitNA_ripe1 = ifelse(fruit_status == 999 & ripe.fruit_status == 1,
                                1, 0))

# Table summarizing problems with flower/open flower status by calendar year
flower_probs <- statusw %>%
  group_by(common_name, yr) %>%
  summarize(n = n(),
            n_flower0 = sum(flower_status == 0),
            n_flower1 = sum(flower_status == 1),
            n_flowerNA = sum(flower_status == 999),
            n_open0 = sum(open.flower_status == 0),
            n_open1 = sum(open.flower_status == 1),
            n_openNA = sum(open.flower_status == 999),
            n_flower0_openNot0 = sum(flower0_openNot0 == 1),
            n_flowerNA_open1 = sum(flowerNA_open1 == 1),
            .groups = "keep") %>%
  filter(n_flower0_openNot0 + n_flowerNA_open1 > 0) %>%
  data.frame()
flower_probs

# Table summarizing problems with fruit/ripe fruit status by calendar year
fruit_probs <- statusw %>%
  group_by(common_name, yr) %>%
  summarize(n = n(),
            n_fruit0 = sum(fruit_status == 0),
            n_fruit1 = sum(fruit_status == 1),
            n_fruitNA = sum(fruit_status == 999),
            n_ripe0 = sum(ripe.fruit_status == 0),
            n_ripe1 = sum(ripe.fruit_status == 1),
            n_ripeNA = sum(ripe.fruit_status == 999),
            n_fruit0_ripeNot0 = sum(fruit0_ripeNot0 == 1),
            n_fruitNA_ripe1 = sum(fruitNA_ripe1 == 1),
            .groups = "keep") %>%
  filter(n_fruit0_ripeNot0 + n_fruitNA_ripe1 > 0) %>%
  data.frame()
fruit_probs

# Table summarizing problems with flower or fruiting status since 1 Oct 2024
status_probs <- statusw %>%
  filter(yr == 2025 | (yr == 2024 & month >= 10)) %>%
  group_by(common_name) %>%
  summarize(n_flower0_openNot0 = sum(flower0_openNot0 == 1),
            n_flowerNA_open1 = sum(flowerNA_open1 == 1),
            n_fruit0_ripeNot0 = sum(fruit0_ripeNot0 == 1),
            n_fruitNA_ripe1 = sum(fruitNA_ripe1 == 1),
            .groups = "keep") %>%
  filter(n_flower0_openNot0 + n_flowerNA_open1 + n_fruit0_ripeNot0 + n_fruitNA_ripe1 > 0) %>%
  data.frame()
status_probs

# Reporting no prior to yes ---------------------------------------------------###########################

# Focusing on observations made this water year (1 Oct 2024 - present), so we
# can download individual phenometrics data through rnpn

ip_filename <- "data/ttr-ipdata-20242025.csv"

if (!file.exists(ip_filename) | update_data == TRUE) {
  
  # Download flowering, fruiting data for Oct 2024 - April 2025 and Oct 2023 - 
  # Sep 2024 (phenometrics aggregated over water year).
  # After some experimenting, it looks like when requesting data by water year, 
  # the years argument corresponds to the start of each year in October. So we
  # need to use years = 2023:2024 for this call.
  ip_dl <- npn_download_individual_phenometrics(
    request_source = "erinz",
    years = 2023:2024,
    period_start = "10-01",
    period_end = "09-30",
    species_ids = spp_list$species_id,
    phenophase_ids= phenophases,
    states = states4,
    additional_fields = c("observedby_person_id",
                          "partner_group",
                          "site_name", 
                          "species_functional_type"))
  ip_dl <- data.frame(ip_dl)
  
  # Download leafing data for milkweeds in 2024-2025
  milkweeds <- spp_list %>%
    filter(str_detect(common_name, "milkweed")) %>%
    pull(species_id)
  ip_mwleaf_dl <- npn_download_individual_phenometrics(
    request_source = "erinz",
    years = 2023:2024,
    period_start = "10-01",
    period_end = "09-30",
    species_ids = milkweeds,
    phenophase_ids= 488,
    states = states4,
    additional_fields = c("observedby_person_id",
                          "partner_group",
                          "site_name", 
                          "species_functional_type"))
  ip_mwleaf_dl <- data.frame(ip_mwleaf_dl)
  
  # Combine everything and format
  ip_df <- rbind(ip_dl, ip_mwleaf_dl) %>%
    mutate(php = case_when(
      phenophase_id == 500 ~ "flower",
      phenophase_id == 501 ~ "open flower",
      phenophase_id == 516 ~ "fruit",
      phenophase_id == 390 ~ "ripe fruit",
      phenophase_id == 488 ~ "leaves")) %>%
    select(-c(observedby_person_id, elevation_in_meters, genus, species, kingdom,
              phenophase_description, first_yes_month, first_yes_day,
              first_yes_julian_date, last_yes_year, last_yes_month, 
              last_yes_day, last_yes_julian_date, numdays_until_next_no)) %>%
    rename(func_type = species_functional_type,
           lat = latitude,
           lon = longitude,
           prior_no = numdays_since_prior_no)
  
  # Write to file
  write.csv(ip_df, ip_filename, row.names = FALSE)
  # Remove objects
  rm(ip_df, ip_dl, ip_mwleaf_dl)
  
}

ip <- read.csv(ip_filename)


# Summary of data collected ---------------------------------------------------#

# Map (just for 2025? both years?)
# Species (table or map or both) - organize by year and functional type?
# Number of observers (though not everyone may be associated with TTR)

