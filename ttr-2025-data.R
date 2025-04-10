################################################################################
# Summarize and evaluate quality of 2025 TTR data
# Erin Zylstra
# 2025-04-08
################################################################################

require(rnpn)
require(lubridate)
require(stringr)
require(tidyr)
library(ggplot2)
library(maps)
# library(mapdata)
require(dplyr)

rm(list = ls())

# Logical indicating whether to re-download data
update_data <- FALSE

# Create custom rounding function ---------------------------------------------#

# This is the same as plyr::round_any(), but loading the plyr package often 
# causes conflicts with function names when dplyr is loaded

round_any <- function(x, accuracy) {
  round(x / accuracy) * accuracy
}

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
  
# Duplicate observations ------------------------------------------------------#

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

# Problems related to intensity values ----------------------------------------#

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

# Phenophase status inconsistencies -------------------------------------------#

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

# Sites monitored this water year ---------------------------------------------#

# Plot locations where flowering/fruiting or milkweed leaves observed in current
# water year (Oct 2024 and present)
status <- status %>%
  mutate(current_wy = ifelse(obsdate >= "2024-10-01", 1, 0)) %>%
  mutate(ind_date = paste0(individual_id, "_", obsdate)) 

locs <- status %>%
  filter(current_wy == 1) %>%
  group_by(site_id, lat, lon, state) %>%
  summarize(nspp = n_distinct(common_name),
            nplants = n_distinct(individual_id),
            nobservers = n_distinct(person_id),
            # nobs: Number of observations, where an observations is all
            # data (all phenophases) submitted for a plant on given date
            nobs = n_distinct(ind_date), 
            .groups = "keep") %>%
  data.frame()

# Get state, county boundaries from map, mapdata packages
state <- map_data("state")
counties <- map_data("county")
ttr_states <- subset(state, region %in% c("louisiana", "new mexico", 
                                          "oklahoma", "texas"))
ttr_counties <- subset(counties, region %in% c("louisiana", "new mexico", 
                                               "oklahoma", "texas"))

ggplot(data = ttr_counties, aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3) +
  geom_polygon(color = "gray", fill = NA, linewidth = 0.2) + 
  geom_polygon(data = ttr_states, color = "black", fill = NA) +
  geom_point(data = locs, aes(x = lon, y = lat, group = NA, size = nplants), 
             color = "blue") +
  labs(size = "No. plants") +
  guides(size = guide_legend(position = "inside")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.justification.inside = c(0.05, 0.05))

# Plants monitored this water year --------------------------------------------#

plants <- status %>%
  filter(current_wy == 1) %>%
  group_by(common_name, func_type, individual_id, site_id, state) %>%
  summarize(nvisits = n_distinct(obsdate),
            nobservers = n_distinct(person_id),
            .groups = "keep") %>%
  data.frame()

plantspp <- plants %>%
  group_by(common_name, func_type) %>%
  summarize(nplants = n(),
            nsites = n_distinct(site_id),
            nplants_LA = sum(state == "LA"),
            nplants_NM = sum(state == "NM"),
            nplants_OK = sum(state == "OK"),
            nplants_TX = sum(state == "TX"),
            mean_no_visits = round(mean(nvisits), 1),
            .groups = "keep") %>%
  data.frame() %>%
  arrange(func_type, desc(nplants))

# Observers submitting data this water year -----------------------------------#

# Why are there some negative person_ids? (associated with NEON).
# Remove them from summary? Some make observations in two states.

observer_spp <- status %>%
  filter(current_wy == 1) %>%
  group_by(person_id, common_name) %>%
  summarize(nplants = n_distinct(individual_id),
            nvisits = n_distinct(ind_date),
            .groups = "keep") %>%
  data.frame()

observers <- observer_spp %>%
  group_by(person_id) %>%
  summarize(nspp = n_distinct(common_name),
            totalplants = sum(nplants),
            plants_per_spp = round(mean(nplants), 2),
            totalvisits = sum(nvisits)) %>%
  data.frame()

# Considering an observer different in different state
observers_state <- status %>%
  filter(current_wy == 1) %>%
  mutate(observer = paste0(state, person_id)) %>%
  group_by(observer, state, common_name) %>%
  summarize(nplants = n_distinct(individual_id),
            nvisits = n_distinct(ind_date),
            .groups = "keep") %>%
  data.frame() %>%
  group_by(observer, state) %>%
  summarize(nspp = n_distinct(common_name),
            totalplants = sum(nplants),
            plants_per_spp = round(mean(nplants), 2),
            totalvisits = sum(nvisits),
            .groups = "keep") %>%
  data.frame()

observers_state_s <- observer_state %>%
  group_by(state) %>%
  summarize(n_observers = n_distinct(observer),
            spp_per_obs_mn = round(mean(nspp), 2),
            spp_per_obs_max = max(nspp),
            plants_per_obs_mn = round(mean(totalplants), 2),
            plants_per_obs_max = max(totalplants),
            plants_per_spp_mn = round(mean(plants_per_spp), 2),
            plants_per_spp_max = max(plants_per_spp)) %>%
  data.frame()
observers_state_s

# Histograms
obs_spp_bar <- count(observers, nspp)
p_nspp <- ggplot(obs_spp_bar) +
  geom_col(aes(x = nspp, y = n), width = 0.75, fill = "steelblue3") +
  scale_x_continuous(breaks = seq(min(obs_spp_bar$nspp), max(obs_spp_bar$nspp))) +
  labs(x = "Number of species each observer monitored",
       y = "Number of observers") +
  theme_bw()
p_nspp

obs_plants_bar <- count(observer_spp, nplants)
p_plantsperspp <- ggplot(obs_plants_bar) +
  geom_col(aes(x = nplants, y = n), width = 0.75, fill = "steelblue3") +
  scale_x_continuous(breaks = seq(min(obs_plants_bar$nplants), 
                                  max(obs_plants_bar$nplants))) +
  labs(x = "Number of plants per species",
       y = "Number of observer-species combinations") +
  theme_bw()
p_plantsperspp


# Evaluating observation frequency (no before first yes) ---------------------#

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

# Will look at observations from:
# 1 Nov 2024 - present (since that'll allow for first yeses with prior no within 30 days) 
# and for comparison: 1 Nov 2023 - 30 Sep 2024

ip <- ip %>%
  # Create "season" (Nov 2024 - April 2025 = 2025)
  mutate(yes_date = parse_date_time(paste(first_yes_year, first_yes_doy),
                                    orders = "Y j"),
         yes_month = month(yes_date),
         season = case_when(
           yes_month %in% 11:12 ~ first_yes_year + 1,
           yes_month %in% 1:9 ~ first_yes_year,
           .default = NA))
# check
# count(ip, first_yes_year, yes_month, season)

# Proportion of first yeses that are preceded by a prior no within X days, 
# by phenophase and season
obsfreq_php <- ip %>%
  filter(!is.na(season)) %>%
  mutate(php = factor(php, levels = c("flower", "open flower", "fruit", 
                                      "ripe fruit", "leaves"))) %>%
  group_by(season, php) %>%
  summarize(nobs = n(),
            nobs3 = sum(!is.na(prior_no) & prior_no <= 3),
            nobs7 = sum(!is.na(prior_no) & prior_no <= 7),
            nobs14 = sum(!is.na(prior_no) & prior_no <= 14),
            nobs30 = sum(!is.na(prior_no) & prior_no <= 30),
            .groups = "keep") %>%
  mutate(prop3 = round(nobs3/nobs, 2),
         prop7 = round(nobs7/nobs, 2),
         prop14 = round(nobs14/nobs, 2),
         prop30 = round(nobs30/nobs, 2)) %>%
  data.frame()
obsfreq_php

# Make a bar chart for this (minus milkweed leaf observations because few of them)
obsfreq_phpl <- ip %>%
  filter(!is.na(season)) %>%
  filter(php != "leaves") %>%
  mutate(php = factor(php, levels = c("flower", "open flower", 
                                      "fruit", "ripe fruit"))) %>%
  group_by(season, php) %>%
  summarize(nobs03 = sum(!is.na(prior_no) & prior_no %in% 1:3),
            nobs07 = sum(!is.na(prior_no) & prior_no %in% 4:7),
            nobs14 = sum(!is.na(prior_no) & prior_no %in% 8:14),
            nobs30 = sum(!is.na(prior_no) & prior_no %in% 15:30),
            nobs99 = sum((!is.na(prior_no) & prior_no > 30) | is.na(prior_no)),
            .groups = "keep") %>%
  pivot_longer(cols = nobs03:nobs99,
               names_to = "cat",
               values_to = "nobs") %>%
  data.frame()
obsfreq_phpl

ggplot(obsfreq_phpl) +
  geom_bar(aes(x = php, y = nobs, fill = cat),
           position = "stack",
           stat = "identity") +
  scale_fill_manual(values = c("#7fc97f", "#beaed4", "#fdc086",
                               "#ffff99", "#386cb0"),
                    labels = c("1-3 days", "4-7 days", "8-14 days",
                               "15-30 days", ">30 days")) +
  facet_grid( ~ season) +
  labs(x = "Phenophase", y = 'Number of first "yes" observations',
       fill = 'Prior "no"') +
  theme_bw()

# Proportion of first yeses that are preceded by a prior no within X days, 
# by species and season
obsfreq_spp <- ip %>%
  filter(!is.na(season)) %>%
  group_by(common_name, season) %>%
  summarize(nobs = n(),
            nobs3 = sum(!is.na(prior_no) & prior_no <= 3),
            nobs7 = sum(!is.na(prior_no) & prior_no <= 7),
            nobs14 = sum(!is.na(prior_no) & prior_no <= 14),
            nobs30 = sum(!is.na(prior_no) & prior_no <= 30),
            .groups = "keep") %>%
  mutate(prop3 = round(nobs3/nobs, 2),
         prop7 = round(nobs7/nobs, 2),
         prop14 = round(nobs14/nobs, 2),
         prop30 = round(nobs30/nobs, 2)) %>%
  data.frame()
obsfreq_spp %>% filter(nobs >= 10)

# Make a bar chart for this
# Only using species that had at least 10 first yeses in one of the years:
spp10 <- unique(obsfreq_spp$common_name[obsfreq_spp$nobs >= 10])
obsfreq_sppl <- ip %>%
  filter(!is.na(season)) %>%
  filter(php != "leaves") %>%
  filter(common_name %in% spp10) %>%
  mutate(common_name = factor(common_name)) %>%
  group_by(season, common_name) %>%
  summarize(nobs03 = sum(!is.na(prior_no) & prior_no %in% 1:3),
            nobs07 = sum(!is.na(prior_no) & prior_no %in% 4:7),
            nobs14 = sum(!is.na(prior_no) & prior_no %in% 8:14),
            nobs30 = sum(!is.na(prior_no) & prior_no %in% 15:30),
            nobs99 = sum((!is.na(prior_no) & prior_no > 30) | is.na(prior_no)),
            .groups = "keep") %>%
  pivot_longer(cols = nobs03:nobs99,
               names_to = "cat",
               values_to = "nobs") %>%
  data.frame()
obsfreq_sppl

ggplot(obsfreq_sppl) +
  geom_bar(aes(x = common_name, y = nobs, fill = cat),
           position = "stack",
           stat = "identity") +
  scale_fill_manual(values = c("#7fc97f", "#beaed4", "#fdc086",
                               "#ffff99", "#386cb0"),
                    labels = c("1-3 days", "4-7 days", "8-14 days",
                               "15-30 days", ">30 days")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
  facet_grid(rows = vars(season)) +
  labs(x = "Species", y = 'Number of first "yes" observations',
       fill = 'Prior "no"') +
  theme_bw() +
  theme(legend.position = "bottom")
