################################################################################
# Model variation in flowering onset dates

# Erin Zylstra
# 2024-07-24
################################################################################

require(dplyr)
require(lubridate)
require(stringr)
require(tidyr)
require(ggplot2)

rm(list = ls())

# Load relevant data ----------------------------------------------------------#

  # Flowering onset data, by plant-wateryr
  onset <- read.csv("data/flowering-onsets-2013-2023.csv")
  
  # Species-level summaries of flowering onset data (for species-phenophase
  # combinations that had a minimum of 30 plant-yrs of data)
  spp_onsets <- read.csv("output/flowering-onsets/flower-onset-data-summaries.csv")
 
# Select parameters for data inclusion/exclusion ------------------------------#
  
  # Set the maximum number of days prior to a "yes" that we'll need a "no" in 
  # order to use the observation in onset models
  daysprior <- 14    

# Reformat and subset onset data ----------------------------------------------#
  
  # Extract plant-wateryr data for flower onset
  onset_fl <- onset %>%
    select(-c(contains("fo"), first_obs, last_obs)) %>%
    rename_with(str_replace, pattern = "fl_", replacement = "") %>%
    rename_with(str_replace, pattern = "_fl", replacement = "") %>%
    mutate(phenophase = "flower") %>%
    relocate(phenophase, .after = "n_obs")
  
  # Extract plant-wateryr data for open flower onset
  onset_fo <- onset %>%
    select(-c(contains("fl"), first_obs, last_obs)) %>%
    rename_with(str_replace, pattern = "fo_", replacement = "") %>%
    rename_with(str_replace, pattern = "_fo", replacement = "") %>%
    mutate(phenophase = "open flower") %>%
    relocate(phenophase, .after = "n_obs") 
  
  # Combine the two dataframes. Now each row contains onset data for a unique 
  # combination of plant, wateryr, and phenophase
  onset <- rbind(onset_fl, onset_fo)
  
  # Remove plant-yrs when:
  # 1) flowers/open flowers weren't observed (firstyes = NA),
  # 2) "yes" to phenophase status wasn't preceded by a "no" (dayssince_lastno = NA)
  # 3) preceding "no" occurred more than XX days ago (dayssince_lastno > daysprior)
  onset <- onset %>%
    filter(!is.na(firstyes)) %>%
    filter(!is.na(dayssince_lastno)) %>%
    filter(dayssince_lastno <= daysprior)
  # count(onset, phenophase)

# Prioritize species to use in models -----------------------------------------#

  # Priority 1: 2 or more locations in SC region
  # Priority 2: 1 location in SC region
  # Priority 3: nearest plant location < 200 km from SC region
  # If a species-phenophase doesn't meet any of these criteria, probably won't 
    # create a model for onset date (and leaving priority = NA)

  spp_onsets <- spp_onsets %>%
    rowwise() %>%
    mutate(SC = sum(c_across(LA:TX))) %>%
    mutate(priority = case_when(
      SC > 1 ~ 1,
      SC == 1 ~ 2,
      dist_to_SC_km < 200 ~ 3,
      .default = NA
    )) %>%
    data.frame()

  # Remove species from onset data if they're not priority 1-3 (Will keep onset
  # data for both phenophases if the species is considered priority 1-3 for
  # at least one phenophase)
  spp_all <- spp_onsets %>%
    filter(!is.na(priority)) %>%
    select(common_name) %>%
    distinct() %>%
    unlist(use.names = FALSE)
  
  onset <- onset %>%
    filter(common_name %in% spp_all)
  
# Download daymet data (if needed) --------------------------------------------#  
  
  sites <- onset %>%
    select(site_id, lon, lat) %>%
    distinct() 
  
  nsites <- nrow(sites)
  yrs <- (min(onset$wateryr) - 1):max(onset$wateryr)
  
  daymet_csv <- paste0("data/climate/daymet-", nsites, "sites-", 
                       first(yrs), "-", last(yrs), ".csv")
  daymet_zip <- str_replace(daymet_csv, ".csv", ".zip")
  
  if (!file.exists(daymet_zip)) {
    source("download-daymet.R")
  }
  
  # TODO ###########################
  # Then unzip, load csv, remove csv
  # Summarize climate data by season...
  
  
# Test workflow/run for one species -------------------------------------------#
  
  # Important to remember that all day numbers in the onset dataset are actually
  # day-of-WATER-year (dowy). Creating simple table to look at important day 
  # numbers (no leap year).
  dowy <- data.frame(first_date = c(as.Date(paste0("2021-", 10:12, "-01")),
                                    as.Date(paste0("2022-", 1:9, "-01"))))
  dowy <- dowy %>%
    mutate(mon = month(first_date),
           n_days = days_in_month(first_date),
           dowy = as.numeric(first_date - make_date(2021, 9, 30)))
  dowy
    # Jan 1 = 93
    # Mar 1 = 152
    # Jun 1 = 244
    # Sep 1 = 336
  
  phenophases <- c("flower", "open flower")
  priority_levels <- 1:3
  
  # For now, picking common buttonbush, open flower (but structure is here
  # to do things in a loop)
  i = 2
  j = 1 
  pheno <- phenophases[i]
  spp_list <- spp_onsets %>%
    filter(phenophase == pheno) %>%
    filter(priority == priority_levels[j]) %>%
    select(common_name) %>%
    unlist(use.names = FALSE)
  
  spp <- spp_list[2]
  
  onsetsub <- onset %>%
    filter(phenophase == pheno) %>%
    filter(common_name == spp)

  # Identify when the plant flowers so we know what climate data to grab (to
  # ensure it's data preceding the event and not after). For the purposes of 
  # describing seasonal climate conditions, NPN has previously defined Spring = 
  # Mar-May; Summer data = Jun-Aug.
  median_onset <- median(onsetsub$firstyes)
  # If onset is before 1 Apr (dowy < 183), then use spring and summer weather 
  # from previous calendar year.
  # If onset is 1 Apr - 30 June (dowy = 183 - 273), then use spring weather in 
  # current year and summer weather in previous year
  # If onset is on or after 1 Jul (dowy = 274), then use spring and summer 
  # weather from current calendar year.
  spring <- ifelse(median_onset < 183, "previous", "current")
  summer <- ifelse(median_onset %in% 183:273, "previous", "current")
