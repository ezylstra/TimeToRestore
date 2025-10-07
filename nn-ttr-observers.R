library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(tidyr)

# Need to access the NN database to get summaries of observation effort for 
# TTR partipants in Texas. This is a bit challenging because observers are not
# easily identified (ie, there is no LPP or similar group for people associated
# with TTR)

# Restrict summaries to water year (WY) 2024 and 2025 (Oct 2023 - Sep 2025)
# and TTR priority species in Texas
# Restrict to persons who have made observations of TTR priority species

# All the following for each WY and across both WYs?
# Number of records (record = status/int of plant, phenophase on given date)
# Number of observations (obs = all status/int for plant on given date)
# Number of plants
# Number of sites
# Mean observation frequency (for a given plant)

# Faster to do all the queries in R with downloaded data and then extract 
# observer table (with emails, names) from database?

# Specify how we'll obtain table with observer info ---------------------------#
# If access_via_R = FALSE, will download table via DBeaver and import. If TRUE
# using packages in R to communicate directly with database
access_via_R <- TRUE

# Download TTR data for Texas, 2023-2025 (if not done already) ----------------#
tx_file <- "data/ttr-data-2023-2025-tx.csv"

if(!file.exists(tx_file)) {

  library(rnpn)
  library(sf)
  library(terra)
  
  states <- vect("data/states/cb_2017_us_state_500k.shp")
  
  # Load and format TTR priority species
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
  # count(spp_list, species_id) %>% filter(n > 1)
  # filter(spp_list, species_id == 90)
  # In priority list, Sambucus nigra listed as both black and common elderberry
  # filter(spp_list, species_id == 182)
  # In priority list, Passiflora incarnata listed as both purple passionflower and Maypop
  
  # Remove entries for TTR priority species whose common name doesn't match NPN
  spp_list <- spp_list %>%
    filter(!ttr_common_name %in% c("Maypop", "Common elderberry"))
  
  # Simplify
  spp_simple <- spp_list %>%
    select(ttr_common_name, common_name, scientific_name, species_id)
  
  # Download status-intensity data for priority species, years
  
  # Note: can't download data for just Texas right now since some recent state 
  # assignments and missing in the database. So this download includes all
  # observations of priority TTR species in calendar years 2023-2025
  
  status_dl <- npn_download_status_data(
    request_source = "erinz",
    years = 2023:2025,
    species_ids = spp_simple$species_id,
    additional_fields = c("observedby_person_id",
                          "site_name"))
  status_dl <- data.frame(status_dl) %>%
    rename(person_id = observedby_person_id,
           lat = latitude,
           lon = longitude)
  
  # Filter observations by state (just keep TX)
  state_fill <- status_dl %>%
    select(site_id, lon, lat, state) %>%
    distinct()
  state_fillv <- vect(state_fill, 
                      geom = c("lon", "lat"), 
                      crs = "epsg:4326")
  state_new <- terra::extract(states, state_fillv)
  state_fill <- cbind(state_fill, state_new = state_new$STUSPS)
  # check:
  # count(state_fill, state, state_new) %>%
  #   mutate(same = ifelse(state == state_new, 1, 0)) %>%
  #   arrange(same)
  
  # Attach new state labels and keep TX observations only
  tx_df <- status_dl %>%
    left_join(select(state_fill, site_id, state_new), by = "site_id") %>%
    select(-state) %>%
    rename(state = state_new) %>%
    filter(!is.na(state) & state == "TX")
  
  # Clean up and save to file
  tx <- tx_df %>%
    select(-c(observation_id, update_datetime, elevation_in_meters, genus,
              species, common_name, kingdom, abundance_value)) %>%
    left_join(spp_simple, by = "species_id") %>%
    rename(nn_common_name = common_name)
  
  # Write to file
  write.csv(tx, "data/ttr-data-2023-2025-tx.csv", row.names = FALSE)
}

# Summarize data by observer --------------------------------------------------#

# Load Texas data
tx <- read.csv(tx_file)

# Remove any duplicate rows
tx <- tx[!duplicated(tx),]

# Create water year variable
tx <- tx %>%
  mutate(obsdate = ymd(observation_date),
         yr = year(obsdate),
         month = month(obsdate)) %>%
  mutate(wy = ifelse(month %in% 10:12, yr + 1, yr)) %>%
  select(-c(observation_date, day_of_year, intensity_category_id))

# Remove observations with negative values for person_ID because they're NEON
# personnel (and not associated with TTR)
tx <- tx %>%
  filter(person_id > 0)

# Are there any plants observed by the same person on the same day (that aren't 
# true duplicates)? Yes, so will keep only the observation with the most 
# advanced status or non-NA intensity value
tx <- tx %>%
  arrange(person_id, individual_id, phenophase_id, obsdate, 
          desc(phenophase_status), intensity_value) %>%
  distinct(person_id, individual_id, phenophase_id, obsdate, .keep_all = TRUE)

# Number of records per person
recs <- tx %>%
  filter(wy %in% 2024:2025) %>%
  group_by(person_id) %>%
  summarize(records = n(),
            records_2024 = sum(wy == 2024),
            records_2025 = sum(wy == 2025),
            sites = n_distinct(site_id)) %>%
  data.frame()

# Create dataframe with summarized info for person-plant-date:
obs <- tx %>%
  group_by(person_id, nn_common_name, individual_id, obsdate, wy) %>%
  summarize(n_php = n_distinct(phenophase_id), .groups = "keep") %>%
  data.frame() %>%
  filter(wy != 2026) %>%
  arrange(person_id, individual_id, obsdate)

# Calculate intervals between observations of same plant
obs$interval_prior <- NA 
for (i in 2:nrow(obs)) {
  obs$interval_prior[i] <- ifelse(
    obs$person_id[i] == obs$person_id[i-1] & 
      obs$individual_id[i] == obs$individual_id[i-1],
    as.numeric(obs$obsdate[i] - obs$obsdate[i-1]), NA
  )
}

# Summarize information by person-plant-water year (ppw)
ppw <- obs %>%
  filter(wy %in% 2024:2025) %>%
  group_by(person_id, nn_common_name, individual_id, wy) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  data.frame()
ppw <- ppw %>%
  group_by(person_id, nn_common_name, individual_id, wy, n_obs) %>%
  summarize(interval_mn = ifelse(n_obs[1] == 1, NA, mean(interval_prior, na.rm = TRUE)),
            interval_mx = ifelse(n_obs[1] == 1, NA, max(interval_prior, na.rm = TRUE)),
            .groups = "keep") %>%
  data.frame()

# Summarize across plants each year (pw = person-wateryr)
pw <- ppw %>%
  group_by(person_id, wy) %>%
  summarize(spp = n_distinct(nn_common_name),
            plants = n_distinct(individual_id),
            mn_dates_per_plant = round(mean(n_obs), 1),
            mn_interval = round(mean(interval_mn, na.rm = TRUE), 1),
            .groups = "keep") %>%
  data.frame()
# Put information into wide form and join with totals across years
pw_wide <- pw %>%
  pivot_wider(id_cols = person_id, 
              names_from = wy,
              values_from = c(spp, plants, mn_dates_per_plant, mn_interval)) %>%
  select(person_id, 
         spp_2024, plants_2024, mn_dates_per_plant_2024, mn_interval_2024,
         spp_2025, plants_2025, mn_dates_per_plant_2025, mn_interval_2025) %>%
  mutate(across(spp_2024:plants_2024, ~replace_na(., 0))) %>%
  mutate(across(mn_dates_per_plant_2024:mn_interval_2024, ~replace_na(., NA))) %>%
  mutate(across(spp_2025:plants_2025, ~replace_na(., 0))) %>%
  mutate(across(mn_dates_per_plant_2025:mn_interval_2025, ~replace_na(., NA))) %>%
  data.frame()

# Summarize across plants and years (person_summaries)
person_summaries <- ppw %>%
  group_by(person_id) %>%
  summarize(spp = n_distinct(nn_common_name),
            plants = n_distinct(individual_id),
            mn_dates_per_plant = round(mean(n_obs), 1),
            mn_interval = round(mean(interval_mn, na.rm = TRUE), 1),
            .groups = "keep") %>%
  mutate(mn_interval = replace_na(mn_interval, NA)) %>%
  data.frame()

# Join everything together:
persons <- recs %>%
  left_join(person_summaries, by = "person_id") %>%
  left_join(pw_wide, by = "person_id") %>%
  relocate(records_2024, .before = "spp_2024") %>%
  relocate(records_2025, .before = "spp_2025")

# Attach observer email/name --------------------------------------------------#
if (access_via_R) {

  # Load necessary packages
  library(DBI)
  library(odbc)
  library(dbplyr)

  con <- DBI::dbConnect(odbc::odbc(),
                        Driver   = "MySQL ODBC 9.4 ANSI Driver",
                        Server   = "usanpn-databases-prod.c0xzlo6s7duc.us-west-2.rds.amazonaws.com",
                        UID      = rstudioapi::askForPassword("Database user"),
                        PWD      = rstudioapi::askForPassword("Database password"),
                        database = "usanpn2",
                        Port     = 3306)

  person_info <- con %>%
    tbl("Person") %>% 
    collect()
  
  person_info <- person_info %>%
    # Remove NEON people
    filter(Person_ID > 0) %>%
    select(Person_ID, First_Name, Middle_Name, Last_Name, Create_Date, email) %>%
    rename(person_id = Person_ID,
           First = First_Name,
           Middle = Middle_Name,
           Last = Last_Name, 
           create_date = Create_Date) %>%
    data.frame()
  
  person_tx <- person_info %>%
    filter(person_id %in% persons$person_id)
  
  
  # Save to csv
  today <- str_remove_all(Sys.Date(), "-")
  write.csv(person_tx,
            paste0("data/ttr-tx-observers-", today, ".csv"),
            row.names = FALSE)

} else {
    
    person_tx <- read.csv("data/ttr-tx-observers-20251006.csv") 

}

persons <- person_tx %>%
  select(person_id, First, Last, email) %>%
  left_join(persons, by = "person_id")

# Write to file
write.csv(persons, 
          "data/ttr-tx-observer-summary-wy20242025.csv",
          row.names = FALSE)
