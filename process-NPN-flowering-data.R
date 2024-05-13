################################################################################
# Download and process NPN flowering data (status and intensities)

# Downloads and cleans NPN flowering status data for priority species,
# summarizes data for each plant and date when any observations were made.  
# Generates an estimated number of open flowers when one or more observers 
# recorded the number of flowers and one or more observers recorded the 
# proportion of flowers that were open. 

# Much of code is derived from that originally developed by Alyssa Rosemartin, 
# Hayley Limes, and Jeff Oliver. 
# See: https://github.com/alyssarosemartin/time-to-restore

# Erin Zylstra
# 2024-05-06
################################################################################

require(rnpn)
require(dplyr)
require(lubridate)
require(tidyr)
require(maps)

rm(list = ls())

# Download and format NPN data for priority plant species ---------------------#

# Load csv with priority plant species
spp <- read.csv("data/priority-species.csv")
  # species_id: NPN species ID
  # LA, OK, NM, TX: 1/0 indicating whether the species got one or more votes in 
    # that state based on googlesheets (note that TX votes were only included
    # in the most recent "revisited" sheet)
  # priority: priority levels (1 = highest priority, as determined by NPN-TTR
    # team; 2 = includes species that received votes from multiple states; 3 =
    # species included in 16-priority list used previously)

# Load information (from other sources) about flowering season
spp_flower <- read.csv("data/priority-species-flowering-season.csv")
  # first_flower_ttr: first flowering month based on information on priority
    # species in TTR google sheets (not sure of original source)
  # first_flower_UT: first flowering month based on U. Texas - Lady Bird
    # Johnson wildflower center (https://www.wildflower.org)
  # first_flower_season: current classification of each species, as spring or 
    # summer flowering. There are a few that still need to be classified (noted
    # as spring/summer).

# Merge the two:
spp <- left_join(spp, spp_flower, by = c("common_name", "genus", "species"))

# Remove species that don't have an NPN species_id and format flowering columns
spp <- spp %>%
  filter(!is.na(species_id)) %>%
  mutate(first_flower_ttr = match(first_flower_ttr, month.abb),
         first_flower_UT = match(first_flower_UT, month.abb))  

# Identify phenophases of interest
phenophases <- c(500, 501)
  # 500 = Flowers or flower buds
  # 501 = Open flowers

# Download NPN data
df <- npn_download_status_data(request_source = "erin",
                               years = 2013:2023,
                               species_ids = spp$species_id,
                               phenophase_ids= phenophases,
                               additional_fields = c("observedby_person_id"),
                               climate_data = FALSE)
df <- data.frame(df)

# Quick look at relative amount of data for each species
  # df %>%
  #   count(species_id, genus, species, common_name) %>%
  #   arrange(desc(n))
# Note that Eryngium yuccifolium is listed as "rattlesnake master" in priority
# lists, but "button eryngo" in NPN database. And Cirsium undulatum (wavy-leaved
# thistle) doesn't have any records.

# Extract year from date and remove a few columns we don't need
df <- df %>%
  mutate(year = lubridate::year(observation_date)) %>%
  select(-c(update_datetime, kingdom, intensity_category_id, abundance_value))

# Count number of people who contributed observations of a given plant in a year
# df <- df %>%
#   group_by(individual_id, year) %>%
#   mutate(n_observers = n_distinct(observedby_person_id)) %>%
#   ungroup() %>%
#   data.frame()

# Remove rows that are duplicates (ie, same plant, date, phenophase status, 
# intensity_value), regardless of who reported the data (observedby_person_id)
df <- df %>% 
  distinct(individual_id, observation_date, phenophase_description, 
           phenophase_status, intensity_value, .keep_all = TRUE)

# Clean up phenophase status and state data -----------------------------------#

# Remove any observations with phenophase_status = -1 (unknown), unless there's
# an intensity value, which suggest the phenophase_status was just mis-entered
df <- df %>%
  mutate(intensity_value = ifelse(intensity_value == "-9999", 
                                  NA, intensity_value)) %>%
  mutate(phenophase_status = ifelse(phenophase_status == -1 & 
                                      !is.na(intensity_value),
                                    1, phenophase_status)) %>%
  filter(phenophase_status != -1)

# Remove a couple observations with errant state names
df <- filter(df, state != "Liaoning Sheng")

# Look at geographic coordinates when state = "-9999"
df %>% filter(state == "-9999") %>% select(latitude, longitude) %>% summary()
df %>% filter(state == "-9999") %>% count(common_name) 
  # All observations of golden crownbeard on Midway Atoll. 
  # Deleting since phenology could be very different there.
df <- filter(df, state != "-9999")

# Check if any coordinates are problematic otherwise, using the maps package
# to identify state by lat/lon
df <- df %>%
  mutate(state_name = maps::map.where("state", df$longitude, df$latitude)) %>%
  left_join(select(maps::state.fips, polyname, abb),
            by = c("state_name" = "polyname"))
df %>%
  filter(!is.na(abb) & state != abb) %>%
  count(state, abb)
  # There are some ND/MN and OH/KY differences, but that's probably not a big
  # deal since states are adjacent. However, will remove a couple entries with
  # NPN database listing NM, but lat/lon location is in SC
df <- df %>% 
  filter(!(state == "NM" & abb == "SC")) %>%
  select(-c(state_name, abb))

# Deal with observations of the same plant, phenophase on the same date -------#

# First, converting intensity value ranges to midpoints so they're easier to 
# work with
df <- df %>%
  mutate(midpoint = case_when(
    intensity_value == "Less than 3" ~ 2, 
    intensity_value == "3 to 10" ~ 7, 
    intensity_value == "11 to 100" ~ 56, 
    intensity_value == "101 to 1,000" ~ 551, 
    intensity_value == "1,001 to 10,000" ~ 5510, 
    intensity_value == "More than 1,000" ~ 1001, 
    intensity_value == "More than 10,000" ~ 10001, 
    intensity_value == "Less than 5%" ~ 0.05, 
    intensity_value == "5-24%" ~ 0.15, 
    intensity_value == "25-49%" ~ 0.37, 
    intensity_value == "50-74%" ~ 0.62, 
    intensity_value == "75-94%" ~ 0.85, 
    intensity_value == "95% or more" ~ 0.95,
    .default = NA)
  )
# Check that status equals 1 if there's an intensity value
count(df, phenophase_status, midpoint)

# Identify when an observer made multiple observations of a plant phenophase in
# one day
df <- df %>%
  group_by(observedby_person_id, individual_id, observation_date) %>%
  mutate(n_records = n(),
         n_fl = sum(phenophase_id == 500),
         n_fo = sum(phenophase_id == 501)) %>%
  ungroup() %>%
  data.frame()
# count(df, n_records, n_fl, n_fo)

# When an observer appears to have made two observations of a plant phenophase 
# in one day, will keep the observation with more advanced phenophase, higher 
# intensity value, or more information. Will do this by sorting observations and 
# keeping only the first.
obsdatep <- df %>%
  group_by(observedby_person_id, individual_id, observation_date, phenophase_id) %>%
  summarize(n_obs = n(),
            .groups = "keep") %>%
  data.frame()
obsdatep$obsnum <- 1:nrow(obsdatep)

df <- df %>%
  arrange(observedby_person_id, individual_id, observation_date, phenophase_id,
          desc(phenophase_status), desc(midpoint)) %>%
  left_join(select(obsdatep, -n_obs), 
            by = c("observedby_person_id", "individual_id", "observation_date", 
                   "phenophase_id")) %>%
  # Create "dups" column, where dups > 1 indicates that the observation can be
  # removed since there's another observation that same day with more advanced
  # phenology or more information.
  mutate(dups = sequence(rle(as.character(obsnum))$lengths))

# Remove extra observations and unnecessary columns
df <- df %>%
  filter(dups == 1) %>%
  select(-c(obsnum, dups, n_records, n_fl, n_fo)) %>%
  arrange(species_id, observation_date, observedby_person_id, phenophase_id)

# Deal with multiple observations on the same date of the same plant by
# different observers (~0.3% of plant-phenophase-date combinations)
datep <- df %>%
  group_by(individual_id, observation_date, phenophase_id) %>%
  summarize(n_obs = n(),
            .groups = "keep") %>%
  data.frame()
datep$obsnum <- 1:nrow(datep)

df <- df %>%
  arrange(individual_id, observation_date, phenophase_id, 
          desc(phenophase_status), desc(midpoint)) %>%
  left_join(select(datep, -n_obs), 
            by = c("individual_id", "observation_date", "phenophase_id")) %>%
  mutate(dups = sequence(rle(as.character(obsnum))$lengths))

# We can:
# 1) just select one observation with the most advanced phenology or most info
# 2) use status = 1 and average over midpoints
# 3) identify observations when both intensity values were reported and then 
#    only use those observations to calculate number of open flowers (and 
#    average across them if there were more than one). Note that it would be 
#    easiest to have data in wide form to do this.
# For now, we'll go with option 2 because it better mirrors the approach that 
# others used previously. 

# Remove "duplicates" with status == 0
df <- df %>% filter(!(dups > 1 & phenophase_status == 0))

# Recreate obsnum column now that we've removed observations
datep <- df %>%
  group_by(individual_id, observation_date, phenophase_id) %>%
  summarize(n_obs = n(),
            .groups = "keep") %>%
  data.frame()
datep$obsnum <- 1:nrow(datep)

df <- df %>%
  select(-c(obsnum, dups)) %>%
  left_join(datep, by = c("individual_id", "observation_date", "phenophase_id"))

# Calculate mean of intensity midpoints, where reported
mean_midpoints <- df %>%
  filter(n_obs > 1) %>%
  group_by(obsnum) %>%
  summarize(midpoint_mn = if_else(all(is.na(midpoint)), NA, 
                                  mean(midpoint, na.rm = TRUE))) %>%
  data.frame()
df <- df %>% 
  left_join(mean_midpoints, by = "obsnum") %>%
  mutate(midpoint = if_else(is.na(midpoint_mn), midpoint, midpoint_mn)) %>%
  select(-c(midpoint_mn, observation_id, observedby_person_id, 
            intensity_value, obsnum)) %>%
  distinct()
# Keeping n_obs column to remember how many observers recorded info about that
# plant and phenophase on that date, but removing the obsnum column

# Now we're left with a single "observation" (status, mean intensity value if
# at least one observer reported it) for each plant, phenophase, and year

# Aggregate data for each plant-date combination ------------------------------#

# To make it easier to find inconsistencies in the data and to calculate the 
# number of open flowers when both intensity values were provided, we'll create 
# a new dataframe where each row contains all information about an individual 
# plant on a given day

# Convert status information into wide form (fl_s = status of flowers 
# phenophase; fo_s = status of open flowers phenophase)
pstatus <- df %>%
  select(-c(phenophase_id, midpoint, n_obs)) %>%
  pivot_wider(names_from = phenophase_description,
              values_from = phenophase_status) %>%
  data.frame() %>%
  rename(fl_s = Flowers.or.flower.buds,
         fo_s = Open.flowers)

# Convert intensity information into wide form (fl_i = number of flowers
# midpoint; fo_i = proportion of flowers open midpoint)
pint <- df %>%
  select(-c(phenophase_id, phenophase_status, n_obs)) %>%
  pivot_wider(names_from = phenophase_description,
              values_from = midpoint) %>%
  data.frame() %>%
  rename(fl_i = Flowers.or.flower.buds,
         fo_i = Open.flowers)

# Convert number of observers for each plant-phenophase-date into wide form 
# (n_obs_fl = number of observers of flower phenophase; n_obs_fo = number of
# observers of open flower phenophase)
pobs <- df %>%
  select(-c(phenophase_id, phenophase_status, midpoint)) %>%
  pivot_wider(names_from = phenophase_description,
              values_from = n_obs) %>%
  data.frame() %>%
  rename(n_obs_fl = Flowers.or.flower.buds,
         n_obs_fo = Open.flowers)

# Merge status, intensity, observer dataframes
obs <- left_join(pstatus, pint) %>%
  left_join(pobs)

# Address inconsistencies -----------------------------------------------------#

# If there was no flowers status data (= NA), but open flowers = 1, change 
# flowers status to 1
obs$fl_s[is.na(obs$fl_s) & obs$fo_s == 1] <- 1

# If observers said no to flowers (0) but yes to open flowers (1):
  # Make flower status = 1 if there's an intensity value for open flowers or if
  # the observation occurred between day 100 and 300. Delete other observations
  # because the open flower status is questionable.
  
  # Check before running this:
    count(obs, fl_s, fo_s, is.na(fo_i), day_of_year %in% 100:300)
    # Very few observations left where the open flower status is questionable
  
obs <- obs %>%
  mutate(fl_s = case_when(
    fl_s == 0 & fo_s == 1 & (!is.na(fo_i) | day_of_year %in% 100:300) ~ 1,
    .default = fl_s)
  )
obs <- obs %>%
  filter(!(!is.na(fl_s) & !is.na(fo_s) & fl_s == 0 & fo_s == 1))

# Looks like no observations with fl_s = NA and fo_s = 0 have intensity values 
# for open flowers, so we can leave as is.

# Calculate estimated number of open flowers ----------------------------------#

# For each observation where observers confirmed that flowers and open flowers
# were present and where intensity values were provided, we can calculate the
# number of open flowers as the product of the number of flowers (midpoint) and 
# proportion of flowers open (midpoint)
obs <- obs %>%
  mutate(num_open = if_else(is.na(fl_i) | is.na(fo_i), NA, fl_i * fo_i))

# Write this dataset to file --------------------------------------------------#
write.csv(obs, 
          "data/flower-status-intensities-priorityspp.csv", 
          row.names = FALSE)
