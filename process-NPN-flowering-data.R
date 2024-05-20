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
# 2024-05-20
################################################################################

require(rnpn)
require(dplyr)
require(lubridate)
require(tidyr)
require(maps)

# Download and format NPN data for priority plant species ---------------------#

# Load csv with priority plant species (only those with specific epithets)
spp <- read.csv("data/priority-species.csv")
  # species_id: NPN species ID
  # votes_ST: number of votes each species received from each state (LA, NM, OK) 
    # in original working groups
  # revisit_votes_ST: number of votes each species received from each state in 
    # "revisited" spreadsheet
  # rvote (1/0): indicates whether species got >=1 vote in revisited spreadsheet
  # nvotes_orig: sum of votes for each species in original working groups
  # prev16 (1/0): indicator whether species was part of previous 16-spp group
  # priority (1-4): indicator of priority class:
    # 1: original 8 priority spp
    # 2: received one or more vote in revisited spreadsheet
    # 3: in original 16 priority species
    # 4: received >1 vote (in one or across multiple states) in orignal spreadsheets
  # species_id: NPN database species ID. NA if species not in database
  # common_name: common name in NPN database
  # common_name_states: common name in spreadsheets

# Evaluate how many priority species do/don't appear in NPN database
count(spp, is.na(species_id), priority)
  # 11/52 species are not in database (priority 2 and 4)
  # Left with 41 spp (8 priority 1; 13 priority 2; 6 priority 3; 14 priority 4)

# Remove species that don't have an NPN species_id
spp <- spp %>% 
  filter(!is.na(species_id))

# Identify phenophases of interest
phenophases <- c(500, 501)
  # 500 = Flowers or flower buds
  # 501 = Open flowers

# Download NPN data
df <- npn_download_status_data(request_source = "erin",
                               years = 2013:2023,
                               species_ids = spp$species_id,
                               phenophase_ids= phenophases,
                               additional_fields = c("observedby_person_id",
                                                     "partner_group",
                                                     "site_id",
                                                     "site_name"),
                               climate_data = FALSE)
df <- data.frame(df)

# Quick look at relative amount of data for each species
  df %>%
    count(species_id, genus, species, common_name) %>%
    arrange(desc(n))
  
  dplyr::setdiff(spp$common_name, unique(df$common_name))
  # No records for:
  # wavyleaf thistle (Cirsium undulatum)
  # swamp sunflower (Helianthus angustifolius)
  # lyreleaf sage (Salvia lyrata)

# Extract year from date and remove a few columns we don't need
df <- df %>%
  mutate(year = lubridate::year(observation_date)) %>%
  select(-c(update_datetime, kingdom, intensity_category_id, abundance_value))

# Clean up phenophase status and state data -----------------------------------#

# Remove any observations with phenophase_status = -1 (unknown) unless there's
# an intensity value, which suggest the phenophase_status was just mis-entered
df <- df %>%
  mutate(intensity_value = ifelse(intensity_value == "-9999", 
                                  NA, intensity_value)) %>%
  mutate(phenophase_status = ifelse(phenophase_status == -1 & 
                                      !is.na(intensity_value),
                                    1, phenophase_status)) %>%
  filter(phenophase_status != -1)

# Remove observations that are outside the US based on "state" entry (but 
# leaving missing values for the moment). Removing HI and AK observations since
# phenology is likely to be very different in those places.
df <- df %>%
  filter(state %in% c(state.abb, "-9999")) %>%
  filter(!state %in% c("HI", "AK"))

# Look at geographic coordinates when state = "-9999"
df %>% filter(state == "-9999") %>% select(latitude, longitude) %>% distinct()
df %>% filter(state == "-9999") %>% count(common_name) 
  # Just 5 locations, 4 of which are of golden crownbeard on Midway Atoll
  # Other is in India (might be FL typo, but doesn't really matter)
  # Deleting since phenology could be very different there.
df <- filter(df, state != "-9999")

# Creating bounding box for lower 48 states and remove any locations beyond that
# even if a state was provided
extreme_lats <- c(24, 50)
extreme_lons <- c(-125, -66)
df <- df %>%
  filter(latitude > extreme_lats[1] & latitude < extreme_lats[2]) %>%
  filter(longitude > extreme_lons[1] & longitude < extreme_lons[2])

# Check if any coordinates are problematic otherwise, using the maps package
# to identify state by lat/lon
df <- df %>%
  mutate(state_name = maps::map.where("state", df$longitude, df$latitude)) %>%
  left_join(select(maps::state.fips, polyname, abb),
            by = c("state_name" = "polyname"))
df %>%
  filter(!is.na(abb) & state != abb) %>%
  count(state, abb)
  # There are instances where the coordinates fall in different states (based on
  # the maps package) than the state that was listed in the NPN database, but 
  # most of the time, it's probably not a big deal since the mismatched states 
  # are geographically adjacent.  However, there are a few instances where 
  # either the state listed or coords must be wrong, so will remove those obs.
df <- df %>% 
  filter(!(!is.na(abb) & state == "NM" & abb == "SC")) %>%
  filter(!(!is.na(abb) & state == "MN" & abb == "CA")) %>%
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

# When an observer appears to have made multiple observations of a plant 
# phenophase in one day, will keep the observation with more advanced 
# phenophase, higher intensity value, or more information. Will do this by 
# sorting observations and keeping only the first.
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
# different observers
datep <- df %>%
  group_by(individual_id, site_id, site_name, partner_group, 
           observation_date, phenophase_id) %>%
  summarize(n_obs = n(),
            .groups = "keep") %>%
  data.frame()
datep$obsnum <- 1:nrow(datep)

# Quantify the extent to which this is occurring:
sum(datep$n_obs > 1) / nrow(datep) * 100
  # 6.8% of plant-date-phenophase combinations are associated with observations 
  # from multiple people. 

# Look at some extreme examples where we have >5 observations of a plant 
# phenophase in a day
datep_high <- filter(datep, n_obs > 5)
count(datep_high, partner_group, site_name)
  # Lots of site_name = Oak Hill Campus OR Oak Hill High School 2022
  oh <- filter(df, grepl("Oak Hill", site_name))
  dim(oh) # 2307 plant phenophase observations
  count(oh, individual_id, common_name, year)
  # Lots of partner_group = UNCO BIO 111 & site_name = Red Maple site
  unco <- filter(df, grepl("UNCO", partner_group))
  dim(unco) # 27730 plant phenophase observations
  count(unco, individual_id, common_name, year)
  # Other groups/sites with many instances of multiple observations per day:
  # New York Botanical Garden Forest Phenology (group)
  # Mohonk Preserve (group)
# Conferred with Erin P on this. She suggested removing data from Oak Hill and 
# UNCO BIO (which they've had other issues with). NY Botanical Garden and Mohonk 
# Preserve are reliable groups so keep their data in.
df <- df %>%
  filter(!grepl("Oak Hill", site_name)) %>%
  filter(!grepl("UNCO", partner_group))

# Recreate summary dataframe and recaluate the extent to which we have multiple
# observations of same plant on same date
datep <- df %>%
  group_by(individual_id, site_id, site_name, partner_group, 
           observation_date, phenophase_id) %>%
  summarize(n_obs = n(),
            .groups = "keep") %>%
  data.frame()
datep$obsnum <- 1:nrow(datep)
sum(datep$n_obs > 1) / nrow(datep) * 100
# Now, 5.3% of plant-date-phenophase combinations are associated with observations 
# from multiple people. 

# Isolate plant-year combinations where one or more times, a phenophase was 
# observed by multiple people on the same day.
yearp <- datep %>%
  mutate(year = year(observation_date)) %>%
  group_by(individual_id, year) %>%
  summarize(multobs = 1 * any(n_obs > 1), .groups = "keep") %>%
  data.frame()
yearp_mult <- filter(yearp, multobs == 1)

# For all observers who collected data on a plant in a given year, rank them
# by the number of observations made.
for (i in 1:nrow(yearp_mult)) {
  dfsub <- df %>%
    filter(individual_id == yearp_mult$individual_id[i]) %>%
    filter(year == yearp_mult$year[i])
  obs_count <- count(dfsub, observedby_person_id) %>%
    arrange(desc(n)) %>%
    mutate(observer_rank = row_number()) %>%
    mutate(individual_id = yearp_mult$individual_id[i], 
           year = yearp_mult$year[i]) 
  if (i == 1) {
    obs_rank <- obs_count
  } else {
    obs_rank <- rbind(obs_rank, obs_count)
  }
}
# Attach information about duplicate observations to the original dataframe
df <- df %>%
  left_join(select(datep, -c(site_id, site_name, partner_group)),
            by = c("individual_id", "observation_date", "phenophase_id")) %>%
  left_join(select(obs_rank, -n), 
            by = c("individual_id", "year", "observedby_person_id"))
# Check that there's always an observer rank if there are multiple observations
# of that plant phenophase on that date:
summary(df$observer_rank[df$n_obs > 1])
summary(df$n_obs[is.na(df$observer_rank)])

# Now we have a couple options:
# 1) retain one observation of a plant phenophase per day, preferentially
#    selecting those by people who submitted consistent observations.
# 2) retain all data, use status = 1 observations, and average over midpoints
# For now will use Option 1 even though Option 2 is most similar to the approach 
# that was used previously

# Remove all but one observation of a plant phenophase per day
df <- df %>%
  # Replace NAs in observer rank with 0 for easier sorting
  mutate(observer_rank = replace_na(observer_rank, 0)) %>%
  # Sort dataframe by plant-phenophase-date and then by observer rank
  arrange(obsnum, observer_rank) %>%
  # Remove all but one observation and remove unnecessary columns
  mutate(dups = sequence(rle(as.character(obsnum))$lengths)) %>%
  filter(dups == 1) %>%
  select(-c(n_obs, obsnum, observer_rank, dups))

# Now we're left with a single observation (status, intensity value if
# observer reported it) for each plant, phenophase, and date

# Aggregate data for each plant-date combination ------------------------------#

# To make it easier to find inconsistencies in the data and to calculate the 
# number of open flowers when both intensity values were provided, we'll create 
# a new dataframe where each row contains all information about an individual 
# plant on a given day

# First remove columns that won't be needed moving forward
df <- df %>%
  select(-c(observation_id, observedby_person_id, 
            phenophase_id, intensity_value))

# Convert status information into wide form (fl_s = status of flowers 
# phenophase; fo_s = status of open flowers phenophase)
pstatus <- df %>%
  select(-midpoint) %>%
  pivot_wider(names_from = phenophase_description,
              values_from = phenophase_status) %>%
  data.frame() %>%
  rename(fl_s = Flowers.or.flower.buds,
         fo_s = Open.flowers)

# Convert intensity information into wide form (fl_i = number of flowers
# midpoint; fo_i = proportion of flowers open midpoint)
pint <- df %>%
  select(-phenophase_status) %>%
  pivot_wider(names_from = phenophase_description,
              values_from = midpoint) %>%
  data.frame() %>%
  rename(fl_i = Flowers.or.flower.buds,
         fo_i = Open.flowers)

# Merge status, intensity, observer dataframes
obs <- left_join(pstatus, pint)

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

count(obs, fl_s, fo_s, is.na(fl_i), is.na(fo_i))
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
