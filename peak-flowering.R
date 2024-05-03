################################################################################
# Estimating peak flowering

# Much of code is derived from that originally developed by Alyssa Rosemartin, 
# Hayley Limes, and Jeff Oliver. 
# See: https://github.com/alyssarosemartin/time-to-restore

# Erin Zylstra
# 2024-05-03
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
npn_orig <- npn_download_status_data(request_source = "erin",
                               years = 2013:2023,
                               species_ids = spp$species_id,
                               phenophase_ids= phenophases,
                               additional_fields = c("observedby_person_id"),
                               climate_data = FALSE)
npn_orig <- data.frame(npn_orig)

df <- npn_orig

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
df <- df %>%
  group_by(individual_id, year) %>%
  mutate(n_observers = n_distinct(observedby_person_id)) %>%
  ungroup() %>%
  data.frame()

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

# Check if any coordinate are problematic otherwise, using the maps package
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

# Put data into wide form -----------------------------------------------------#

# To make it easier to find data entry problems and calculate the number of open
# flowers when both intensity values were provided, we'll create a new dataframe
# where each row contains all information an observer collected about an
# individual plant on a given day

# Convert intensity value ranges to midpoints so they're easier to work with
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

# Identify when observers made multiple observations of a plant phenophase 
# in one day
df <- df %>%
  group_by(observedby_person_id, individual_id, observation_date) %>%
  mutate(n_records = n(),
         n_fl = sum(phenophase_id == 500),
         n_fo = sum(phenophase_id == 501)) %>%
  ungroup() %>%
  data.frame()
count(df, n_records, n_fl, n_fo)
  # Vast majority have 2 records (1 flower, 1 open flower obs), as expected
  # But there are 29 instances where an observer made 2 observations of a 
  # phenophase in a day and >2000 instances with an open flower obs and no 
  # corresponding flower obs

# When an observer appears to have made two observations of a plant phenophase 
# in one day, will keep the observation with more advanced phenophase, higher 
# intensity value, or more information. Will do this by sorting observations and 
# keeping only the first
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

# Remove these extra observations and unnecessary columns
df <- df %>%
  filter(dups == 1) %>%
  select(-c(obsnum, dups, n_records, n_fl, n_fo)) %>%
  arrange(species_id, observation_date, observedby_person_id, phenophase_id)

####################### NEED TO PICK UP HERE ###################################
# I think we need to do something here to deal with multiple phenophase 
# observations of the same plant on the same date by different observers.

# SEE NOTES




# Convert status information into wide form, with each row containing all 
# data recorded by an observer of a given plant on a given day
pstatus <- df %>%
  select(-c(observation_id, phenophase_id, intensity_value, midpoint)) %>%
  pivot_wider(names_from = phenophase_description,
              values_from = phenophase_status) %>%
  data.frame() %>%
  rename(flowers_s = Flowers.or.flower.buds,
         flowers_open_s = Open.flowers)

# Convert intensity information into wide form, with each row containing all 
# data recorded by an observer of a given plant on a given day
pint <- df %>%
  select(-c(observation_id, phenophase_id, phenophase_status, intensity_value)) %>%
  pivot_wider(names_from = phenophase_description,
              values_from = midpoint) %>%
  data.frame() %>%
  rename(flowers_i = Flowers.or.flower.buds,
         flowers_open_i = Open.flowers)

# Merge status and intensity dataframes
obs <- left_join(pstatus, pint)

# Address any inconsistencies -------------------------------------------------#

# If there was no flowers status obs (= NA), but open flowers = 1, change 
# flowers status to 1
obs$flowers_s[is.na(obs$flowers_s) & obs$flowers_open_s == 1] <- 1

# If observers said no to flowers (0) but yes to open flowers (1):
  # Make flower status = 1 if there's an intensity value for open flowers or if
  # the observation occurred between day 100 and 300. Delete other observations
  # because the open flower status is questionable.
  
  # Check before running this:
    count(obs, flowers_s, flowers_open_s, is.na(flowers_open_i),
          day_of_year %in% 100:300)
  # Looks like they'll be just two observations left where the open flower 
  # status is questionable.
  
obs <- obs %>%
  mutate(flowers_s = case_when(
    flowers_s == 0 & flowers_open_s == 1 & 
      (!is.na(flowers_open_i) | day_of_year %in% 100:300) ~ 1,
    .default = flowers_s)
  )
obs <- obs %>%
  filter(!(!is.na(flowers_s) & !is.na(flowers_open_s) & 
             flowers_s == 0 & flowers_open_s == 1))

# Looks like no observations with flowers_s = NA and flowers_open_s = 0 have
# intensity values for open flowers, so we can leave as is.

# Need to identify where we have multiple observations of the same plant on the
# same date (by different observers, since we've already taken care of multiple
# observations by the same observer). Where status is 1, we could calculate the
# mean of intensity values. However, I think it's a problem if different 
# observers reported different phenophase statuses for the same plant on the 
# same date, whether we're characterizing peak flowering or open flower 
# phenophases.  

# Create indicators for entries that have intensity values
obs <- obs %>%
  arrange(species_id, individual_id, observation_date, 
          desc(flowers_s), desc(flowers_open_s), desc(flowers_open_i)) %>%
  mutate(fl_i_bin = if_else(is.na(flowers_i), 0, 1),
         fo_i_bin = if_else(is.na(flowers_open_i), 0, 1),
         both_i = fl_i_bin * fo_i_bin) %>%
  group_by(individual_id, observation_date) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  data.frame()

# Look at these extreme examples to highlight variation in observations! Does
# it even make sense to average intensity values? And if we do take averages, 
# are we averaging over set of flower intensity values and then averaging over 
# set of open flower intensity values even if they came from different people?
# If we're going to calculate the nubmer of open flowers, then should we do
# that first for any observation that had both intensity values, then average?
filter(obs, n_obs == 10)
filter(obs, n_obs == 9)
filter(obs, n_obs == 6)
filter(obs, n_obs == 5)



plantdate <- obs %>%
  group_by(individual_id, observation_date) %>%
  filter(n() > 1) %>%
  summarize(fl_s_min = if_else(all(is.na(flowers_s)), NA, 
                               min(flowers_s, na.rm = TRUE)),
            fl_s_max = if_else(all(is.na(flowers_s)), NA, 
                               max(flowers_s, na.rm = TRUE)),
            .groups = "keep") %>%
  data.frame()



            fl_s_min = if_else(all(is.na(flowers_s)), NA, 
                               min(flowers_s, na.rm = TRUE)),
            fl_s_max = if_else(all(is.na(flowers_s)), NA, 
                               max(flowers_s, na.rm = TRUE)),
            fo_s_min = if_else(all(is.na(flowers_open_s)), NA, 
                               min(flowers_s, na.rm = TRUE)),
            fo_s_max = if_else(all(is.na(flowers_open_s)), NA, 
                               max(flowers_s, na.rm = TRUE)),
            .groups = "keep") %>%
  data.frame() %>%
  mutate(flowers_same = if_else(fl_s_min == fl_s_max, 1, 0),
         flowers_open_same = if_else(fo_s_min == fo_s_max, 1, 0))

  
  
  
  # Calculate number of observations
  mutate(n_obs = n()) %>%
  # Create indicator to see whether flower/open flower status is consistent
  
  
  mutate(flowers_s_same = ifelse(length(unique(flowers_s)) == 1, 1, 0),
         flowers_open_s_same = ifelse(length(unique(flowers_open_s)) == 1, 1, 0)) %>%
  ungroup() %>%
  data.frame()




