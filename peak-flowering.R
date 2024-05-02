################################################################################
# Estimating peak flowering

# Much of code is derived from that originally developed by Alyssa Rosemartin, 
# Hayley Limes, and Jeff Oliver. 
# See: https://github.com/alyssarosemartin/time-to-restore

# Erin Zylstra
# 2024-05-02
################################################################################

require(rnpn)
require(dplyr)
require(lubridate)
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




# Check if observers report No to Flowers or Flower buds but Yes to Open Flowers
# First, see if observers made multiple observations of a plant in a day
inddate <- df %>%
  group_by(observedby_person_id, individual_id, observation_date) %>%
  summarize(n_records = n(),
            n_fl = sum(phenophase_id == 500),
            n_fo = sum(phenophase_id == 501),
            .groups = "keep") %>%
  data.frame()
count(inddate, n_records, n_fl, n_fo)
  # Vast majority have 2 records (1 flower, 1 flower open), as expected
  # 1246 have only the flower record (no data for flower open)
  # 2047 have only the flower open records (no data for flower)
  # 29 have 2 records for flowers, flowers open, or both

  
  
  mutate(i_flower = ifelse(phenophase_id == 500 & phenophase_status == 1, 1, 0),
         i_open = ifelse(phenophase_id == 501 & phenophase_status == 1, 1, 0),
         conflict = ifelse(i_flower == 0 & i_open == 1, 1, 0))



# Not sure the grouping works because there are occasionally entries with same
# observer, individual, date, phenophase...
df <- df %>% 
  group_by(observedby_person_id, individual_id, observation_date) %>% 
  mutate(FlowersNo = ifelse(phenophase_status == 0 & phenophase_id == 500, 1, NA)) %>% 
  mutate(OpenFlowersYes= ifelse(phenophase_status == 1 & phenophase_id ==501, 1, NA)) %>% 
  mutate(FlowersError = ifelse(FlowersNo == 1 & OpenFlowersYes == 1, 1, NA)) %>%
  ungroup() %>%
  data.frame()

# No inconsistencies, so we'll drop these fields
df <- subset(df, select = -c(FlowersNo,OpenFlowersYes,FlowersError))



# AT SOME POINT, NEED TO CHECK FOR DUP OBSERVATIONS, WHERE INTENSITY VALUES
# OR STATUS VALUES DIFFER.




# Quick look at the number of individuals monitored, by species
df %>%
  group_by(common_name) %>%
  summarize(n_individs = n_distinct(individual_id)) %>%
  data.frame() %>%
  arrange(desc(n_individs))
# Species with notably few individuals  
# green antelopehorn, pinkscale blazing star (1)
# aquatic milkweed (6)
# broadleaf milkweed (7)
# All other species have >= 16 individuals





# Check if observers report No to Flowers or Flower buds but Yes to Open Flowers
# (I might have more thorough checks in FlowerForBats script...)

# Not sure the grouping works because there are occasionally entries with same
# observer, individual, date, phenophase...
df <- df %>% 
  group_by(observedby_person_id, individual_id, observation_date) %>% 
  mutate(FlowersNo = ifelse(phenophase_status == 0 & phenophase_id == 500, 1, NA)) %>% 
  mutate(OpenFlowersYes= ifelse(phenophase_status == 1 & phenophase_id ==501, 1, NA)) %>% 
  mutate(FlowersError = ifelse(FlowersNo == 1 & OpenFlowersYes == 1, 1, NA)) %>%
  ungroup() %>%
  data.frame()

# No inconsistencies, so we'll drop these fields
df <- subset(df, select = -c(FlowersNo,OpenFlowersYes,FlowersError))

# Identify repeated observations on the same date/individual plant/phenophase
df <- df %>% 
  group_by(individual_id, observation_date, phenophase_id) %>% 
  mutate(same_date_ind_pp = as.integer(n() > 1))  %>% 
  ungroup() %>%
  data.frame()

# View the number of records with the same date/ind/pp - presumed conflicts
same_date_ind_pp_summary <- df %>%
  count(phenophase_description, same_date_ind_pp) %>%
  group_by(phenophase_description) %>%
  mutate(freq = n / sum(n)) %>%
  data.frame()

# About 1% of records for Flowers/Flower buds and 0.7% for Open Flowers have the 
# same date/phenophase/individual plant
same_date_records <- df %>% 
  subset(same_date_ind_pp == 1)










# Code midpoints to values to facilitate calculating the number of flowers
df <- df %>%
  mutate(year = lubridate::year(observation_date)) %>%
  mutate(intensity_value = na_if(intensity_value, "-9999")) %>%
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

# Get a count of observers who contributed to observing a given ind plant in a year
# Added ungroup and data.frame options since I don't want a grouped tbl
df <- df %>%
  group_by(individual_id, year) %>%
  mutate(number_observers = n_distinct(observedby_person_id)) %>%
  ungroup() %>%
  data.frame()

# Remove unneeded columns (I think I've figured it out based on their indices)
df <- select(df, -c(update_datetime, genus, species, kingdom, abundance_value))

# Remove rows that are duplicates (ie, same plant, date, phenophase status, 
# intensity_value), regardless of who reported the data (observedby_person_id)
df <- df %>% 
  distinct(individual_id, observation_date, phenophase_description, 
           phenophase_status, intensity_value, .keep_all = TRUE)

# I checked that each individual plant is only associated with a single site


