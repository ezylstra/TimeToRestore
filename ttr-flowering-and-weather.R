# Exploring options for models of weather effects on flowering onset
# 27 July 2026

require(lubridate)
require(stringr)
require(dplyr)
require(tidyr)
require(terra)
# require(sf)
require(ggplot2)
library(RColorBrewer)


# Load functions --------------------------------------------------------------#

functions <- list.files("functions", full.names = TRUE)
for (f in functions) {
  source(f)
}

# Load/format observation data ------------------------------------------------#

df <- read.csv("data/status-intensity-flowers-May2026.csv")
df <- df %>%
  mutate(obsdate = ymd(obsdate))

# This includes data from 2013 - May 2026 and 12 states (SC region plus
# surrounding states). Includes 51 species: 6 that have <100 observations, 21 
# species with >1000 observations. Top 3 species have way more observations 
# than others: red maple (50974), eastern redbud (12281), beautyberry (9999).

# Update series function
# Create series
# Use series dataset to find first yeses with a prior no (within 14 days?)
# May want to use different period for red maple (or others species?) -- want
  # to evaluate this in some systematic way...

# Identify what constitutes "sufficient" data for a species to be included in
# subsequent analyses
# Need a minimum number of unique combos: year, lat/lon (need the lat/lon to 
# have different weather data though, so need a minimum distance between
# locations). Could attach cell ID from a PRISM raster to identify unique
# "weather locs".

# What climate/weather data are we going to use? PRISM?
# I didn't download the status observations with daymet data, and even if I had
# done that we wouldn't have daymet data for most recent observations (or for
# future observations moving forward). I do have daily PRISM data on my computer
# from 2013-2025, so would just need to download data for 2026 (through June?)
# and 2012 data (or just use status data from 2014-2026). 

# What climate variables are we going to calculate?
# Fixed seasonal variables using 3-month NPN delineations (eg, mean spring temp)
# Weather summaries immediately preceding event (or average date of event)?
# For a given species, could find average date in SC region that flowering is
# first observed and then use 3-months prior (eg, if average first observation
# date is 14 April, then use mean temperatures for Jan-Mar)


# Limit to species that have at least 100 observations (bare minimum!)
df <- df %>% 
  group_by(spp) %>%
  mutate(nobs = n()) %>%
  filter(nobs > 100) %>%
  ungroup() %>%
  select(-nobs) %>%
  data.frame()

# If there's more than one observation of a plant on a given date, select one 
# with positive status and (higher) intensity value
df <- df %>%
  arrange(id, obsdate, php_id, desc(status), desc(midpoint)) %>%
  distinct(id, obsdate, php_id, .keep_all = TRUE)

