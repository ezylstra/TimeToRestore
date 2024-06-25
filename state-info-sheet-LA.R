################################################################################
# Create material for state info sheet: LA

# Much of code is derived from that originally developed by Alyssa Rosemartin, 
# Hayley Limes, and Jeff Oliver. 
# See: https://github.com/alyssarosemartin/time-to-restore

# Erin Zylstra
# 2024-06-24
################################################################################

require(dplyr)
require(lubridate)
require(stringr)
require(tidyr)
require(ggplot2)
library(terra)
library(tidyterra)
library(mgcv)
library(gridExtra)

rm(list = ls())

# Set figure parameters -------------------------------------------------------#
alphaline <- 0.3
alphapoly <- 0.4
figw <- 6.5
figh <- 3

# Load csv with sample sizes for weekly proportion of plants with open flowers #
ofss <- read.csv("data/openflower_weeklyobs_samplesizes.csv")

# Identify priority species with sufficient data ------------------------------#
species <- ofss %>%
  filter(LA == 1) %>%
  filter(mn_wkobs_sc >= 6) %>%
  select(!contains("z789"))
  # Note: no restrictions based on number of plant-years

# Load iNat data --------------------------------------------------------------#
# Save for later. May need to do a new download. Want to grab flowering tag if
# possible.

# Load processed NPN status/intensity data and format -------------------------#
df <- read.csv("data/flower-status-intensities-priorityspp.csv")

# Rename/remove columns where necessary
df <- df %>%
  rename(lat = latitude,
         lon = longitude,
         plant_id = individual_id) %>%
  select(-c(person_id, n_observations))

# Summarize data by plant-year
samples <- df %>%
  group_by(common_name, plant_id, site_id, state, year) %>%
  summarize(n_obs = n(),                      # No. of daily observations
            n_status_fl = sum(!is.na(status_fl)),  # No. of flower status obs
            n_status_fo = sum(!is.na(status_fo)),  # No. of open flower status obs
            n_value_fl = sum(!is.na(midpoint_fl)), # No. of flower intensity values
            n_value_fo = sum(!is.na(midpoint_fo)), # No. of open flower intensity values
            .groups = "keep") %>%
  data.frame()

# Remove any plant year with no open flower status recorded
samples <- samples %>% 
  filter(n_status_fo > 0) %>%
  mutate(plant_yr = paste(plant_id, year, sep = "_"))
df <- df %>%
  mutate(plant_yr = paste(plant_id, year, sep = "_")) %>%
  filter(plant_yr %in% samples$plant_yr)

# Add week to data, so we can calculate weekly proportions
  # We'll be creating wk_doy column to assign each week with a day of the year.
  # If we want this date to be the start of each week (eg, date for week 1 would 
  # be Jan 1), then set day_of_week = 1. If we want the date to be mid week 
  # (Thursday) then set day_of_week = 4. 
  day_of_week <- 1
  df <- df %>%
    mutate(wk = week(observation_date),
           wk_date = parse_date_time(paste(2024, wk, day_of_week, sep = "/"), 
                                          "Y/W/w"),
           wk_date = as.Date(wk_date),
           wk_doy = yday(wk_date))

# Just keep one observation of each plant, each week. Sort so the most
# advanced phenophase gets kept (if more than one value in a week)
df1 <- df %>%
  arrange(common_name, plant_id, year, wk, 
          desc(status_fl), desc(status_fo)) %>%
  distinct(plant_id, year, wk, .keep_all = TRUE) %>%
  # Remove week 53
  filter(wk != 53) %>%
  # Add an indicator for plants southcentral states
  mutate(sc = 1 * state %in% c("LA", "NM", "OK", "TX"))

# Create heat map -------------------------------------------------------------#
# Eventually, we'll want to automate this and cycle through species, but for
# now try going through the process for 1-2 species

i = 1

# Extract data for a species
spp <- species$common_name[i]
sppdat <- df1 %>% filter(common_name == spp)

# Isolate location information
sppsites <- sppdat %>%
  select(site_id, lat, lon, state, sc) %>%
  distinct() %>%
  mutate(sc = as.factor(sc))

# Plot locations of monitored plants
sppsitesv <- vect(sppsites, geom = c("lon", "lat"), crs = "epsg:4269")
states <- vect("data/states/cb_2017_us_state_500k.shp")
states <- subset(states, 
                 !states$STUSPS %in% c("HI", "AK", "VI", "MP", "GU", "PR", "AS"))
states$sc <- as.factor(ifelse(states$STUSPS %in% c("LA", "NM", "OK", "TX"), 1, 0))

ggplot(data = states) +
  geom_spatvector(aes(fill = sc), color = "gray30") +
  scale_fill_discrete(type = c("white", "gray80"), guide = "none") +
  geom_spatvector(data = sppsitesv, aes(color = sc), size = 1) +
  scale_color_discrete(type = c("steelblue3", "red"), guide = "none")

# Calculate how much data available in region, by year and across years
spp_nobs <- sppdat %>%
  group_by(year) %>%
  summarize(nplants = length(unique(plant_id)),
            nplants_sc = length(unique(plant_id[sc == 1]))) %>%
  data.frame()
spp_nobsperwk <- sppdat %>%
  filter(wk %in% 10:40) %>%
  group_by(year, wk) %>%
  summarize(nobs_open = n(), 
            nobs_open_sc = length(year[sc == 1]),
            .groups = "keep") %>%
  group_by(year) %>%
  summarize(nobs_weeklymn = mean(nobs_open),
            nobs_weeklymn_sc = mean(nobs_open_sc)) %>%
  data.frame()
spp_nobs <- spp_nobs %>%
  left_join(spp_nobsperwk, by = "year")

nobs_open <- sppdat %>%
  filter(wk %in% 10:40) %>%
  group_by(wk) %>%
  summarize(nobs_open = n(),
            nobs_open_sc = length(wk[sc == 1])) %>%
  data.frame()

overall <- data.frame(
  year = "all",
  nplants = length(unique(sppdat$plant_id)),
  nplants_sc = length(unique(sppdat$plant_id[sppdat$sc == 1])),
  nobs_weeklymn = mean(nobs_open$nobs_open),
  nobs_weeklymn_sc = mean(nobs_open$nobs_open_sc))

spp_nobs <- rbind(spp_nobs, overall)
spp_nobs

# For buttonbush, mean number of observations each week in the SC region is 
# < 5 for every year, but 10.3 across all years, so I think we only have that
# option. Will need to automate this in some way......

# Set threshold for mean number of plants observed each week in weeks 10-40
weekly_nobs_min <- 8

sppyears <- spp_nobs %>% filter(nobs_weeklymn_sc >= weekly_nobs_min)
if (nrow(sppyears) == 0) {
  warning("Data are insufficient to characterize open flower phenophase for ", 
          spp)
} else{
  if (nrow(sppyears) == 1) {
    message("Can only characterize open flower phenophase for ", spp, 
            " by combining data across all years")
  } else {
    message("Data are sufficient to characterize open flower phenophase for " ,
            spp, " in ", paste(sppyears$year[-nrow(sppyears)], collapse = ", "))
  }
}

