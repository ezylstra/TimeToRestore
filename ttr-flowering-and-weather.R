# Exploring options for models of weather effects on flowering onset
# 11 Aug 2026

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

# Create series dataset
dfs <- create_series(df)

# How much data do we have on each species?
dfs %>% filter(state4 == 1) %>% count(spp) %>% arrange(desc(n))

# For now, just keep species that have at least 50 "yes" series (combined total
# across two phenophases) in the 4-state region
dfs <- dfs %>%
  group_by(spp) %>%
  mutate(nseries = sum(state4)) %>%
  ungroup() %>%
  filter(nseries >= 50) %>%
  select(-nseries) %>% 
  data.frame()

# Series data, no restrictions on prior no
series <- dfs

# Series data, prior no within 14 days
series14 <- dfs %>%
  filter(!is.na(days_prior_no) & days_prior_no <= 14)
filter(series14, state4 == 1) %>% count(spp) %>% arrange(desc(n))

# Series data, prior no within 30 days
series30 <- dfs %>%
  filter(!is.na(days_prior_no) & days_prior_no <= 30)
filter(series30, state4 == 1) %>% count(spp) %>% arrange(desc(n))

# Need to find first yeses in some annual period (calendar for most, but 
# not all?). Lumping two phenophases together, but I think that's okay here.
ggplot(filter(series30, state4 == 1)) +
  geom_histogram(aes(x = first_yes_doy), binwidth = 2) +
  facet_wrap(~spp, ncol = 2)
# Red maple looks like the one species where water year would be better than
# calendar year (though many don't have much of a seasonal signal....)

# Extract first yeses in calendar year all species but red maple:
series <- series %>%
  mutate(first_yes_month = month(first_yes_date)) %>%
  mutate(first_yes_wateryr = ifelse(first_yes_month < 10,
                                    first_yes_year, 
                                    first_yes_year + 1)) %>%
  arrange(spp_id, id, php_id, first_yes_date)
firsts <- series %>%
  filter(spp != "red maple") %>%
  distinct(id, php_id, first_yes_year, .keep_all = TRUE)
firsts_rema <- series %>%
  filter(spp == "red maple") %>%
  distinct(id, php_id, first_yes_wateryr, .keep_all = TRUE)
firsts <- rbind(firsts, firsts_rema)

series14 <- series14 %>%
  mutate(first_yes_month = month(first_yes_date)) %>%
  mutate(first_yes_wateryr = ifelse(first_yes_month < 10,
                                    first_yes_year, 
                                    first_yes_year + 1)) %>%
  arrange(spp_id, id, php_id, first_yes_date)
firsts14 <- series14 %>%
  filter(spp != "red maple") %>%
  distinct(id, php_id, first_yes_year, .keep_all = TRUE)
firsts_rema14 <- series14 %>%
  filter(spp == "red maple") %>%
  distinct(id, php_id, first_yes_wateryr, .keep_all = TRUE)
firsts14 <- rbind(firsts14, firsts_rema14)

series30 <- series30 %>%
  mutate(first_yes_month = month(first_yes_date)) %>%
  mutate(first_yes_wateryr = ifelse(first_yes_month < 10,
                                    first_yes_year, 
                                    first_yes_year + 1)) %>%
  arrange(spp_id, id, php_id, first_yes_date)
firsts30 <- series30 %>%
  filter(spp != "red maple") %>%
  distinct(id, php_id, first_yes_year, .keep_all = TRUE)
firsts_rema30 <- series30 %>%
  filter(spp == "red maple") %>%
  distinct(id, php_id, first_yes_wateryr, .keep_all = TRUE)
firsts30 <- rbind(firsts30, firsts_rema30)

# Look at number of first yeses by species, phenophase
firsts30 %>%
  filter(state4 == 1) %>% 
  count(spp, php) %>%
  pivot_wider(id_cols = spp,
              names_from = php,
              values_from = n) %>%
  data.frame() %>%
  rename(open_flower = open.flower) %>%
  arrange(desc(flower + open_flower))

# TODO: 
# Download PRISM 2012, 2016 data
# Experiment with eastern redbud and beautyberry, then red maple
# Create simple dashboard with map and table with number of locations, 
# number of years, etc...
  # Identify what constitutes "sufficient" data for a species to be included in
  # subsequent analyses. Need a minimum number of unique combos: year, lat/lon 
  # (need the lat/lon to have different weather data though, so need a minimum 
  # distance between locations). Could attach cell ID from a PRISM raster to 
  # identify unique "weather locs".
# Try some climate variables:
  # Fixed 3-month mean temperatures, total precip (% of 30-year mean?)
  # variables summarized over previous 3 months
  # GDD better than temperature?
  # Include adjustment for latitude in all models

# Option to select variables OR be given "best model"?
# If best model option available, need to create a full model after removing 
# those that are highly correlated...

# What climate variables are we going to calculate?
# Fixed seasonal variables using 3-month NPN delineations (eg, mean spring temp)
# Weather summaries immediately preceding event (or average date of event)?
# For a given species, could find average date in SC region that flowering is
# first observed and then use 3-months prior (eg, if average first observation
# date is 14 April, then use mean temperatures for Jan-Mar)
