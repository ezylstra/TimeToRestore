library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
# library(shiny)
# library(bslib)
# library(sf)
library(ggplot2)
# library(leaflet)
# library(RColorBrewer)
# library(sortable)
# library(grid)
# library(cowplot)
library(terra)
library(tidyterra)
library(rnpn)

# Load/format observation data ------------------------------------------------#

df <- read.csv("data/status-intensity-flowers-May2026.csv")
df <- df %>%
  mutate(obsdate = ymd(obsdate)) %>%
  mutate(region = case_when(
    region == "Dallas" ~ "Dallas/Ft. Worth",
    region == "Austin" ~ "Austin/San Antonio",
    .default = region
  ))
# For now limit to 4 state area and 2025-2026
df <- df %>%
  filter(state4 == 1) %>%
  filter(yr >= 2025)

# If there's more than one observation of a plant on a given date, select one 
# with positive status and (higher) intensity value
df <- df %>%
  arrange(id, obsdate, php_id, desc(status), desc(midpoint)) %>%
  distinct(id, obsdate, php_id, .keep_all = TRUE)

# Put all data for a plant, date in the same row (wide form)
dfw <- df %>%
  select(-c(observation_id, php)) %>%
  pivot_wider(names_from = c(php_id),
              values_from = c(status, intensity_value, midpoint)) %>%
  data.frame()

# If flower = 0, then open flower must be 0. Remove observations with open
# flower = 1 and change any open flower = NA to 0.
dfw <- dfw %>%
  filter(!(!is.na(status_500) & !is.na(status_501) & status_500 == 0 & status_501 == 1)) %>%
  mutate(status_501 = ifelse(!is.na(status_500) & status_500 == 0, 0, status_501))

# If open flower = 1, then flower must be 1. Change any that are NA.
dfw <- dfw %>%
  mutate(status_500 = ifelse(!is.na(status_501) & status_501 == 1, 1, status_500))

# Finally, change any midpoint values to 0 when status is 0
dfw <- dfw %>%
  mutate(midpoint_500 = ifelse(!is.na(status_500) & status_500 == 0, 
                               0, midpoint_500)) %>%
  mutate(midpoint_501 = ifelse(!is.na(status_501) & status_501 == 0,
                               0, midpoint_501))

# Create table with location options
loc_options <- distinct(dfw, state, region) %>%
  mutate(region = ifelse(is.na(region), "Statewide", region)) %>%
  rbind(data.frame(state = "SC states (4)", 
                   region = "Entire SC region"))

# Desired minimum sample size for calculating proportions
min_sample_size <- 5

# Load iNat (via Phenobase) data ----------------------------------------------#
# iNat observations with flower annotation (don't have separate open flower
# annotation)

phenofolder <- "data/iNat/phenobase"
phenozips <- list.files(phenofolder, full.names = TRUE)

for (i in seq_along(phenozips)) {
  filename <- unz(phenozips[i], "data.csv")
  pheno1 <- read.csv(filename)
  if (i == 1) {
    inatfull <- pheno1
  } else {
    inatfull <- bind_rows(inatfull, pheno1)
  }
}
rm(pheno1)

# Simplify data
inat <- inatfull %>%
  mutate(date = ymd(eventDate)) %>%
  rename(coordU = coordinateUncertaintyInMeters) %>%
  select(scientificName, year, latitude, longitude, date, coordU)

# Look at what points have tons of geographic uncertainty
# inat %>%
#   mutate(uncertain = ifelse(is.na(coordU) | coordU > 10000, 1, 0)) %>%
#   ggplot() +
#   geom_point(aes(x = longitude, y = latitude, color = factor(uncertain)))
# No obvious patterns

# What does it mean if the coordU field is NA?
# For now will delete these, along with records that have very high values
# (location uncertainty > 10 km)
inat <- inat %>%
  filter(!is.na(latitude)) %>%
  filter(!is.na(coordU) & coordU <= 10000) %>%
  select(-coordU)

# Simplify species names (iNat/phenobase provides subspecies for some)
inat <- inat %>%
  separate_wider_delim(scientificName, 
                       delim = " ", 
                       names = c("genus", "species", "subspecies"),
                       too_few = "align_start")

# Get NPN species_ids and common names
npn_spp <- npn_species() %>%
  select(species_id, common_name, genus, species)
inat <- inat %>%
  left_join(npn_spp, by = c("genus", "species")) %>%
  data.frame()

# Assign states for each plant location
states <- vect("data/states/cb_2017_us_state_500k.shp")
states <- project(states, "epsg:4326")
states <- select(states, STUSPS)
inatlocs <- inat %>%
  select(longitude, latitude) %>%
  vect(crs = "epsg:4326")
statefill <- terra::extract(states, inatlocs)
inat$state <- statefill$STUSPS

# Extract observations for 4-state area
inat4 <- inat %>%
  filter(state %in% c("NM", "OK", "TX", "LA"))

# Create month, week, bi-weekly assignments
inat4 <- inat4 %>%
  mutate(month = month(date),
         wk = week(date)) %>%
  mutate(wk2 = ifelse(wk%%2 == 0, wk - 1, wk))

inat_month <- inat4 %>%
  group_by(common_name, month) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  group_by(common_name) %>%
  mutate(nmonths = n_distinct(month),
         totalobs = sum(nobs)) %>%
  ungroup() %>%
  filter(totalobs > 100) %>%
  data.frame()

# ggplot(inat_month, aes(x = month, y = nobs)) +
#   geom_smooth(method = "loess", span = 0.35, se = FALSE) +
#   geom_point() +
#   scale_x_continuous(breaks = seq(1, 11, by = 2)) +
#   scale_y_log10() +
#   facet_wrap(~common_name, scales = "free_y")

smoothed_month <- inat_month %>%
  group_by(common_name) %>%
  group_modify(function(.x, .y) {
    fit <- loess(log(nobs) ~ month, data = .x, span = 0.35)
    tibble(month = seq(min(.x$month), max(.x$month), length.out = 100)) %>%
      mutate(nobs = exp(predict(fit, newdata = .)))
  })

ggplot(inat_month, aes(x = month, y = nobs)) +
  geom_point() +
  geom_line(data = smoothed_month, color = "steelblue") +
  facet_wrap(~common_name, scales = "free_y")

# Biweekly
inat_biweek <- inat4 %>%
  group_by(common_name, wk2) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  group_by(common_name) %>%
  mutate(nperiods = n_distinct(wk2),
         totalobs = sum(nobs)) %>%
  ungroup() %>%
  filter(totalobs > 100) %>%
  data.frame()

smoothed_biweek <- inat_biweek %>%
  group_by(common_name) %>%
  group_modify(function(.x, .y) {
    fit <- loess(log(nobs) ~ wk2, data = .x, span = 0.25)
    tibble(wk2 = seq(min(.x$wk2), max(.x$wk2), length.out = 100)) %>%
      mutate(nobs = exp(predict(fit, newdata = .)))
  })

ggplot(inat_biweek, aes(x = wk2, y = nobs)) +
  geom_point() +
  geom_line(data = smoothed_biweek, color = "steelblue") +
  facet_wrap(~common_name, scales = "free_y")

# Weekly (seems like there's too much noise at weekly timescale)
# inat_week <- inat4 %>%
#   group_by(common_name, wk) %>%
#   summarize(nobs = n(), .groups = "drop") %>%
#   group_by(common_name) %>%
#   mutate(nweeks = n_distinct(wk),
#          totalobs = sum(nobs)) %>%
#   ungroup() %>%
#   filter(totalobs > 100) %>%
#   data.frame()
# 
# smoothed_week <- inat_week %>%
#   group_by(common_name) %>%
#   group_modify(function(.x, .y) {
#     fit <- loess(log(nobs) ~ wk, data = .x, span = 0.35)
#     tibble(wk = seq(min(.x$wk), max(.x$wk), length.out = 100)) %>%
#       mutate(nobs = exp(predict(fit, newdata = .)))
#   })
# 
# ggplot(inat_week, aes(x = wk, y = nobs)) +
#   geom_point() +
#   geom_line(data = smoothed_week, color = "steelblue") +
#   facet_wrap(~common_name, scales = "free_y")
