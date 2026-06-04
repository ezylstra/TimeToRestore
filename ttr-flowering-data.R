# Exploring options for flowering data dashboard
# 27 May 2026

require(lubridate)
require(stringr)
require(dplyr)
require(tidyr)
# require(terra)
require(sf)
require(ggplot2)
library(leaflet)
library(RColorBrewer)

# Load/format observation data ------------------------------------------------#

df <- read.csv("data/status-intensity-flowers-May2026.csv")
df <- df %>%
  mutate(obsdate = ymd(obsdate))

# Subset observation data -----------------------------------------------------#

# Are there any species-region combinations with very high number of plants?
df2526 <- df %>%
  filter(yr %in% 2025:2026) %>%
  filter(!is.na(region)) %>%
  group_by(spp, region) %>%
  summarize(n_plants = n_distinct(id),
            n_500 = sum(!is.na(status[php_id == 500])),
            n_501 = sum(!is.na(status[php_id == 501])),
            .groups = "drop") %>%
  data.frame()
filter(df2526, region == "Austin") %>% arrange(desc(n_plants))
filter(df2526, region == "Dallas") %>% arrange(desc(n_plants))
filter(df2526, region == "Houston") %>% arrange(desc(n_plants))

# Set focus region
df <- df %>%
  # filter(region == "Dallas") %>%
  filter(state == "TX")

# What is the year and species distribution like?
df %>% 
  group_by(yr) %>%
  summarize(n_sites = n_distinct(site),
            n_species = n_distinct(spp),
            n_500 = sum(php_id == 500),
            n_501 = sum(php_id == 501)) %>%
  data.frame()
table(df$spp, df$yr)

# This is what I originally did for Dallas area:
  # # Species-wise, there are only a couple worth looking at prior to 2025:
  #   # American beautyberry has decent number of observations for 2022-2026
  #   # Eastern redbud has decent number of observations in 2019-2025 and 2024-2026
  # # Otherwise, limit to species that have at least 50 observations in 2025-2026
  # # (note: this is probably way too generous, since counts combine 2 phenophases)
  # df <- df %>%
  #   group_by(spp) %>%
  #   mutate(nobs_2526 = sum(yr %in% 2025:2026)) %>% 
  #   ungroup() %>%
  #   mutate(keep = case_when(
  #     spp == "American beautyberry" & yr >= 2022 ~ 1,
  #     spp == "eastern redbud" & yr >= 2019 ~ 1,
  #     nobs_2526 >= 50 & yr >= 2025 ~ 1,
  #     .default = 0
  #   )) %>%
  #   filter(keep == 1) %>%
  #   select(-c(nobs_2526, keep))
  # table(df$spp, df$yr) 
  # # Total of 20 species left
  # 
  # df %>%
  #   filter(yr >= 2025) %>%
  #   group_by(spp) %>%
  #   summarize(nplants = n_distinct(id),
  #             n500 = sum(php_id == 500),
  #             n501 = sum(php_id == 501),
  #             n500_2025 = sum(php_id == 500 & yr == 2025),
  #             n501_2025 = sum(php_id == 501 & yr == 2025)) %>%
  #   data.frame() %>%
  #   arrange(desc(nplants))
  # 
  # # Going to play around with data for 9 species 
  # df9 <- df %>%
  #   filter(spp %in% c("American beautyberry",
  #                     "wax mallow", 
  #                     "eastern purple coneflower",
  #                     "green antelopehorn",
  #                     "wild bergamot",
  #                     "eastern redbud",
  #                     "mealycup sage",
  #                     "blue mistflower",
  #                     "Texas lupine")) %>%
  #   data.frame()

# Need to figure out what species are worth looking at. For each species and
# year, find number of plants or number of plants with at least 3 or 5
# observations?

sppyr <- df %>%
  group_by(id, php, yr) %>%
  mutate(nobs = n()) %>%
  ungroup() %>%
  group_by(spp, yr, php) %>%
  summarize(nplants = n_distinct(id),
            nplants3 = n_distinct(id[nobs >= 3]),
            nplants5 = n_distinct(id[nobs >= 5]),
            .groups = "drop") %>%
  data.frame()
# checks:
sppyr %>%
  filter(nplants5 >= 10) %>%
  arrange(php, spp, yr)

# Identify species that have at least 10 individuals with 5 or more observations
# in one or more years (10 seems reasonable when calculating proportions, though
# many weeks sample sizes may be much less than that). Only one species/year 
# before 2025, so will just grab 2025-2026 data for all (for now)
species10 <- sppyr %>%
  filter(nplants5 >= 10) %>%
  distinct(spp) %>% 
  pull(spp)

# Extract species data
df <- df %>%
  filter(spp %in% species10) %>%
  filter(yr >= 2025)

# Remove any duplicate observations (same plant, date, php) -------------------#

# If there's more than one observation, select one with positive status and 
# (higher) intensity value
df <- df %>%
  arrange(id, obsdate, php_id, desc(status), desc(midpoint)) %>%
  distinct(id, obsdate, php_id, .keep_all = TRUE)

# Remove any inconsistencies between phenophases/intensity values -------------#

# Put all data for a plant, date in the same row (wide form)
dfw <- df %>%
  select(-c(observation_id, php)) %>%
  pivot_wider(names_from = c(php_id),
              values_from = c(status, intensity_value, midpoint)) %>%
  data.frame()
# count(dfw, status_500, status_501, is.na(midpoint_500), is.na(midpoint_501))

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

# Map with site locations -----------------------------------------------------#

# Map to see where 2025-2026 observations were made?
locs <- dfw %>%
  group_by(site, site_name, lat, lon, spp) %>%
  summarize(plants = n_distinct(id), .groups = "drop") %>%
  data.frame()
# Jittering locations (slightly) in map data to see things better ####
# map_frame <- reactive({
  locs <- locs %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
    st_jitter(factor = 0.0001)
  locs <- locs %>%
    cbind(st_coordinates(locs)) %>%
    rename(lon = X, lat = Y) %>%
    st_drop_geometry()

# Get colors
nspp <- length(unique(locs$spp))
colors <- brewer.pal(nspp, "Paired")

# output$map <- renderLeaflet({
  # leaflet(data = map_frame()) %>%
   m <- leaflet(data = locs) %>%
    addTiles(options = tileOptions(opacity = 0.6))
   
  for (i in 1:nspp) {
    m <- m %>% addCircleMarkers(
      lng = ~lon, 
      lat = ~lat,
      data = filter(locs, spp == species10[i]),
      group = species10[i],
      radius = 5,
      color = "black",
      fillColor = colors[i], 
      fillOpacity = 0.9,
      stroke = TRUE,
      weight = 2,
      popup = ~paste0(spp, "<br>",
                      "Site ID: ", site, "<br>",
                      "No. plants: ", plants))
  }
  m %>%
    addLayersControl(overlayGroups = species10,
                             options = layersControlOptions(collapse = FALSE)) %>%
    addLegend(position = "bottomright",
              colors = colors,
              labels = species10,
              opacity = 0.9)
# })

# Data visualizations/summaries -----------------------------------------------#  
  
# When were (open flower) observations made (combining info across years)?
# Using open flowers observations because there are fewer flower observations
dfw %>%
  select(-contains("intensity_")) %>%
  select(-contains("midpoint_")) %>%
  pivot_longer(cols = status_500:status_501,
               names_to = "php",
               values_to = "status") %>%
  mutate(php = ifelse(php == "status_500", "flowers", "open flowers")) %>%
  filter(php == "open flowers") %>%
  ggplot() +
  geom_histogram(aes(x = doy), bins = 50, fill = "steelblue3") +
  facet_wrap(~spp, ncol = 2) +
  labs(x = "Day of year", y = "No. observations", title = "Observation dates")

# When did phenophases occur (combining info across years)?
dfw %>%
  select(-contains("intensity_")) %>%
  select(-contains("midpoint_")) %>%
  pivot_longer(cols = status_500:status_501,
               names_to = "php",
               values_to = "status") %>%
  mutate(php = ifelse(php == "status_500", "flowers", "open flowers")) %>%
  filter(status == 1) %>%
  ggplot() +
  geom_histogram(aes(x = doy, fill = php), bins = 50) +
  facet_grid(php ~ spp) +
  scale_fill_manual(values = c("coral2", "darkorchid3")) +
  labs(x = "Day of year", y = "No. observations", 
       title = "Phenophase occurrence (yes dates)") +
  theme(legend.position = "none")

# Select yeartype -------------------------------------------------------------#
##### Save this for later since it isn't important for weekly proportions 
##### Will be critical for thinking about onsets using first yes dates

# Is there a way to automate this?
# From Daijiang Li: We identified the largest annual temporal gap between unique 
# flowering dates and designated the midpoint of this gap as the "phenological 
# break." Observations occurring after the calendar year reset but before this 
# break were re-indexed by adding 365 to their original DOY (e.g., a flowering 
# event on Day 10 of the following year would be rendered as Day 375). This 
# linearization of circular data ensured that reproductive events were treated 
# as continuous distributions rather than fragmented segments.

# There are a couple choices to make:
# 1) Summarize by individual, species, or across species
# 2) Use date (including year) or day of year (ignoring year)

# # Individual, date
# gapsfl <- dfw %>%
#   filter(yr >= 2025) %>%
#   filter(!is.na(status_500) & status_500 == 1) %>%
#   group_by(id) %>%
#   mutate(bothyrs = ifelse(n_distinct(yr) == 2, 1, 0)) %>%
#   ungroup() %>%
#   filter(bothyrs == 1) %>%
#   arrange(id, obsdate) %>%
#   mutate(gap = NA)
# for (i in 2:nrow(gapsfl)) {
#   gapsfl$gap[i] <- ifelse(gapsfl$id[i] != gapsfl$id[i-1], NA,
#                           as.numeric(gapsfl$obsdate[i] - gapsfl$obsdate[i-1]))
# }
# # Identify max gap for each plant and when it occurred
# gapsfl <- gapsfl %>%
#   group_by(id) %>%
#   mutate(max_gap = ifelse(n() == 1, NA, max(gap, na.rm = TRUE))) %>%
#   ungroup() %>%
#   data.frame()
# # Ignore any instances where max gap < 14 days
# gapsfl <- gapsfl %>%
#   filter(max_gap >= 14)
# # Calculate previous doy (and note whether this was in last calendar year
# # or current day of year)
# gapsfl <- gapsfl %>%
#   filter(gap == max_gap) %>%
#   select(spp, id, doy, obsdate, gap) %>%
#   arrange(spp, id) %>%
#   mutate(doy_prior = doy - gap) %>%
#   mutate(before_gap = ifelse(doy_prior < 0, "last year", "current year"))
# gapsfl
# # Summarize by species
# gapsfl_spp <- gapsfl %>%
#   group_by(spp) %>%
#   summarize(nplants = n_distinct(id),
#             newyeargap = sum(before_gap == "last year")) %>%
#   mutate(prop_calyear = round(newyeargap / nplants, 2)) %>%
#   data.frame()
# 
# # Species, date
# gapsfl <- dfw %>%
#   filter(yr >= 2025) %>%
#   filter(!is.na(status_500) & status_500 == 1) %>%
#   group_by(spp) %>%
#   mutate(bothyrs = ifelse(n_distinct(yr) == 2, 1, 0)) %>%
#   ungroup() %>%
#   filter(bothyrs == 1) %>%
#   arrange(spp, obsdate) %>%
#   mutate(gap = NA)
# for (i in 2:nrow(gapsfl)) {
#   gapsfl$gap[i] <- ifelse(gapsfl$spp[i] != gapsfl$spp[i-1], NA,
#                           as.numeric(gapsfl$obsdate[i] - gapsfl$obsdate[i-1]))
# }
# # Identify max gap for each species and when it occurred
# gapsfl <- gapsfl %>%
#   group_by(spp) %>%
#   mutate(max_gap = ifelse(n() == 1, NA, max(gap, na.rm = TRUE))) %>%
#   ungroup() %>%
#   data.frame()
# # Calculate previous doy (and note whether this was in last calendar year
# # or current day of year)
# gapsfl_spp <- gapsfl %>%
#   filter(gap == max_gap) %>%
#   select(spp, doy, obsdate, gap) %>%
#   arrange(spp) %>%
#   mutate(doy_prior = doy - gap) %>%
#   mutate(before_gap = ifelse(doy_prior < 0, "last year", "current year"))

# Species, day of year
gapsfl <- dfw %>%
  filter(yr >= 2025) %>%
  filter(!is.na(status_500) & status_500 == 1) %>%
  arrange(spp, doy) %>%
  distinct(spp, doy)
for (i in 1:length(species10)) {
  gapsfl_temp <- filter(gapsfl, spp == species10[i])
  gapsfl_temp$doyprior <- c(last(gapsfl_temp$doy),
                            gapsfl_temp$doy[1:(nrow(gapsfl_temp) - 1)])
  gapsfl_temp$doy2 <- gapsfl_temp$doy
  gapsfl_temp$doy2[1] <- gapsfl_temp$doy[1] + 365
  gapsfl_temp$gapprior <- gapsfl_temp$doy2 - gapsfl_temp$doyprior
  if (i == 1) {
    gapsfl2 <- gapsfl_temp
  } else {
    gapsfl2 <- rbind(gapsfl2, gapsfl_temp)
  }
}
# Identify max gap for each species and when it occurred
gapsfl_spp <- gapsfl2 %>%
  group_by(spp) %>%
  mutate(max_gap = ifelse(n() == 1, NA, max(gapprior, na.rm = TRUE))) %>%
  ungroup() %>%
  filter(gapprior == max_gap) %>%
  mutate(overlap_jan1 = ifelse(doy < doyprior, 1, 0)) %>%
  mutate(breakdate = ceiling((doy2 - doyprior)/2) + doyprior) %>%
  # If gap spans, Jan 1, use Jan 1 as break date; if break date > 365, 
  # substract 365 (don't think I need the 2nd part, but keeping in case)
  mutate(breakdate = case_when(
    overlap_jan1 == 1 ~ 1,
    breakdate > 365 ~ breakdate - 365,
    .default = breakdate
  )) %>%
  select(-c(doy2, gapprior)) %>%
  data.frame()
gapsfl_spp

# # When did phenophases occur (combining info across species and years)?
# dfw %>%
#   select(-contains("intensity_")) %>%
#   select(-contains("midpoint_")) %>%
#   pivot_longer(cols = status_500:status_501,
#                names_to = "php",
#                values_to = "status") %>%
#   mutate(php = ifelse(php == "status_500", "flowers", "open flowers")) %>%
#   filter(status == 1) %>%
#   ggplot() +
#   geom_histogram(aes(x = doy, fill = php), bins = 50) +
#   facet_wrap(~php, ncol = 1) +
#   scale_fill_manual(values = c("coral2", "darkorchid3")) +
#   labs(x = "Day of year", y = "No. observations", 
#        title = "Phenophase occurrence (yes dates)") +
#   theme(legend.position = "none")
# 
# # All species, day of year
# # If we want to do this, probably best to summarize with weekly data since there
# # will be potentially many days with very few or 1 observations...

#### Print message about yeartype with phenophase occurrence histograms?

# Assign week number to each observation --------------------------------------#

# wk = number of complete 7-day periods since Jan 1 (so Jan 7 always = wk 1)
# wk_doy1 = start of each week (eg, date for week 1 would be Jan 1)
# wk_doy4 = middle of each week (eg, date for week 1 would be Jan 4)
dfw <- dfw %>%
  mutate(wk = week(obsdate)) %>%
  # Remove observations in week 53 (Dec 31 [and Dec 30 in leap years])
  filter(wk < 53) %>%
  # Create wk_doy columns
  mutate(wk_doy1 = (wk * 7) - 6,
         wk_doy4 = (wk * 7) - 3) %>%
  mutate(wk_date1 = parse_date_time(x = paste(2025, wk_doy1), orders = "yj"),
         wk_date4 = parse_date_time(x = paste(2025, wk_doy4), orders = "yj"))

# Assign bi-weekly period 
dfw <- dfw %>%
  mutate(wk2 = ifelse(wk%%2 == 0, wk - 1, wk))

# Set parameters --------------------------------------------------------------#

# Period 
# weekly = TRUE: calculate proportion of plants in phase each week
# weekly = FALSE: calculate proportion of plants in phase every 2-week period
weekly <- FALSE

# Year or multi-year calculations
# yearly = TRUE: calculate separate proportions each year
# yearly = FALSE: calculate proportions by combining data across years 
yearly <- FALSE
#### Instead of this, could set yearpick = "2025", "2026", or "Combined"
#### If we did this, would need to create a year column that works if NOT  
#### using calendar year (eg, calendar year at start of annual "period")

# Calculate proportions -------------------------------------------------------#

# Set period as week or bi-week
dfw <- dfw %>%
  rowwise() %>%
  mutate(period = ifelse(weekly == TRUE, wk, wk2)) %>%
  data.frame()

# Separate out flowers, open flowers data
flowers <- dfw %>%
  select(-contains("501")) %>%
  rename(status = status_500,
         intensity_value = intensity_value_500,
         midpoint = midpoint_500) %>%
  filter(!is.na(status)) 
open <- dfw %>%
  select(-contains("500")) %>%
  rename(status = status_501,
         intensity_value = intensity_value_501,
         midpoint = midpoint_501) %>%
  filter(!is.na(status)) 

# Keep just one observation of each plant, each week or 2-week period. Sort so
# the most advanced phenophase gets kept (if more than one value in period)
flowersp <- flowers %>%
  arrange(id, yr, period, desc(status), desc(midpoint)) %>%
  distinct(id, yr, period, .keep_all = TRUE)
openp <- open %>%
  arrange(id, yr, period, desc(status), desc(midpoint)) %>%
  distinct(id, yr, period, .keep_all = TRUE)

# Calculate proportions
if (yearly == TRUE) {
  flowerprops <- flowersp %>%
    group_by(spp, yr, period, wk_doy1, wk_doy4, wk_date1, wk_date4) %>%
    summarize(nplants = n(),
              nyes = sum(status),
              prop = nyes/nplants,
              .groups = "drop") %>%
    mutate(year = as.character(yr)) %>%
    data.frame()
  openprops <- openp %>%
    group_by(spp, yr, period, wk_doy1, wk_doy4, wk_date1, wk_date4) %>%
    summarize(nplants = n(),
              nyes = sum(status),
              prop = nyes/nplants,
              .groups = "drop") %>%
    mutate(year = as.character(yr)) %>%
    data.frame()
} else {
  flowerprops <- flowersp %>%
    group_by(spp, period, wk_doy1, wk_doy4, wk_date1, wk_date4) %>%
    summarize(nplants = n(),
              nyes = sum(status),
              prop = nyes/nplants,
              .groups = "drop") %>%
    mutate(year = "Combined") %>%
    data.frame()
  openprops <- openp %>%
    group_by(spp, period, wk_doy1, wk_doy4, wk_date1, wk_date4) %>%
    summarize(nplants = n(),
              nyes = sum(status),
              prop = nyes/nplants,
              .groups = "drop") %>%
    mutate(year = "Combined") %>%
    data.frame()
}

# First look. Sample sizes sufficient?
statusfills <- c("No" = "gray", "Yes" = "steelblue3")
flowerprops %>%
  ggplot(aes(x = period, y = nplants)) +
  geom_bar(stat = "identity", width = 0.5, aes(fill = "No")) +
  geom_bar(aes(y = nyes, fill = "Yes"), stat = "identity", width = 0.5) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "coral3", linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dotted", color = "coral3", linewidth = 1) +
  facet_wrap(~spp) +
  scale_fill_manual(values = statusfills) +
  labs(x = "Week", y = "No. plants", fill = "Status report",
       title = "Flowers") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85, 0.2))

openprops %>%
  ggplot(aes(x = period, y = nplants)) +
  geom_bar(stat = "identity", width = 0.5, aes(fill = "No")) +
  geom_bar(aes(y = nyes, fill = "Yes"), stat = "identity", width = 0.5) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "coral3", linewidth = 0.5) +
  geom_hline(yintercept = 10, linetype = "dotted", color = "coral3", linewidth = 0.5) +
  facet_wrap(~spp) +
  scale_fill_manual(values = statusfills) +
  labs(x = "Week", y = "No. plants", fill = "Status report",
       title = "Open flowers") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85, 0.2),
        panel.grid = element_blank())

# Species with largest bi-weekly sample sizes?
# American beautyberry
# eastern redbud
# ~wax mallow
# ~eastern purple coneflower
  
# Other ways to visualize these data (offer users an option?)
  # Heatmap 
  # Bubble figure with various options for average/smoothed curve

# Calculate summary statistics (estimated "onset/offset/duration/peak" based
# on these visualizations/analyses)? 

# Combine phenophase into to get open flower counts?

# Compare with iNat observations?



