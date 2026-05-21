# Exploring options for flowering data dashboard
# 21 May 2026

require(lubridate)
require(stringr)
require(dplyr)
require(tidyr)
# require(terra)
require(ggplot2)

# Load observation data -------------------------------------------------------#

df <- read.csv("data/status-intensity-flowers-May2026.csv")

# Subset observation data -----------------------------------------------------#

# Will focus on Dallas region for now
df <- df %>%
  filter(region == "Dallas")

# What is the year and species distribution like?
df %>% 
  group_by(yr) %>%
  summarize(n_sites = n_distinct(site),
            n_species = n_distinct(spp),
            n_500 = sum(php_id == 500),
            n_501 = sum(php_id == 501)) %>%
  data.frame()
table(df$spp, df$yr)

# Species-wise, there are only a couple worth looking at prior to 2025:
  # American beautyberry has decent number of observations for 2022-2026
  # Eastern redbud has decent number of observations in 2019-2025 and 2024-2026
# Otherwise, limit to species that have at least 50 observations in 2025-2026
# (note: this is probably way too generous, since counts combine 2 phenophases)
df <- df %>%
  group_by(spp) %>%
  mutate(nobs_2526 = sum(yr %in% 2025:2026)) %>% 
  ungroup() %>%
  mutate(keep = case_when(
    spp == "American beautyberry" & yr >= 2022 ~ 1,
    spp == "eastern redbud" & yr >= 2019 ~ 1,
    nobs_2526 >= 50 & yr >= 2025 ~ 1,
    .default = 0
  )) %>%
  filter(keep == 1) %>%
  select(-c(nobs_2526, keep))
table(df$spp, df$yr) 
# Total of 20 species left

df %>%
  filter(yr >= 2025) %>%
  group_by(spp) %>%
  summarize(n500 = sum(php_id == 500),
            n501 = sum(php_id == 501),
            n500_2025 = sum(php_id == 500 & yr == 2025),
            n501_2025 = sum(php_id == 501 & yr == 2025)) %>%
  data.frame() %>%
  arrange(desc(n500))

# Going to play around with data for 5 species 
df5 <- df %>%
  filter(spp %in% c("American beautyberry",
                    "wax mallow", 
                    "eastern purple coneflower",
                    "green antelopehorn",
                    "wild bergamot")) %>%
  data.frame()

# Remove any duplicate observations (same plant, date, php) -------------------#

# If there's more than one observation, select one with positive status and 
# (higher) intensity value
df5 <- df5 %>%
  arrange(id, obsdate, php_id, desc(status), desc(midpoint)) %>%
  distinct(id, obsdate, php_id, .keep_all = TRUE)

# Remove any inconsistencies between phenophases/intensity values -------------#

# Put all data for a plant, date in the same row (wide form)
df5w <- df5 %>%
  select(-c(observation_id, php)) %>%
  pivot_wider(names_from = c(php_id),
              values_from = c(status, intensity_value, midpoint)) %>%
  data.frame()

count(df5w, status_500, status_501, is.na(midpoint_500), is.na(midpoint_501))

# If flower = 0, then open flower must be 0. Remove observations with open
# flower = 1 and change any open flower = NA to 0.
df5w <- df5w %>%
  filter(!(!is.na(status_500) & !is.na(status_501) & status_500 == 0 & status_501 == 1)) %>%
  mutate(status_501 = ifelse(!is.na(status_500) & status_500 == 0, 0, status_501))

# If open flower = 1, then flower must be 1. Change any that are NA.
df5w <- df5w %>%
  mutate(status_500 = ifelse(!is.na(status_501) & status_501 == 1, 1, status_500))

# Finally, change any midpoint values to 0 when status is 0
df5w <- df5w %>%
  mutate(midpoint_500 = ifelse(!is.na(status_500) & status_500 == 0, 
                               0, midpoint_500)) %>%
  mutate(midpoint_501 = ifelse(!is.na(status_501) & status_501 == 0,
                               0, midpoint_501))

# Data visualizations/summaries -----------------------------------------------#

# Simple histograms to see when observations were made
# Map to see where observations were made?

# Using all status data (not just first yes/onset info)
  # Simple histograms with flower/open flower dates
    # Use this to decide on yeartype?
  # Proportion of positive observations have intensity values (by phenophase)
  # Proportion of plants in phase (by year or combined across years)
    # Heatmap 
    # Bubble figure with various options for average/smoothed curve

# Calculate summary statistics (estimated "onset/offset/duration/peak" based
# on these visualizations/analyses)?


