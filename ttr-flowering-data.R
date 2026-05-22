# Exploring options for flowering data dashboard
# 22 May 2026

require(lubridate)
require(stringr)
require(dplyr)
require(tidyr)
# require(terra)
require(ggplot2)
library(leaflet)

# Load/format observation data ------------------------------------------------#

df <- read.csv("data/status-intensity-flowers-May2026.csv")
df <- df %>%
  mutate(obsdate = ymd(obsdate))

# Subset observation data -----------------------------------------------------#

# Are are any species-region combinations with notable number of plants?
df2526 <- df %>%
  filter(yr %in% 2025:2026) %>%
  filter(!is.na(region)) %>%
  group_by(spp, region) %>%
  summarize(n_plants = n_distinct(id),
            n_500 = sum(!is.na(status[php_id == 500])),
            n_501 = sum(!is.na(status[php_id == 501])),
            .groups = "drop") %>%
  data.frame()
arrange(df2526, region, desc(n_plants))

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
  summarize(nplants = n_distinct(id),
            n500 = sum(php_id == 500),
            n501 = sum(php_id == 501),
            n500_2025 = sum(php_id == 500 & yr == 2025),
            n501_2025 = sum(php_id == 501 & yr == 2025)) %>%
  data.frame() %>%
  arrange(desc(nplants))

# Going to play around with data for 9 species 
df9 <- df %>%
  filter(spp %in% c("American beautyberry",
                    "wax mallow", 
                    "eastern purple coneflower",
                    "green antelopehorn",
                    "wild bergamot",
                    "eastern redbud",
                    "mealycup sage",
                    "blue mistflower",
                    "Texas lupine")) %>%
  data.frame()

# Remove any duplicate observations (same plant, date, php) -------------------#

# If there's more than one observation, select one with positive status and 
# (higher) intensity value
df9 <- df9 %>%
  arrange(id, obsdate, php_id, desc(status), desc(midpoint)) %>%
  distinct(id, obsdate, php_id, .keep_all = TRUE)

# Remove any inconsistencies between phenophases/intensity values -------------#

# Put all data for a plant, date in the same row (wide form)
df9w <- df9 %>%
  select(-c(observation_id, php)) %>%
  pivot_wider(names_from = c(php_id),
              values_from = c(status, intensity_value, midpoint)) %>%
  data.frame()
# count(df9w, status_500, status_501, is.na(midpoint_500), is.na(midpoint_501))

# If flower = 0, then open flower must be 0. Remove observations with open
# flower = 1 and change any open flower = NA to 0.
df9w <- df9w %>%
  filter(!(!is.na(status_500) & !is.na(status_501) & status_500 == 0 & status_501 == 1)) %>%
  mutate(status_501 = ifelse(!is.na(status_500) & status_500 == 0, 0, status_501))

# If open flower = 1, then flower must be 1. Change any that are NA.
df9w <- df9w %>%
  mutate(status_500 = ifelse(!is.na(status_501) & status_501 == 1, 1, status_500))

# Finally, change any midpoint values to 0 when status is 0
df9w <- df9w %>%
  mutate(midpoint_500 = ifelse(!is.na(status_500) & status_500 == 0, 
                               0, midpoint_500)) %>%
  mutate(midpoint_501 = ifelse(!is.na(status_501) & status_501 == 0,
                               0, midpoint_501))

# Map with site locations -----------------------------------------------------#

# Map to see where 2025-2026 observations were made?
locs9 <- df9w %>%
  filter(yr >= 2025) %>%
  group_by(site, site_name, lat, lon, spp) %>%
  summarize(plants = n_distinct(id), .groups = "drop") %>%
  data.frame()

# output$map <- renderLeaflet({
  # leaflet(data = map_frame()) %>%
  leaflet(data = locs9) %>%
    addTiles(options = tileOptions(opacity = 0.7)) %>%
    addCircleMarkers(
      lng = ~lon, 
      lat = ~lat,
      # data = filter(map_frame(), spp == "American beautyberry"),
      data = filter(locs9, spp == "American beautyberry"),
      group = "American beautyberry",
      radius = 5,
      fillColor = "blue", 
      fillOpacity = 0.7,
      stroke = FALSE,
      popup = ~paste0(spp, "<br>",
                      "Site ID: ", site, "<br>",
                      "Site name: ", site_name, "<br>",
                      "No. plants: ", plants)
    ) %>%
    addCircleMarkers(
      lng = ~lon, 
      lat = ~lat,
      # data = filter(map_frame(), spp == "wax mallow"),
      data = filter(locs9, spp == "wax mallow"),
      group = "wax mallow",
      radius = 5,
      fillColor = "orange", 
      fillOpacity = 0.7,
      stroke = FALSE,
      popup = ~paste0(spp, "<br>",
                      "Site ID: ", site, "<br>",
                      "Site name: ", site_name, "<br>",
                      "No. plants: ", plants)
    ) %>%
    addCircleMarkers(
        lng = ~lon, 
        lat = ~lat,
        # data = filter(map_frame(), spp == "eastern purple coneflower"),
        data = filter(locs9, spp == "eastern purple coneflower"),
        group = "eastern purple coneflower",
        radius = 5,
        fillColor = "purple", 
        fillOpacity = 0.7,
        stroke = FALSE,
        popup = ~paste0(spp, "<br>",
                        "Site ID: ", site, "<br>",
                        "Site name: ", site_name, "<br>",
                        "No. plants: ", plants)
      ) %>%
    addCircleMarkers(
      lng = ~lon, 
      lat = ~lat,
      # data = filter(map_frame(), spp == "wild bergamot"),
      data = filter(locs9, spp == "wild bergamot"),
      group = "wild bergamot",
      radius = 5,
      fillColor = "red", 
      fillOpacity = 0.7,
      stroke = FALSE,
      popup = ~paste0(spp, "<br>",
                      "Site ID: ", site, "<br>",
                      "Site name: ", site_name, "<br>",
                      "No. plants: ", plants)
    ) %>%
    addCircleMarkers(
      lng = ~lon, 
      lat = ~lat,
      # data = filter(map_frame(), spp == "green antelopehorn"),
      data = filter(locs9, spp == "green antelopehorn"),
      group = "green antelopehorn",
      radius = 5,
      fillColor = "green", 
      fillOpacity = 0.7,
      stroke = FALSE,
      popup = ~paste0(spp, "<br>",
                      "Site ID: ", site, "<br>",
                      "Site name: ", site_name, "<br>",
                      "No. plants: ", plants)
    ) %>%
    addCircleMarkers(
      lng = ~lon, 
      lat = ~lat,
      # data = filter(map_frame(), spp == "Texas lupine"),
      data = filter(locs9, spp == "Texas lupine"),
      group = "Texas lupine",
      radius = 5,
      fillColor = "yellow", 
      fillOpacity = 0.7,
      stroke = FALSE,
      popup = ~paste0(spp, "<br>",
                      "Site ID: ", site, "<br>",
                      "Site name: ", site_name, "<br>",
                      "No. plants: ", plants)
    ) %>%
    addCircleMarkers(
      lng = ~lon, 
      lat = ~lat,
      # data = filter(map_frame(), spp == "blue mistflower"),
      data = filter(locs9, spp == "blue mistflower"),
      group = "blue mistflower",
      radius = 5,
      fillColor = "skyblue", 
      fillOpacity = 0.7,
      stroke = FALSE,
      popup = ~paste0(spp, "<br>",
                      "Site ID: ", site, "<br>",
                      "Site name: ", site_name, "<br>",
                      "No. plants: ", plants)
    ) %>%
    addCircleMarkers(
      lng = ~lon, 
      lat = ~lat,
      # data = filter(map_frame(), spp == "eastern redbud"),
      data = filter(locs9, spp == "eastern redbud"),
      group = "eastern redbud",
      radius = 5,
      fillColor = "pink", 
      fillOpacity = 0.7,
      stroke = FALSE,
      popup = ~paste0(spp, "<br>",
                      "Site ID: ", site, "<br>",
                      "Site name: ", site_name, "<br>",
                      "No. plants: ", plants)
    ) %>%
    addCircleMarkers(
      lng = ~lon, 
      lat = ~lat,
      # data = filter(map_frame(), spp == "mealycup sage"),
      data = filter(locs9, spp == "mealycup sage"),
      group = "mealycup sage",
      radius = 5,
      fillColor = "gray", 
      fillOpacity = 0.7,
      stroke = FALSE,
      popup = ~paste0(spp, "<br>",
                      "Site ID: ", site, "<br>",
                      "Site name: ", site_name, "<br>",
                      "No. plants: ", plants)
    ) %>%
    addLayersControl(overlayGroups = c("American beautyberry",
                                       "wax mallow",
                                       "eastern purple coneflower",
                                       "wild bergamot",
                                       "green antelopehorn",
                                       "Texas lupine",
                                       "blue mistflower",
                                       "eastern redbud",
                                       "mealycup sage"),
                     options = layersControlOptions(collapse = FALSE)) %>%
    addLegend(position = "bottomright",
              colors = c("blue", "orange", "purple", "red", "green",
                         "yellow", "skyblue", "pink", "gray"),
              labels = c("American beautyberry",
                         "wax mallow",
                         "eastern purple coneflower",
                         "wild bergamot",
                         "green antelopehorn",
                         "Texas lupine",
                         "blue mistflower",
                         "eastern redbud",
                         "mealycup sage"),
              opacity = 1)
# })
  
# Note that we're not worrying about species overlap right now (ie, species
# will completely overlap when at the same site)

# Data visualizations/summaries -----------------------------------------------#  
  
# When were (open flower) observations made (combining info across years)?
# Using open flowers observations because there are fewer flower observations
df9w %>%
  select(-contains("intensity_")) %>%
  select(-contains("midpoint_")) %>%
  pivot_longer(cols = status_500:status_501,
               names_to = "php",
               values_to = "status") %>%
  mutate(php = ifelse(php == "status_500", "flowers", "open flowers")) %>%
  ggplot() +
  geom_histogram(aes(x = doy, fill = php), bins = 50, fill = "steelblue3") +
  facet_wrap(~spp, ncol = 2) +
  labs(x = "Day of year", y = "No. observations", title = "Observation dates")

# When did phenophases occur (combining info across years)?
df9w %>%
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
# By species? For all species in region?

# Is there a way to automate this?
# From Daijiang Li: We identified the largest annual temporal gap between unique 
# flowering dates and designated the midpoint of this gap as the "phenological 
# break." Observations occurring after the calendar year reset but before this 
# break were re-indexed by adding 365 to their original DOY (e.g., a flowering 
# event on Day 10 of the following year would be rendered as Day 375). This 
# linearization of circular data ensured that reproductive events were treated 
# as continuous distributions rather than fragmented segments.

# Try this at an individual level?

# Look at gaps for plants that were observed in 2025 and 2026
gapsfl <- df9w %>%
  filter(yr >= 2025) %>%
  filter(!is.na(status_500) & status_500 == 1) %>%
  group_by(id) %>%
  mutate(bothyrs = ifelse(n_distinct(yr) == 2, 1, 0)) %>%
  ungroup() %>%
  filter(bothyrs == 1) %>%
  arrange(id, obsdate) %>%
  mutate(gap = NA)
for (i in 2:nrow(gapsfl)) {
  gapsfl$gap[i] <- ifelse(gapsfl$id[i] != gapsfl$id[i-1], NA,
                          as.numeric(gapsfl$obsdate[i] - gapsfl$obsdate[i-1]))
}
# Identify max gap for each plant and when it occurred
gapsfl <- gapsfl %>%
  group_by(id) %>%
  mutate(max_gap = ifelse(n() == 1, NA, max(gap, na.rm = TRUE))) %>%
  ungroup() %>%
  data.frame()
# Calculate previous doy (and note whether this was in last calendar year
# or current day of year)
gapsfl %>% 
  filter(gap == max_gap) %>%
  select(spp, id, doy, obsdate, gap) %>%
  arrange(spp, id) %>%
  mutate(doy_prior = doy - gap) %>%
  mutate(before_gap = ifelse(doy_prior < 0, "last year", "current year"))

# Only 20 plants, but all 9 species represented, and biggest gap for all 
# spanned new year, so calendar year good.

# Or, should we do this at a species level?
  # gapsfl2 <- df9w %>%
  #   filter(yr >= 2025) %>%
  #   filter(!is.na(status_500) & status_500 == 1) %>%
  #   distinct(spp, obsdate, doy) %>%
  #   arrange(spp, obsdate) %>%
  #   mutate(gap = NA)
  # for (i in 2:nrow(gapsfl2)) {
  #   gapsfl2$gap[i] <- ifelse(gapsfl2$spp[i] != gapsfl2$spp[i-1], NA,
  #                           as.numeric(gapsfl2$obsdate[i] - gapsfl2$obsdate[i-1]))
  # }
  # # Identify max gap for each species and when it occurred
  # gapsfl2 <- gapsfl2 %>%
  #   group_by(spp) %>%
  #   mutate(max_gap = ifelse(n() == 1, NA, max(gap, na.rm = TRUE))) %>%
  #   ungroup() %>%
  #   data.frame()
  # # Calculate previous doy (and note whether this was in last calendar year
  # # or current day of year)
  # gapsfl2 %>%
  #   filter(gap == max_gap) %>%
  #   select(spp, doy, obsdate, gap) %>%
  #   arrange(spp) %>%
  #   mutate(doy_prior = doy - gap) %>%
  #   mutate(before_gap = ifelse(doy_prior < 0, "last year", "current year"))
# Not sure this is great at the species level since one anomalous positive 
# observation near the turn of the year can minimize the new-year gap and cause
# the maximum gap for that species to occur at another time.

# Weekly proportions of plants in phase ---------------------------------------#

# We'll create wk and wk_doy columns:
# wk = number of complete 7-day periods since Jan 1 (so Jan 7 always = wk 1)
# wk_doy1 = start of each week (eg, date for week 1 would be Jan 1)
# wk_doy4 = middle of each week (eg, date for week 1 would be Jan 4)
df9w <- df9w %>%
  mutate(wk = week(obsdate)) %>%
  # Remove observations in week 53 (Dec 31 [and Dec 30 in leap years])
  filter(wk < 53) %>%
  # Create wk_doy columns
  mutate(wk_doy1 = (wk * 7) - 6,
         wk_doy4 = (wk * 7) - 3) %>%
  mutate(wk_date1 = parse_date_time(x = paste(2025, wk_doy1), orders = "yj"),
         wk_date4 = parse_date_time(x = paste(2025, wk_doy4), orders = "yj"))

# Separate out flowers, open flowers data
flowers <- df9w %>%
  select(-contains("501")) %>%
  rename(status = status_500,
         intensity_value = intensity_value_500,
         midpoint = midpoint_500) %>%
  filter(!is.na(status)) 
open <- df9w %>%
  select(-contains("500")) %>%
  rename(status = status_501,
         intensity_value = intensity_value_501,
         midpoint = midpoint_501) %>%
  filter(!is.na(status)) 

# Keep just one observation of each plant, each week. Sort so the most advanced 
# phenophase gets kept (if more than one value in a week)
flowerswk <- flowers %>%
  arrange(id, yr, wk, desc(status), desc(midpoint)) %>%
  distinct(id, yr, wk, .keep_all = TRUE)
openwk <- open %>%
  arrange(id, yr, wk, desc(status), desc(midpoint)) %>%
  distinct(id, yr, wk, .keep_all = TRUE)

# Calculate weekly proportions - ACROSS YEARS
flowerswk <- flowerswk %>%
  group_by(spp, wk, wk_doy1, wk_doy4, wk_date1, wk_date4) %>%
  summarize(nplants = n(),
            nyes = sum(status),
            prop = nyes/nplants,
            .groups = "drop") %>%
  data.frame()
openwk <- openwk %>%
  group_by(spp, wk, wk_doy1, wk_doy4, wk_date1, wk_date4) %>%
  summarize(nplants = n(),
            nyes = sum(status),
            prop = nyes/nplants,
            .groups = "drop") %>%
  data.frame()

# Are sample sizes sufficient?
statusfills <- c("No" = "gray", "Yes" = "steelblue3")
flowerswk %>%
  ggplot(aes(x = wk, y = nplants)) +
  geom_bar(stat = "identity", width = 0.5, aes(fill = "No")) +
  geom_bar(aes(y = nyes, fill = "Yes"), stat = "identity", width = 0.5) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "coral3", linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dotted", color = "coral3", linewidth = 1) +
  facet_wrap(~spp) +
  scale_fill_manual(values = statusfills) +
  labs(x = "Week", y = "No. plants", fill = "Status report",
       title = "Flowers, all years combined") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85, 0.2))

openwk %>%
  ggplot(aes(x = wk, y = nplants)) +
  geom_bar(stat = "identity", width = 0.5, aes(fill = "No")) +
  geom_bar(aes(y = nyes, fill = "Yes"), stat = "identity", width = 0.5) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "coral3", linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dotted", color = "coral3", linewidth = 1) +
  facet_wrap(~spp) +
  scale_fill_manual(values = statusfills) +
  labs(x = "Week", y = "No. plants", fill = "Status report",
       title = "Open flowers, all years combined") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85, 0.2))
# Not really.....

# Would it make a difference if we aggregated bi-weekly?
  flowers2wk <- flowers %>%
    mutate(wk2 = ifelse(wk%%2 == 0, wk - 1, wk)) %>%
    arrange(id, yr, wk2, desc(status), desc(midpoint)) %>%
    distinct(id, yr, wk2, .keep_all = TRUE)
  open2wk <- open %>%
    mutate(wk2 = ifelse(wk%%2 == 0, wk - 1, wk)) %>%
    arrange(id, yr, wk2, desc(status), desc(midpoint)) %>%
    distinct(id, yr, wk2, .keep_all = TRUE)
  flowers2wk <- flowers2wk %>%
    group_by(spp, wk2) %>%
    summarize(nplants = n(),
              nyes = sum(status),
              prop = nyes/nplants,
              .groups = "drop") %>%
    data.frame()
  open2wk <- open2wk %>%
    group_by(spp, wk2) %>%
    summarize(nplants = n(),
              nyes = sum(status),
              prop = nyes/nplants,
              .groups = "drop") %>%
    data.frame()
  summary(flowerswk$nplants)
  summary(flowers2wk$nplants)
  summary(openwk$nplants)
  summary(open2wk$nplants)
  
  statusfills <- c("No" = "gray", "Yes" = "steelblue3")
  flowers2wk %>%
    ggplot(aes(x = wk2, y = nplants)) +
    geom_bar(stat = "identity", width = 0.75, aes(fill = "No")) +
    geom_bar(aes(y = nyes, fill = "Yes"), stat = "identity", width = 0.75) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "coral3", linewidth = 0.8) +
    geom_hline(yintercept = 10, linetype = "dotted", color = "coral3", linewidth = 1) +
    facet_wrap(~spp) +
    scale_fill_manual(values = statusfills) +
    labs(x = "Week", y = "No. plants", fill = "Status report",
         title = "Flowers, all years combined") +
    theme_bw() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.85, 0.2))
  
  open2wk %>%
    ggplot(aes(x = wk2, y = nplants)) +
    geom_bar(stat = "identity", width = 0.75, aes(fill = "No")) +
    geom_bar(aes(y = nyes, fill = "Yes"), stat = "identity", width = 0.75) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "coral3", linewidth = 0.8) +
    geom_hline(yintercept = 10, linetype = "dotted", color = "coral3", linewidth = 1) +
    facet_wrap(~spp) +
    scale_fill_manual(values = statusfills) +
    labs(x = "Week", y = "No. plants", fill = "Status report",
         title = "Open flowers, all years combined") +
    theme_bw() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.85, 0.2))
# This does help a little bit actually...
  
# Species with largest bi-weekly sample sizes?
# American beautyberry (more years)
# eastern redbud (more years)
# ~wax mallow
# ~eastern purple coneflower
  
# Other ways to visualize these data (offer users an option?)
  # Heatmap 
  # Bubble figure with various options for average/smoothed curve

# Calculate summary statistics (estimated "onset/offset/duration/peak" based
# on these visualizations/analyses)? 

# Combine phenophase into to get open flower counts? --------------------------#

