################################################################################
# Explore data on open flower phenophase

# Much of code is derived from that originally developed by Alyssa Rosemartin, 
# Hayley Limes, and Jeff Oliver. 
# See: https://github.com/alyssarosemartin/time-to-restore

# Erin Zylstra
# 2024-06-20
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

# Load processed NPN status/intensity data, species info ----------------------#
df <- read.csv("data/flower-status-intensities-priorityspp.csv")

# Summarize the amount of information available for each plant, year ----------#
samples <- df %>%
  group_by(common_name, individual_id, site_id, state, year) %>%
  summarize(n_obs = n(),                      # No. of daily observations
            n_status_fl = sum(!is.na(status_fl)),  # No. of flower status obs
            n_status_fo = sum(!is.na(status_fo)),  # No. of open flower status obs
            n_value_fl = sum(!is.na(midpoint_fl)), # No. of flower intensity values
            n_value_fo = sum(!is.na(midpoint_fo)), # No. of open flower intensity values
            .groups = "keep") %>%
  data.frame()

# Add species info to samples dataframe (priority level and for each state, an
# indicator of whether that species ever received one vote or more)
spp <- read.csv("data/priority-species.csv")
spp <- spp %>%
  mutate(across(starts_with("votes_"), ~ replace_na(., 0))) %>%
  mutate(across(starts_with("revisit_"), ~ replace_na(., 0))) %>%
  mutate(LA = if_else(votes_LA > 0 | revisit_votes_LA > 0, 1, 0),
         NM = if_else(votes_NM > 0 | revisit_votes_NM > 0, 1, 0),
         OK = if_else(votes_OK > 0 | revisit_votes_OK > 0, 1, 0),
         TX = if_else(revisit_votes_TX > 0, 1, 0)) %>%
  select(common_name, priority, LA, NM, OK, TX)
samples <- samples %>%
  left_join(spp, by = "common_name")

# Summarize data available on open flower phenophase --------------------------#

# Remove any plant year with no open flower status recorded
samples <- samples %>% 
  filter(n_status_fo > 0) %>%
  mutate(ind_yr = paste(individual_id, year, sep = "_"))
df <- df %>%
  mutate(ind_yr = paste(individual_id, year, sep = "_")) %>%
  filter(ind_yr %in% samples$ind_yr)

# How many plants, plant-years with >=1 open flower status, by species?
spp_ss <- samples %>%
  group_by(common_name, priority, LA, NM, OK, TX) %>%
  summarize(n_plantyrs = n(), 
            n_plants = length(unique(individual_id)),
            n_plants_LA = length(unique(individual_id[state == "LA"])),
            n_plants_NM = length(unique(individual_id[state == "NM"])),
            n_plants_OK = length(unique(individual_id[state == "OK"])),
            n_plants_TX = length(unique(individual_id[state == "TX"])),
            .groups = "keep") %>%
  arrange(desc(n_plantyrs)) %>%  
  data.frame()
# LA
spp_ss %>% 
  filter(n_plants_LA > 0) %>%
  arrange(desc(n_plants_LA))
# NM
spp_ss %>% 
  filter(n_plants_NM > 0) %>%
  arrange(desc(n_plants_NM))
# OK
spp_ss %>% 
  filter(n_plants_OK > 0) %>%
  arrange(desc(n_plants_OK))
# TX
spp_ss %>% 
  filter(n_plants_TX > 0) %>%
  arrange(desc(n_plants_TX))

# Most species with very few monitored plants in the 4 states (and this doesn't
# take into account whether there was just a single observation of the plant
# in a year or more). Exception would be red maple in LA (54 plants). 

# Visualize site locations and look at USDA zones -----------------------------#

sites <- df %>%
  select(site_id, latitude, longitude, state) %>%
  rename(lat = latitude,
         lon = longitude) %>%
  distinct()
sitesv <- vect(sites, geom = c("lon", "lat"), crs = "epsg:4269")

# Load polygon layer with USDA hardiness/growing zones (2023)
zones <- vect("data/zones/phzm_us_zones_shp_2023.shp")
zones <- zones[, c("zone", "trange")]

states <- vect("data/states/cb_2017_us_state_500k.shp")
states <- subset(states, 
                 !states$STUSPS %in% c("HI", "AK", "VI", "MP", "GU", "PR", "AS"))

# Visualize zones (specifying levels so zones appear in logical order)
# ggplot(zones) + 
#   geom_spatvector(aes(fill = factor(zone, 
#                                     levels = paste0(rep(3:12, each = 2), 
#                                                     rep(c("a", "b"))))), 
#                   col = NA) +
#   scale_fill_whitebox_d(palette = "muted", name = "Zone") +
#   geom_spatvector(data = states, aes(fill = NA)) +
#   geom_spatvector(data = sitesv, size = 0.3)

# What if we just used zones 7-9?
# ggplot(subset(zones, zones$zone %in% c("7a", "7b", "8a", "8b", "9a", "9b"))) + 
#   geom_spatvector(aes(fill = factor(zone)), col = NA) +
#   scale_fill_whitebox_d(palette = "muted", name = "Zone") +
#   geom_spatvector(data = states, aes(fill = NA)) +
#   geom_spatvector(data = sitesv, size = 0.3)

# Likely want to restrict things by latitude (38-deg N?)
zext38 <- ext(zones)
zext38[4] <- 38
zones38 <- crop(zones, zext38)

# ggplot() +
#   geom_spatvector(data = states, aes(fill = NA)) +
#   geom_spatvector(data = subset(zones38, 
#                          zones38$zone %in% c("7a", "7b", "8a", 
#                                              "8b", "9a", "9b")),
#                   aes(fill = factor(zone)), col = NA) +
#   scale_fill_whitebox_d(palette = "muted", name = "Zone") +
#   geom_spatvector(data = states, aes(fill = NA)) +
#   geom_spatvector(data = sitesv, size = 0.3)

# Extract zone associated with each site
e <- terra::intersect(sitesv, zones)
sites <- as.data.frame(e, geom = "XY") %>%
  rename(lon = x, lat = y)
# Add zone info to observation dataframe
df <- df %>%
  select(-c(latitude, longitude)) %>%
  left_join(sites, by = c("site_id", "state")) %>%
  mutate(zone = str_pad(zone, width = 3, side = "left", pad = 0))

# Add week to data, so we can calculate weekly proportions --------------------#

# We'll be creating doy column to assign each week with a day of the year.
day_of_week <- 1
# Note: if we want the generic date and DOY to be the start of each week (eg,
# date for week 1 would be Jan 1), then set day_of_week = 1. If we want the
# generic date and DOY to be mid week (Thursday) then set day_of_week = 4. 

df <- df %>%
  mutate(wk = week(observation_date),
         date_generic = parse_date_time(paste(2024, wk, day_of_week, sep = "/"), 
                                        "Y/W/w"),
         date_generic = as.Date(date_generic),
         doy = yday(date_generic),
         fyear = factor(year, levels = as.character(2013:2023)))

# Summarize weekly data available, by species ---------------------------------#

# Just keep one observation of each individual, each week. Sort so the most
# advanced phenophase gets kept (if more than one value in a week)
df1 <- df %>%
  arrange(common_name, individual_id, year, wk, 
          desc(status_fl), desc(status_fo)) %>%
  distinct(individual_id, year, wk, .keep_all = TRUE) %>%
  # Remove week 53
  filter(wk != 53) %>%
  mutate(z789 = ifelse(lat <=38 & grepl("7|8|9", zone), 1, 0),
         sc = 1 * state %in% c("LA", "NM", "OK", "TX"))

# Summarize across all plants, years
propdat_wk <- df1 %>%
  group_by(common_name, wk) %>%
  summarize(n_obs = n(), .groups = "keep") %>%
  data.frame()
  # Calculate mean no. observations used to calculate proportion open in weeks 
  # 10-40 (~ March - September)
propdat <- propdat_wk %>%
  filter(wk %in% 10:40) %>%
  group_by(common_name) %>%
  summarize(mn_wkobs = round(mean(n_obs), 1)) %>%
  data.frame()

# Summarize across plants in 4 states, all years
propdat_wk_sc <- df1 %>%
  filter(sc == 1) %>%
  group_by(common_name, wk) %>%
  summarize(n_obs_sc = n(), .groups = "keep") %>%
  data.frame()
propdat_sc <- propdat_wk_sc %>%
  filter(wk %in% 10:40) %>%
  group_by(common_name) %>%
  summarize(mn_wkobs_sc = round(mean(n_obs_sc), 1)) %>%
  data.frame()

# Summarize across plants in zones 7-9 and below 38-deg lat, all years
propdat_wk_z789 <- df1 %>%
  filter(z789 == 1) %>%
  group_by(common_name, wk) %>%
  summarize(n_obs_z789 = n(), .groups = "keep") %>%
  data.frame()
propdat_z789 <- propdat_wk_z789 %>%
  filter(wk %in% 10:40) %>%
  group_by(common_name) %>%
  summarize(mn_wkobs_z789 = round(mean(n_obs_z789), 1)) %>%
  data.frame()

# Get number of plants, plant years
nplants <- df1 %>%
  group_by(common_name) %>%
  summarize(n_plants = length(unique(individual_id)),
            n_plants_sc = length(unique(individual_id[sc == 1])),
            n_plants_z789 = length(unique(individual_id[z789 == 1])),
            n_plantyrs = length(unique(ind_yr)),
            n_plantyrs_sc = length(unique(ind_yr[sc == 1])),
            n_plantyrs_z789 = length(unique(ind_yr[z789 == 1]))) %>%
  data.frame()

# Combine everything
propdata <- spp %>%
  right_join(nplants, by = "common_name") %>%
  left_join(propdat, by = "common_name") %>%
  left_join(propdat_sc, by = "common_name") %>%
  left_join(propdat_z789, by = "common_name") %>%
  arrange(priority, desc(n_plantyrs_sc))

# Look at data summaries for all plants, plants in 4 southcentral states
propdata %>% select(!contains("z789"))

# Write to file
# write.csv(propdata, 
#           "data/openflower_weeklyobs_samplesizes.csv", 
#           row.names = FALSE)

# Explore wild bergamot data --------------------------------------------------#
# Priority 1 species. 107 plants monitored, 207 plant-years.

wb <- df1 %>% filter(common_name == "wild bergamot")

# Where are the plants?
wb_sites <- wb %>%
  select(site_id, lat, lon, state, zone) %>%
  distinct()
wb_sitesv <- vect(wb_sites, geom = c("lon", "lat"), crs = "epsg:4269")
  # 2 sites in LA (zones 9a, 9b); 1 in NM (zone 6b)

ggplot(subset(zones, zones$zone %in% c("7a", "7b", "8a", "8b", "9a", "9b"))) + 
  geom_spatvector(aes(fill = factor(zone)), col = NA) +
  scale_fill_whitebox_d(palette = "muted", name = "Zone") +
  geom_spatvector(data = states, aes(fill = NA)) +
  geom_spatvector(data = wb_sitesv, size = 0.5)

# Limit which years are included by the number of individuals observed?
  # Just use years when 10 or more individuals were observed?
  # Just use years when the average number of individuals observered per week 
    # between weeks 20 and 40 is >= 10
wb_nobs <- wb %>%
  group_by(year) %>%
  summarize(nplants = length(unique(individual_id))) %>%
  data.frame()
wb_nobsperwk <- wb %>%
  filter(wk %in% 20:40) %>%
  group_by(year, wk) %>%
  summarize(nobs_open = n(), .groups = "keep") %>%
  group_by(year) %>%
  summarize(nobs_weeklymn = mean(nobs_open)) %>%
  data.frame()
wb_nobs <- wb_nobs %>%
  left_join(wb_nobsperwk, by = "year")
wb_nobs

# Calculate proportion of flowers open each week, year (all plants; 2019-2023)
wb_all_yrprop <- wb %>%
  filter(year %in% 2019:2023) %>%
  group_by(year, wk, doy, fyear) %>%
  summarize(nobs = n(),
            nobs_open = sum(!is.na(status_fo)),
            n_open = sum(status_fo, na.rm = TRUE),
            .groups = "keep") %>%
  data.frame() %>%
  mutate(prop_open = n_open / nobs_open)

# Calculate proportion of flowers open each week, across years
wb_all_prop <- wb_all_yrprop %>%
  group_by(wk, doy) %>%
  summarize(nobs_open = sum(nobs_open),
            n_open = sum(n_open),
            .groups = "keep") %>%
  data.frame() %>%
  mutate(prop_open = n_open / nobs_open)
# ggplot object
wb_raw_plot <- ggplot() + 
  geom_point(data = wb_all_prop, aes(x = doy, y = prop_open, size = nobs_open),
             alpha = alphaline) +
  geom_line(data = wb_all_prop, aes(x = doy, y = prop_open), 
            alpha = alphaline) +
  labs(y = paste0("Proportion of plants with open flowers"), 
       x = "Day of year") +
  annotate("text", x = 365, y = 0.98, label = "Wild bergamot",
           hjust = 1, vjust = 1, fontface = 2) +
  theme(text = element_text(size = 10),
        legend.text = element_text(size = 8))
wb_raw_plot

# We can use a GAM to model the proportion of plants with open flowers and
# using (subjective) criteria, identify the start/end/duration of open flower 
# phenophase. Can also use the GAM to identify where there's strong evidence of
# annual variation in timing of open flowers

# GAM with the same, cyclic cubic smooth each year 
# see Canizares et al. 2023, 
# https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
wb_fo1 <- gam(prop_open ~ s(doy, bs = "cc", k = 20), weights = nobs_open,
              data = wb_all_yrprop, method = "REML", family = "binomial")

  # Make predictions
  wb_fo1_preds <- data.frame(doy = 1:365)
  wb_fo1_preds <- cbind(wb_fo1_preds,
                        as.data.frame(predict(wb_fo1,
                                              newdata = wb_fo1_preds,
                                              type = "link",
                                              se.fit = TRUE))) %>%
    rename(fitl = fit) %>%
    mutate(lcll = fitl - 1.96 * se.fit,
           ucll = fitl + 1.96 * se.fit,
           fit = exp(fitl) / (1 + exp(fitl)),
           lcl = exp(lcll) / (1 + exp(lcll)),
           ucl = exp(ucll) / (1 + exp(ucll))) %>%
    select(-c(lcll, ucll))

  # ggplot object
  wb_fo1_plot <- ggplot() + 
    geom_point(data = wb_all_yrprop, 
               aes(x = doy, y = prop_open, group = fyear, color = fyear),
               alpha = alphaline, size = 1) +
    geom_line(data = wb_all_yrprop, 
              aes(x = doy, y = prop_open, group = fyear, color = fyear),
              alpha = alphaline) +
    scale_color_discrete(name = "Year") +    
    geom_ribbon(data = wb_fo1_preds, aes(x = doy, ymin = lcl, ymax = ucl), 
                color = "gray", alpha = alphapoly, linetype = 0) +
    geom_line(data = wb_fo1_preds, aes(x = doy, y = fit)) +
    labs(y = paste0("Proportion of plants with open flowers"), 
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "Wild bergamot",
             hjust = 1, vjust = 1, fontface = 2) +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8))
  wb_fo1_plot

# GAM with an independent smooth each year
wb_foy <- gam(prop_open ~ fyear + s(doy, by = fyear, bs = "cc", k = 20), 
              weights = nobs_open, data = wb_all_yrprop, method = "REML", 
              family = "binomial")
# Compare the two models with AIC
AIC(wb_fo1, wb_foy)

  # Make predictions
  wb_yrpreds <- data.frame(fyear = NA)
  for(yr in sort(unique(wb_all_yrprop$year))) {
    min_doy <- min(wb_all_yrprop$doy[wb_all_yrprop$year == yr])
    max_doy <- max(wb_all_yrprop$doy[wb_all_yrprop$year == yr])
    max_doy <- ifelse(max_doy > 350, 365, max_doy)
    doys <- seq(min_doy, max_doy)
    if (yr == min(wb_all_yrprop$year)) {
      wb_yrpreds <- data.frame(fyear = yr, doy = doys)
    } else {
      wb_yrpreds <- rbind(wb_yrpreds, data.frame(fyear = yr, doy = doys))
    }
  }
  wb_yrpreds$fyear <- factor(wb_yrpreds$fyear, levels = as.character(2012:2023))
  wb_foy_preds <- cbind(wb_yrpreds,
                        as.data.frame(predict(wb_foy, 
                                              newdata = wb_yrpreds,
                                              type = "link", 
                                              se.fit = TRUE))) %>%
    rename(fitl = fit) %>%
    mutate(lcll = fitl - 1.96 * se.fit,
           ucll = fitl + 1.96 * se.fit,
           fit = exp(fitl) / (1 + exp(fitl)),
           lcl = exp(lcll) / (1 + exp(lcll)),
           ucl = exp(ucll) / (1 + exp(ucll))) %>%
    select(-c(lcll, ucll))
  
  # ggplot object
  wb_foy_plot <- ggplot(data = wb_foy_preds, 
                        aes(x = doy, y = fit, group = fyear, color = fyear)) +
    scale_color_discrete(name = "Year") +
    geom_line() +
    labs(y = paste0("Proportion of plants with open flowers"), 
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "Wild bergamot",
             hjust = 1, vjust = 1, fontface = 2) +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8))
  wb_foy_plot

# There's not a lot of sites, but we could see whether there's evidence 
# that open flower phenology differs for plants in zones 7-9 verus plants in 
# other zones
  
  # Calculate proportion of flowers open each week, region (calculating
  # proportion across years [2019:2023] because otherwise sample sizes are 
  # really small in most year-wk-region combinations)
  wb_zone_prop <- wb %>%
    mutate(zone789 = ifelse(grepl("7|8|9", zone), 1, 0)) %>%
    mutate(zone789 = factor(zone789)) %>%
    filter(year %in% 2019:2023) %>%
    group_by(wk, doy, zone789) %>%
    summarize(nobs = n(),
              nobs_open = sum(!is.na(status_fo)),
              prop_open = sum(status_fo, na.rm = TRUE) / nobs_open,
              .groups = "keep") %>%
    data.frame()
  
  wb_foz <- gam(prop_open ~ zone789 + s(doy, by = zone789, bs = "cc", k = 20), 
                weights = nobs_open, data = wb_zone_prop, method = "REML", 
                family = "binomial")
  
  # Make predictions
  wb_foz_preds <- expand.grid(doy = 1:365, 
                              zone789 = factor(0:1), 
                              KEEP.OUT.ATTRS = FALSE)
  wb_foz_preds <- cbind(wb_foz_preds,
                        as.data.frame(predict(wb_foz,
                                              newdata = wb_foz_preds,
                                              type = "link",
                                              se.fit = TRUE))) %>%
    rename(fitl = fit) %>%
    mutate(lcll = fitl - 1.96 * se.fit,
           ucll = fitl + 1.96 * se.fit,
           fit = exp(fitl) / (1 + exp(fitl)),
           lcl = exp(lcll) / (1 + exp(lcll)),
           ucl = exp(ucll) / (1 + exp(ucll))) %>%
    select(-c(lcll, ucll))
  
  # ggplot object
  wb_foz_plot <- ggplot() +
    geom_point(data = wb_zone_prop, 
               aes(x = doy, y = prop_open, group = zone789, color = zone789,
                   size = nobs_open), alpha = alphaline) +
    geom_line(data = wb_zone_prop, 
              aes(x = doy, y = prop_open, group = zone789, color = zone789),
              alpha = alphaline) +
    scale_color_discrete(name = "Zone789") +    
    geom_line(data = wb_foz_preds, aes(x = doy, y = fit, group = zone789,
                                       color = zone789), linewidth = 1.5) +
    labs(y = paste0("Proportion of plants with open flowers"), 
         x = "Day of year") +
    annotate("text", x = 365, y = 0.98, label = "Wild bergamot",
             hjust = 1, vjust = 1, fontface = 2) +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8))
  wb_foz_plot
  
  # Can't use AIC to compare models since the zone model is based on a different
  # dataset...
  # But, directly comparing models with regional curves and a single curve (both
  # across years), suggests regional curves fit the data much better. However, 
  # it does look like there's a bunch of noise in the zone 7-9 data, probably
  # because sample sizes are small. 
  wb_foz1 <- gam(prop_open ~ s(doy, bs = "cc", k = 20), 
                 weights = nobs_open, data = wb_zone_prop, method = "REML", 
                 family = "binomial")
  AIC(wb_foz, wb_foz1)

  wb_zone_prop %>% arrange(zone789, wk)
  count(filter(wb_zone_prop, zone789 == 1), nobs)
    # Most weeks (40/52) had 1-5 plant observations. The most obs in any week was 8
  
# Create "heatmaps" for proportion of plants with open flowers ----------------#
# Continue to use wild bergamot for an example. Using the dataset create above
# that contains just one observation of a plant per week
  
# Need to decide how we calculate these proportions....
  # All plants or just those within some region that includes 4 target states?
  # All years combined, a subset of years combined, or individually by year?
  
  # Look at data available by year for wild bergamot:
  wb_nobs
    # Won't use data from 2013-2015, when only 1-2 plants were observed.
  wb <- filter(wb, year %in% 2016:2023)  
  
  # Look where plants are located
  wb_plants <- wb %>%
    select(individual_id, state, site_id, lat, lon, zone) %>%
    distinct() 
  count(wb_plants, grepl("7|8|9", zone))
  count(wb, grepl("7|8|9", zone))
    # Only 23 of 104 plants (22%) and 11% of observations in zones 7-9, 
    # so might be difficult to split things out by region (and the small sample
    # sizes in zone 7-9 would be apparent in a heat map)
  
  # Calculating proportions across all plants, across years (2016-2023)
  wbprop <- wb %>%
    filter(!is.na(status_fo)) %>%
    group_by(wk, doy, date_generic) %>%
    summarize(n_obs = n(),
              sum_open = sum(status_fo),
              prop_open = sum_open / n_obs,
              .groups = "keep") %>%
    mutate(y = 1) %>%
    data.frame()
  wbprop
    # Weekly sample sizes range from 5-83
  
  # Quick look (using geom_tile code from original repo, but will need to spend
  # more time looking into how the call works. Not sure that the limits, labels
  # are all correct. Before I expanded the limits a little I was getting a 
  # warning that some rows were removed)
  g1 <- ggplot(wbprop) +
    geom_line(aes(x = wk, y = prop_open))
  g2 <- ggplot(wbprop) +
    geom_tile(aes(x = wk, y = y, fill = prop_open)) +
    scale_fill_gradient(low = "gray90", high = "blue", (name = "Proportion")) +
    scale_x_continuous(breaks = seq(1, 52, by = 4.25), limits=c(0.5, 52.5),
                       labels = c("Jan", "Feb", "Mar", "Apr", "May", "June",
                                "July", "Aug", "Sep", "Oct", "Nov", "Dec", "")) +
    theme(legend.position = "none")
  grid.arrange(g1, g2, ncol = 1)

  # Calculating proportions across all plants, each year (2016-2023)
  wbprop_yr <- wb %>%
    filter(!is.na(status_fo)) %>%
    group_by(year, wk, doy, date_generic) %>%
    summarize(n_obs = n(),
              sum_open = sum(status_fo),
              prop_open = sum_open / n_obs,
              .groups = "keep") %>%
    mutate(y = 1) %>%
    data.frame()
  head(wbprop_yr, 20)
  count(wbprop_yr, n_obs)
  # Weekly sample sizes range from 1-21 (median = 6)
  # Can't really do a lot with that....
  summary(filter(wbprop_yr, wk %in% 20:40))
  count(filter(wbprop_yr, wk %in% 20:40), n_obs)
  # For weeks 20-40, sample sizes range from 1-19 (median = 9)
  
  