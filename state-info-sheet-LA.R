################################################################################
# Create material for state info sheet: LA

# Much of code is derived from that originally developed by Alyssa Rosemartin, 
# Hayley Limes, and Jeff Oliver. 
# See: https://github.com/alyssarosemartin/time-to-restore

# Erin Zylstra
# 2024-06-25
################################################################################

require(dplyr)
require(lubridate)
require(stringr)
require(tidyr)
require(ggplot2)
library(cowplot)
library(terra)
library(tidyterra)
library(mgcv)

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

# Rename/remove columns where necessary and obs with open flower status = NA
df <- df %>%
  filter(!is.na(status_fo)) %>%
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

# Create figure(s) with weekly proportion of plants with open flowers ---------#
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
# < 5 for every year, but 10.2 across all years, so I think we only have that
# option. Will need to automate this in some way......

# Set threshold for mean number of plants observed each week in weeks 10-40
weekly_nobs_min <- 8

sppyears <- spp_nobs %>% 
  filter(nobs_weeklymn_sc >= weekly_nobs_min)
if (nrow(sppyears) == 0) {
  warning("Data are insufficient to characterize open flower phenophase for ", 
          spp)
  # Add a next here, so we skip the rest of the species loop
} else{
  if (nrow(sppyears) == 1) {
    message("Can only characterize open flower phenophase for ", spp, 
            " by combining data across all years")
  } else {
    message("Data are sufficient to characterize open flower phenophase for " ,
            spp, " in ", paste(sppyears$year[-nrow(sppyears)], collapse = ", "))
  }
}

if (nrow(sppyears) > 1) {
  # Create figures for each year with sufficient data
}

# Format data for figures with all years combined
prop_allyrs <- sppdat %>%
  group_by(sc, wk, wk_doy) %>%
  summarize(n_obs = n(),
            n_open = sum(status_fo),
            .groups = "keep") %>%
  data.frame() %>%
  mutate(prop_open = n_open / n_obs) %>%
  mutate(sc_f = as.factor(ifelse(sc == 0, "Other", "SC")))

# Compare proportions for plants in/out of SC region
prop_plot <- ggplot(data = prop_allyrs,
                    aes(x = wk_doy, y = prop_open, group = sc_f)) + 
  geom_point(aes(color = sc_f, size = n_obs), alpha = alphaline) +
  geom_line(aes(color = sc_f), alpha = alphaline) +
  scale_color_manual(values = c("blue", "red")) +
  labs(y = paste0("Proportion of plants with open flowers"), 
       x = "Day of year") +
  annotate("text", x = 365, y = 0.98, label = spp,
           hjust = 1, vjust = 1, fontface = 2) +
  labs(color = "Region", size = "No. obs") +
  theme(text = element_text(size = 10),
        legend.text = element_text(size = 8))
prop_plot
  # Sample sizes a lot smaller, but it does look like plants in SC region might
  # flower earlier and have more of a secondary peak in fall

# Focusing on SC region only
propSC_allyrs <- prop_allyrs %>%
  filter(sc == 1)

# GAM
gam1 <- gam(prop_open ~ s(wk_doy, bs = "cc", k = 20), weights = n_obs,
            data = propSC_allyrs, method = "REML", family = "binomial")
summary(gam1)
# Make predictions
gam1_preds <- data.frame(wk_doy = 1:365)
gam1_preds <- cbind(gam1_preds,
                    as.data.frame(predict(gam1,
                                          newdata = gam1_preds,
                                          type = "link",
                                          se.fit = TRUE))) %>%
  rename(fitl = fit) %>%
  mutate(lcll = fitl - 1.96 * se.fit,
         ucll = fitl + 1.96 * se.fit,
         fit = exp(fitl) / (1 + exp(fitl)),
         lcl = exp(lcll) / (1 + exp(lcll)),
         ucl = exp(ucll) / (1 + exp(ucll))) %>%
  select(-c(lcll, ucll))

# Figure with GAM predictions on top of raw proportions
gam1_plot <- ggplot(data = propSC_allyrs, aes(x = wk_doy, y = prop_open)) + 
  geom_point(aes(size = n_obs), alpha = alphaline, color = "coral1") +
  geom_line(alpha = alphaline, color = "coral1") +
  geom_ribbon(data = gam1_preds, aes(x = wk_doy, y = fit, ymin = lcl, ymax = ucl), 
              color = "gray", alpha = alphapoly, linetype = 0) +
  geom_line(data = gam1_preds, aes(x = wk_doy, y = fit)) +
  labs(y = paste0("Proportion of plants with open flowers"), 
       x = "Day of year",
       size = "No. obs") +
  annotate("text", x = 365, y = 0.98, label =  spp,
           hjust = 1, vjust = 1, fontface = 2) +
  theme(text = element_text(size = 10),
        legend.text = element_text(size = 8))
gam1_plot

# Heat map
propSC_allyrs <- propSC_allyrs %>%
  mutate(tile_ht = 0.5,
         wk_doy_c = wk_doy + 3)


g1 <- ggplot(propSC_allyrs) +
  geom_tile(aes(x = wk_doy_c, y = tile_ht, fill = prop_open)) +
  scale_fill_gradient(low = "gray90", high = "blue") +
  labs(x = "Day of year", fill = "Proportion \nopen") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.right = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.position = "left",
        legend.margin = margin(c(5, 0, 5, 1)))

g2 <- ggplot(propSC_allyrs) +
  geom_col(aes(x = wk_doy_c, y = n_obs), fill = "gray70") +
  geom_hline(yintercept = weekly_nobs_min, 
             linetype = "dashed", color = "steelblue") +
  labs(x = "Day of year", y = "No. observations") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme(axis.line = element_line(color = "black"),
        axis.line.y.right = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 10))

plot_grid(g1, g2, nrow = 2, rel_heights = c(2, 1), align = "v")

# TODO: standardize colors (SC/other) across figure types
# Change y-axis labels from doy to Month?
# Explore text size given ultimate size of figure on info sheet
