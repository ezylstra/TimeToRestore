################################################################################
# Create material for state info sheet: LA

# Much of code is derived from that originally developed by Alyssa Rosemartin, 
# Hayley Limes, and Jeff Oliver. 
# See: https://github.com/alyssarosemartin/time-to-restore

# Erin Zylstra
# 2024-07-08
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
alphaline <- 0.6 # 0.3
alphapoly <- 0.4
figw <- 6.5
figh <- 3
whitebox_palette <- "soft" # muted, atlas
# See: https://dieghernan.github.io/tidyterra/reference/scale_whitebox.html
cols2 <- whitebox.colors(n = 2, palette = whitebox_palette)

# Load csv with sample sizes for weekly proportion of plants with open flowers #
ofss <- read.csv("data/openflower_weeklyobs_samplesizes.csv")

# Identify priority species with sufficient data ------------------------------#
# Selecting priority species for LA that average 6 or more observations in the
# southcentral region per week (for all years combined)
species <- ofss %>%
  filter(LA == 1) %>%
  filter(mn_wkobs_sc >= 6) %>%
  select(!contains("z789"))
  # Note: no restrictions based on number of plant-years

# Load iNat data --------------------------------------------------------------#
# Save for later. May need to do a new download. Want to grab flowering tag if
# possible.

# Load shapefile with state boundaries ----------------------------------------# 
states <- vect("data/states/cb_2017_us_state_500k.shp")
states <- subset(states, 
                 !states$STUSPS %in% c("HI", "AK", "VI", "MP", "GU", "PR", "AS"))
states$sc <- as.factor(ifelse(states$STUSPS %in% c("LA", "NM", "OK", "TX"), 1, 0))

# Load processed NPN status/intensity data and format -------------------------#
df <- read.csv("data/flower-status-intensities-priorityspp.csv")

# Rename/remove columns where necessary and remove obs with open flower 
# status = NA
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
  # We'll be creating wk_doy columns to assign each week with a day of the year.
  # wk_doy1 = start of each week (eg, date for week 1 would be Jan 1)
  # wk_doy4 = middle of each week (eg, date for week 1 would be Jan 4)
  df <- df %>%
    mutate(wk = week(observation_date)) %>%
    # Remove observations in week 53
    filter(wk < 53) %>%
    # Create wk_doy columns
    mutate(wk_date1 = parse_date_time(paste(2024, wk, 1, sep = "/"), "Y/W/w"),
           wk_date1 =  as.Date(wk_date1),
           wk_doy1 = yday(wk_date1),
           wk_date4 = parse_date_time(paste(2024, wk, 4, sep = "/"), "Y/W/w"),
           wk_date4 =  as.Date(wk_date4),
           wk_doy4 = yday(wk_date4))

# Just keep one observation of each plant, each week. Sort so the most
# advanced phenophase gets kept (if more than one value in a week)
df1 <- df %>%
  arrange(common_name, plant_id, year, wk, 
          desc(status_fl), desc(status_fo)) %>%
  distinct(plant_id, year, wk, .keep_all = TRUE) %>%
  # Add an indicator for plants southcentral states
  mutate(sc = 1 * state %in% c("LA", "NM", "OK", "TX"))

# Create figure(s) with weekly proportion of plants with open flowers ---------#
# Eventually, we'll want to automate this and cycle through species, but for
# now try going through the process for 1-2 species

i = 4

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
ggplot(data = states) +
  geom_spatvector(aes(fill = sc), color = "gray30") +
  scale_fill_discrete(type = c("white", "gray80"), guide = "none") +
  geom_spatvector(data = sppsitesv, aes(color = sc), size = 1) +
  scale_color_discrete(type = c("steelblue3", "red"), guide = "none")

# Calculate how much data available in region, by year and across years
# (focus on weeks 10-40)
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
# Similar for eastern baccharis
# Similar for trumpet honeysuckle, except there we have even less data (average
  # 6.1 plants observed each week, across years)

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
  group_by(sc, wk, wk_date1, wk_doy1, wk_date4, wk_doy4) %>%
  summarize(n_obs = n(),
            n_open = sum(status_fo),
            .groups = "keep") %>%
  data.frame() %>%
  mutate(prop_open = n_open / n_obs) %>%
  mutate(sc_f = as.factor(ifelse(sc == 0, "Other", "SC")))

# To create a figure with ticks between month labels on the x axis, need to do a 
# little extra work. 
  # Dates where we want month labels (15th of month)
  x_lab <- as.Date(paste0("2024-", 1:12, "-15"))
  # Dates where we want ticks (1st of month)
  x_tick <- as.Date(c(paste0("2024-", 1:12, "-01"), "2025-01-01"))
  n_x_tick <- length(x_tick)
  # Will specify axis breaks & ticks at 1st and 15th of month. Make labels on
  # the 1st black and change color of tick marks on the 15th to NA.

# Compare proportions for plants in/out of SC region
prop_plot <- ggplot(data = prop_allyrs,
                    aes(x = wk_date4, y = prop_open, group = sc_f)) + 
  geom_point(aes(color = sc_f, size = n_obs), alpha = alphaline) +
  geom_line(aes(color = sc_f), alpha = alphaline) +
  scale_color_whitebox_d(palette = whitebox_palette) +
  scale_x_continuous(breaks = c(x_lab, x_tick),
                     labels = c(month.abb, rep("", n_x_tick))) +
  labs(y = paste0("Proportion of plants with open flowers")) +
  annotate("text", x = as.Date("2024-12-31"), y = 0.98, label = spp,
           hjust = 1, vjust = 1, fontface = 2) +
  labs(color = "Region", size = "No. obs") +
  theme(text = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.ticks.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                              rep("black", n_x_tick))),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                    rep("white", n_x_tick))),
        axis.title.x = element_blank())
prop_plot  
# Does open flower phenology look different in SC region than in the rest of
# the US?

# Focusing on SC region only
propSC_allyrs <- prop_allyrs %>%
  filter(sc == 1)

# GAM
# Looking at the red maple data illustrated a weird problem that can occur. If
# k is set relatively high then sometimes you get a very wide credible interval
# around predictions for a portion of the year when the curve is flat (no plants
# flowering). In these cases, better to set k much lower, and then check that 
# model fit is adequate.
k_values <- seq(5, 20, by = 5)

for (k_index in 1:4) {
  gam1 <- gam(prop_open ~ s(wk_doy4, bs = "cc", k = k_values[k_index]), weights = n_obs,
              data = propSC_allyrs, method = "REML", family = "binomial")
  k_p <- k.check(gam1)[,"p-value"] 
  if(k_p > 0.1) break()
}
if (k_p <= 0.1 & k_index == 4) {
  warning("k value of 20 still results in poor model fit.")
} 

summary(gam1)
message("k set at ", k_values[k_index])

# Make predictions
gam1_preds <- data.frame(
  wk_doy4 = min(propSC_allyrs$wk_doy4):max(propSC_allyrs$wk_doy4),
  wk_date4 = as.Date(min(propSC_allyrs$wk_date4):max(propSC_allyrs$wk_date4)))
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
gam1_plot <- ggplot(data = propSC_allyrs, 
                    aes(x = wk_date4, y = prop_open)) + 
  geom_point(aes(size = n_obs), alpha = alphaline, color = cols2[2]) +
  geom_line(alpha = alphaline, color = cols2[2]) +
  geom_ribbon(data = gam1_preds, 
              aes(x = wk_date4, y = fit, ymin = lcl, ymax = ucl), 
              color = "gray", alpha = alphapoly, linetype = 0) +
  geom_line(data = gam1_preds, aes(x = wk_date4, y = fit)) +
  scale_x_continuous(breaks = c(x_lab, x_tick),
                     labels = c(month.abb, rep("", n_x_tick))) +
  labs(y = paste0("Proportion of plants with open flowers"), 
       size = "No. obs") +
  annotate("text", x = as.Date("2024-12-31"), y = 0.98, label =  spp,
           hjust = 1, vjust = 1, fontface = 2) +
  theme(text = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.ticks.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                              rep("black", n_x_tick))),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                    rep("white", n_x_tick))),
        axis.title.x = element_blank())
gam1_plot

# Heat map
propSC_allyrs <- propSC_allyrs %>%
  mutate(tile_ht = 0.5)

g1 <- ggplot(propSC_allyrs) +
  geom_tile(aes(x = wk_date4, y = tile_ht, fill = prop_open)) +
  scale_fill_gradient(low = "gray95", high = cols2[2]) +
  labs(fill = "Proportion \nopen") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = c(x_lab, x_tick),
                     labels = c(month.abb, rep("", n_x_tick)),
                     expand = c(0.01, 0.01)) +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.right = element_line(color = "black"),
        axis.ticks.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                              rep("black", n_x_tick))),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.position = "left",
        legend.margin = margin(c(5, 0, 5, 1)))

g2 <- ggplot(propSC_allyrs) +
  geom_col(aes(x = wk_date4, y = n_obs), fill = "gray70") +
  geom_hline(yintercept = weekly_nobs_min, 
             linetype = "dashed", color = "black") +
  labs(y = "No. observations") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = c(x_lab, x_tick),
                     labels = c(month.abb, rep("", n_x_tick)),
                     expand = c(0.01, 0.01)) +
  theme(axis.line = element_line(color = "black"),
        axis.ticks.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                              rep("black", n_x_tick))),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 10))

plot_grid(g1, g2, nrow = 2, rel_heights = c(2, 1), align = "v")

