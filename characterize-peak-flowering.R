################################################################################
# Explore what's possible with number of open flower estimates

# Much of code is derived from that originally developed by Alyssa Rosemartin, 
# Hayley Limes, and Jeff Oliver. 
# See: https://github.com/alyssarosemartin/time-to-restore

# Erin Zylstra
# 2024-05-30
################################################################################

require(dplyr)
require(lubridate)
require(tidyr)
require(ggplot2)

rm(list = ls())

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
            n_ests_open = sum(!is.na(num_open)),   # No. of observations with
                                                      # an estimated number of 
                                                      # open flowers
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

# Summarize data available on "peak" flowering --------------------------------#
# Based on estimated number of open flowers from status-intensity data

# How many plant-years have >1 estimate of the number of open flowers, by 
# species?
spp_ss <- count(filter(samples, n_ests_open > 1), 
                common_name, priority, LA, NM, OK, TX) %>%
  arrange(desc(n))
spp_ss
sum(spp_ss$n)
  # Range = 2-1549 (2-295 if you exclude red maple); median = 20. 
  # 15 species that have >=30 plant-years with 2 or more estimates. But note
    # we haven't restricted things geographically.
  # Key plants for LA: red maple, black elderberry, butterfly milkweed, eastern
    # purple coneflower, wild bergamot, common buttonbush, New England aster,
    # cardinalflower, trumpet honeysuckle
  # Key plants for NM: butterfly milkweed, eastern purple coneflower, rubber
    # rabbitbrush, wild bergamot, blackeyed Susan, cardinalflower, common
    # sunflower, horsetail milkweed
  # Key plants for OK: common milkweed, butterfly milkweed, eastern purple 
    # coneflower, wild bergamot, swamp milkweed, blackeyed Susan, common 
    # buttonbush, New England aster, cardinalflower, common sunflower
  # Key plants for TX: butterfly milkweed

# How many plant-years in each of the 4 target states have >1 estimate?
states4 <- c("LA", "NM", "OK", "TX")
count(filter(samples, n_ests_open > 1 & state %in% states4), 
      state, common_name, priority, LA, NM, OK, TX) %>% 
  arrange(state, desc(n))
# No plants in OK
# LA: two spp with >10 plant-years (red maple [105], buttonbush [18])
# NM: two spp with >10 plant-years (rabbitbrush [21], horsetail milkweed [19])
# TX: two spp with >10 plant-years (sunflower [13], trumpet honeysuckle [12])

# Identify when plants in "peak" flowering ------------------------------------#
# Based on estimated number of open flowers

# Remove plant-year combinations where we have <2 (non-zero) estimates of the 
# number of open flowers
peak_df <- df %>%  
  group_by(individual_id, year) %>%
  filter(sum(!is.na(num_open)) > 1) %>%
  # filter(any(num_open != 0)) %>%  # Use this option if keeping plant-yrs with 1 estimate
  ungroup() %>%
  mutate(sc = ifelse(state %in% c("LA", "NM", "OK", "TX"), 1, 0)) %>%
  data.frame()

# Extract the maximum number of open flowers estimated for each plant-year
peak_df <- peak_df %>%
  group_by(individual_id, year) %>%
  mutate(max_open = max(num_open, na.rm = TRUE)) %>%
  ungroup() %>%
  data.frame()

# For each plant, summarizing interannual variation in the maximum number of 
# open flowers
  # maxfo_annualvar <- peak_df %>%
  #   select(common_name, individual_id, year, max_open) %>%
  #   distinct() %>%
  #   group_by(individual_id, common_name) %>%
  #   mutate(n_yrs = n()) %>%
  #   filter(n_yrs > 1) %>%
  #   summarize(n_yrs = n(),
  #             min = min(max_open),
  #             mn = mean(max_open),
  #             max = max(max_open),
  #             range = max - min,
  #             sd = sd(max_open),
  #             cv = sd / mn * 100,
  #             .groups = "keep") %>%
  #   data.frame()
  # maxfo_annualvar %>%
  #   group_by(common_name) %>%
  #   summarize(n_plants = n(),
  #             mean_cv = mean(cv)) %>%
  #   arrange(desc(n_plants)) %>%
  #   data.frame()
# Looks like a lot of interannual variation in the estimated maximum number 
# of open flowers, so probably shouldn't use maximum for each plant across years

# Looking at a few individuals of high priority species
  # filter(maxfo_annualvar, common_name == "butterfly milkweed")
  # peak_df %>%
  #   filter(individual_id == 14304 & fl_s == 1) %>%
  #   select(common_name, observation_date, year, midpoint_fl, midpoint_fo, 
  #          num_open, max_open)
  # filter(maxfo_annualvar, common_name == "eastern purple coneflower")
  # peak_df %>%
  #   filter(individual_id == 66579) %>%
  #   select(common_name, observation_date, year, midpoint_fl, midpoint_fo, 
  #          num_open, max_open)
# How much variation in max number of open flowers within, among species?
  # max_flowers <- peak_df %>%
  #   select(common_name, individual_id, year, max_open) %>%
  #   distinct() %>%
  #   group_by(common_name) %>%
  #   summarize(n_plants = length(unique(individual_id)),
  #             n_yrs = length(unique(year)),
  #             n_plantyrs = n(),
  #             min = min(max_open),
  #             mn = mean(max_open),
  #             md = median(max_open),
  #             max = max(max_open),
  #             sd = sd(max_open)) %>%
  #   data.frame()
  # max_flowers
# A lot of variation among species

# Select criteria we'll use to delineate whether observation is part of "peak"
# 1) estimated number of open flowers within x% of the maximum value for that
#    plant and year,
# 2) criteria above OR estimated number of open flowers > 500.

peak_def <- "near max"
# peak_def <- "near max or over 500"

# Explore "threshold" amount of open flowers (%, relative to the max), 
# above which, will be considered peak flowering.
thresholds <- c(0.50, 0.67, 0.75)

# Classify observations of plant on particular date as in peak or not
peak_df <- peak_df %>%
  group_by(individual_id, year) %>%
  mutate(near_max50 = if_else(num_open >= max_open * thresholds[1], 1, 0)) %>%
  mutate(near_max67 = if_else(num_open >= max_open * thresholds[2], 1, 0)) %>%
  mutate(near_max75 = if_else(num_open >= max_open * thresholds[3], 1, 0)) %>%
  mutate(over_500 = ifelse(num_open > 500, 1, 0)) %>%
  ungroup() %>%
  data.frame()
# Number of plant-years above three thresholds
count(peak_df, near_max50, near_max67, near_max75)
# 50%: 8159
# 67%: 6872
# 75%: 6384

# Will use 67% threshold for now.
if (peak_def == "near max") {
  peak_df <- peak_df %>%
    mutate(peak = ifelse(near_max67 == 1, 1, 0)) %>%
    select(-starts_with("near_")) %>%
    select(-over_500)
} else {
  peak_df <- peak_df %>%
    mutate(peak = ifelse(near_max67 == 1 | over_500 == 1, 1, 0)) %>%
    select(-starts_with("near_")) %>%
    select(-over_500)
}

# Change estimated number of open flowers (num_open) and peak values to 0 when 
# status for flowers or open flowers is "no"/0
peak_df <- peak_df %>%
  mutate(num_open = case_when(
    !is.na(status_fl) & status_fl == 0 ~ 0,
    !is.na(status_fo) & status_fo == 0 ~ 0,
    .default = num_open)) %>%
  mutate(peak = case_when(
    !is.na(status_fl) & status_fl == 0 ~ 0,
    !is.na(status_fo) & status_fo == 0 ~ 0,
    .default = peak))
# checks:
count(peak_df, status_fl, status_fo, num_open)
count(peak_df, status_fl, status_fo, peak)

# Summarize data for each plant-year ------------------------------------------#

# Create a new dataframe with summaries of data related to peak flowering
# for each plant & year
plantyr <- peak_df %>%
  group_by(site_id, common_name, individual_id, sc, year) %>%
  summarize(n_obs = n(),
            n_status_fl = sum(!is.na(status_fl)),
            n_status_fo = sum(!is.na(status_fo)),
            n_ests_open = sum(!is.na(num_open)),
            first_obs = min(day_of_year),
            last_obs = max(day_of_year),
            .groups = "keep") %>%
  data.frame()

plantyr$open_status_interval_mn <- NA
plantyr$open_status_interval_max <- NA
plantyr$open_est_interval_mn <- NA
plantyr$open_est_interval_max <- NA
plantyr$first_open <- NA
plantyr$numdays_since_closed <- NA
plantyr$last_open <- NA
plantyr$numdays_until_next_close <- NA
plantyr$duration_open <- NA
plantyr$discontinuous_open <- NA
plantyr$num_obs_open_period <- NA
plantyr$num_ests_open_period <- NA
plantyr$first_peak <- NA
plantyr$last_peak <- NA
plantyr$duration_peak <- NA
plantyr$discontinuous_peak <- NA

for (i in 1:nrow(plantyr)) {
  
  # Subset of peak_df where flowers open status is known (ie, not NA)
  open_df <- peak_df %>%
    filter(individual_id == plantyr$individual_id[i],
           year == plantyr$year[i],
           !is.na(status_fo))
  open_est_df <- open_df %>% filter(!is.na(num_open))
  
  # Number of days between observations where open flower status was recorded?
  if (nrow(open_df) > 1) {
    plantyr$open_status_interval_mn[i] <- round(mean(diff(open_df$day_of_year)), 1)
    plantyr$open_status_interval_max[i] <- round(max(diff(open_df$day_of_year)), 1)
  }
  # Number of days between estimates of the number of open flowers?
  if (nrow(open_est_df) > 1) {
    plantyr$open_est_interval_mn[i] <- round(mean(diff(open_est_df$day_of_year)), 1)
    plantyr$open_est_interval_max[i] <- round(max(diff(open_est_df$day_of_year)), 1)
  }
  
  first_open_ind <- first(which(open_df$status_fo == 1))
  first_closed_ind <- first(which(open_df$status_fo == 0))
  # First day of year with open flowers
  plantyr$first_open[i] <- open_df$day_of_year[first_open_ind]
  # Number of days since previous observation without open flowers
  if (first_open_ind > first_closed_ind & !is.na(first_closed_ind)) {
    plantyr$numdays_since_closed[i] <- 
      plantyr$first_open[i] - open_df$day_of_year[first_open_ind - 1]
  } else {
    plantyr$numdays_since_closed[i] <- NA
  }
  
  last_open_ind <- last(which(open_df$status_fo == 1))
  last_closed_ind <- last(which(open_df$status_fo == 0))
  # Last day of year with open flowers 
  plantyr$last_open[i] <- open_df$day_of_year[last_open_ind]
  # Number of days until next observation without open flowers
  if (last_open_ind < last_closed_ind & !is.na(last_closed_ind)) {
    plantyr$numdays_until_next_close[i] <- 
      open_df$day_of_year[last_open_ind + 1] - plantyr$last_open[i]
  } else {
    plantyr$numdays_until_next_close[i] <- NA
  }
  # Number of days between first and last open flower observation
  plantyr$duration_open[i] <- plantyr$last_open[i] - plantyr$first_open[i] + 1
  open_period <- open_df[first_open_ind:last_open_ind,]
  
  # Were there observations without open flowers between first/last date?
  plantyr$discontinuous_open[i] <- ifelse(any(open_period$status_fo == 0), 1, 0)
  
  # Number of observations between first/last date with open flowers (inclusive)
  plantyr$num_obs_open_period[i] <- nrow(open_period)
  # Number of open flower estimates between first/last date with open flowers
  plantyr$num_ests_open_period[i] <- sum(open_period$status_fo == 1)
  
  # Same stats, but for peak
  first_peak_ind <- first(which(open_est_df$peak == 1))
  last_peak_ind <- last(which(open_est_df$peak == 1))
  plantyr$first_peak[i] <- open_est_df$day_of_year[first_peak_ind]
  plantyr$last_peak[i] <- open_est_df$day_of_year[last_peak_ind]
  plantyr$duration_peak[i] <- plantyr$last_peak[i] - plantyr$first_peak[i] + 1
  plantyr$discontinuous_peak[i] <- 
    ifelse(any(open_est_df$peak[first_peak_ind:last_peak_ind] == 0), 1, 0)

}

# Summaries of estimated PEAK duration/continuity:
sum(plantyr$duration_peak == 1) / nrow(plantyr) 
  # 44% of "peaks" have duration == 1
count(filter(plantyr, duration_peak > 1), common_name)
  # 861 (50%) of plant-years with duration > 1 are red maple
sum(plantyr$duration_peak > 200)
  # 2 plant-years with duration > 200
sum(plantyr$discontinuous_peak == 1) / nrow(plantyr)
sum(plantyr$discontinuous_peak[plantyr$duration_peak > 1] == 1) / sum(plantyr$duration_peak > 1)
  # About 13% of peaks discontinous (or 23% of peaks with duration > 1 day)
count(plantyr, duration_peak, discontinuous_peak)
summary(glm(discontinuous_peak ~ duration_peak, 
            data = plantyr, family = "binomial"))
  # As durations increase, so does probability of discontinous peak

# Summaries, by species
plantyr_spp <- plantyr %>%
  group_by(common_name) %>%
  summarize(n_plants = length(unique(individual_id)),
            n_plants_sc = length(unique(individual_id[sc == 1])),
            n_plantyrs = n(),
            n_plantyrs_sc = length(individual_id[sc == 1]),
            nobsperyr = round(mean(n_obs), 1),
            nobsperyr_sc = round(mean(n_obs[sc == 1]), 1),
            estinterval = round(mean(open_est_interval_mn), 1),
            estinterval_sc = round(mean(open_est_interval_mn[sc == 1]), 1)) %>%
  data.frame()

openest_spp <- spp %>%
  right_join(plantyr_spp, by = "common_name") %>%
  arrange(priority, desc(n_plantyrs_sc))

# Write to file
# write.csv(openest_spp, 
#           "data/estopenflower_samplesizes.csv",
#           row.names = FALSE)

# Summaries of estimated OPEN FLOWER duration/continuity:
  # Note that the duration of the open flower period HAS to be > 1 given that 
  # we've only include plant-years with >1 estimate of the number of open flowers
sum(plantyr$duration_open > 200)
  # 41 plant-years with duration > 200
sum(plantyr$discontinuous_open == 1) / nrow(plantyr)
  # About 25% of open flower periods discontinous
count(plantyr, duration_open, discontinuous_open)
summary(glm(discontinuous_open ~ duration_open, 
            data = plantyr, family = "binomial"))
  # As durations increase, so does probability of discontinous open flower period

# Narrow plant-years based on amount of data availabile to characterize -------#
# peak flowering --------------------------------------------------------------#

# Use criteria to identify plant-yrs where it would be reasonable (or not) to 
# characterize "peak" flowering based on the estimated number of open flowers

# Would like to see one or more observations with no open flowers prior to open 
# flower period (note: if numdays_since_closed = NA, then no prior observations)
hist(plantyr$numdays_since_closed, breaks = 50)
sum(plantyr$numdays_since_closed <= 14, na.rm = TRUE) / nrow(plantyr)
  # 70% of plant-years have a "no" within 2 weeks of first yes to open
# Similarly for observations after the open flower period
hist(plantyr$numdays_until_next_close, breaks = 50)
sum(plantyr$numdays_until_next_close <= 14, na.rm = TRUE) / nrow(plantyr)
  # 67% of plant years have a "no" within 2 weeks of last yes to open

plantyr2 <- plantyr %>% 
  filter(numdays_since_closed <= 14) %>%
  filter(numdays_until_next_close <= 14)

# Filter by a minimum number of observations 
hist(plantyr2$n_status_fo, breaks = 50)
hist(plantyr2$n_ests_open, breaks = 50)
  # Minimum of 5 observations with estimated number of flowers (can include 0)

plantyr2 <- plantyr2 %>% filter(n_ests_open >=5)

# Want plants that were observed regularly 
hist(plantyr2$open_est_interval_mn, breaks = 50)
hist(plantyr2$open_est_interval_max, breaks = 50)
  # Average 2 weeks or less between consecutive estimates of open flowers

plantyr2 <- plantyr2 %>% filter(open_est_interval_mn <= 14)

# What are we left with?
plantyr2_spp <- plantyr2 %>%
  group_by(common_name) %>%
  summarize(n_plants = length(unique(individual_id)),
            n_plants_sc = length(unique(individual_id[sc == 1])),
            n_plantyrs = n(),
            n_plantyrs_sc = length(individual_id[sc == 1]),
            nobsperyr = round(mean(n_obs), 1),
            nobsperyr_sc = round(mean(n_obs[sc == 1]), 1),
            estinterval = round(mean(open_est_interval_mn), 1),
            estinterval_sc = round(mean(open_est_interval_mn[sc == 1]), 1)) %>%
  data.frame()

openest2_spp <- spp %>%
  right_join(plantyr2_spp, by = "common_name") %>%
  arrange(priority, desc(n_plantyrs_sc))

# Write to file
# write.csv(openest2_spp,
#           "data/estopenflower_restrictions_samplesizes.csv",
#           row.names = FALSE)

# Priority 1 species with >= 25 plant-years:
  # Wild bergamot (38; LA, NM, OK)
  # Eastern purple coneflower (31; LA, NM, OK)
  # Common buttonbush (27, LA, OK)
# Priority 2 species with >= 25 plant-years:
  # Red maple (714; LA)
  # Black elderberry (201; LA)
  # Butterfly milkweed (62; LA, NM, OK, TX)
# Priority 3 species with >= 25 plant-years:
  # Swamp milkweed (26; OK)
# Priority 4 species with >= 25 plant-years:
  # Common milkweed (158; OK)
  # Rubber rabbitbrush (40; NM)
  # New England aster (26; LA, OK)

# What do the red maple data look like?
plantyr2_rm <- plantyr2 %>%
  filter(common_name == "red maple")
plantyr2_rm %>%
  group_by(sc) %>%
  summarize(n_plants = length(unique(individual_id)),
            n_plantyrs = n(),
            mn_n_obs = round(mean(n_obs), 1),
            mn_obs_interval = round(mean(open_est_interval_mn), 1),
            n_plantyrs_interval7 = sum(open_est_interval_mn <= 7)) %>%
  data.frame()

# Attach an indicator to peak_df, identifying plant years with sufficient data
plantyr2$peak_sufficient <- 1
peak_df <- peak_df %>%
  left_join(select(plantyr2, individual_id, year, peak_sufficient), 
            by = c("individual_id", "year"))

# Visualize estimates for wild bergamot ---------------------------------------#

wb <- peak_df %>%
  filter(common_name == "wild bergamot" & peak_sufficient == 1) %>%
  select(-c(latitude, longitude, species_id, genus, species,
            n_observations, max_open, peak_sufficient)) %>%
  arrange(individual_id, year)
count(wb, site_id, individual_id, state) 
count(wb, individual_id, year)
  # 18 individual plants, none of which are in target states.
  # 15-79 observations of a plant in a year

wb_summary <- plantyr2 %>%
  filter(common_name == "wild bergamot")

for (i in 1:nrow(wb_summary)) {
  plot_data <- wb %>% 
    filter(individual_id == wb_summary$individual_id[i], 
           year == wb_summary$year[i])
  
  g <- ggplot() + 
    geom_line(data = plot_data, aes(x = day_of_year, y = num_open)) +
    geom_point(data = filter(plot_data, peak == 1),
               aes(x = day_of_year, y = num_open), 
               shape = 21, fill = "black", size = 2) +
    geom_point(data = filter(plot_data, peak == 0),
               aes(x = day_of_year, y = num_open), 
               shape = 21, fill = "white", size = 2) +
    geom_hline(yintercept = 0.67 * max(plot_data$num_open, na.rm = TRUE),
               linetype = 2) +
    xlim(1, 365) +
    labs(x = "Day of the year", y = "Estimated number of open flowers", 
         title = paste0("Wild bergamot individual ", 
                        plot_data$individual_id[1], " (", 
                        plot_data$state[1], "), ",
                        plot_data$year[1], " [i = ", i, "]")) +
    theme_classic()
  print(g)
}
# Quality of data varies a LOT

# Visualize estimates for purple coneflower -----------------------------------#

epc <- peak_df %>%
  filter(common_name == "eastern purple coneflower" & peak_sufficient == 1) %>%
  select(-c(latitude, longitude, species_id, genus, species,
            n_observations, max_open, peak_sufficient)) %>%
  arrange(individual_id, year)
count(epc, site_id, individual_id, state) 
count(epc, individual_id, year)
# 18 individual plants, none of which are in target states.
# 16-259 observations of a plant in a year

epc_summary <- plantyr2 %>%
  filter(common_name == "eastern purple coneflower")

for (i in 1:nrow(epc_summary)) {
  plot_data <- epc %>% 
    filter(individual_id == epc_summary$individual_id[i], 
           year == epc_summary$year[i])
  
  g <- ggplot() + 
    geom_line(data = plot_data, aes(x = day_of_year, y = num_open)) +
    geom_point(data = filter(plot_data, peak == 1),
               aes(x = day_of_year, y = num_open), 
               shape = 21, fill = "black", size = 2) +
    geom_point(data = filter(plot_data, peak == 0),
               aes(x = day_of_year, y = num_open), 
               shape = 21, fill = "white", size = 2) +
    geom_hline(yintercept = 0.67 * max(plot_data$num_open, na.rm = TRUE),
               linetype = 2) +
    xlim(1, 365) +
    labs(x = "Day of the year", y = "Estimated number of open flowers", 
         title = paste0("Eastern purple coneflower individual ", 
                        plot_data$individual_id[1], " (", 
                        plot_data$state[1], "), ",
                        plot_data$year[1], " [i = ", i, "]")) +
    theme_classic()
  print(g)
}
