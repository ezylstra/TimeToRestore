################################################################################
# Characterize flowering phenophases (onset, duration, end)

# Much of code is derived from that originally developed by Alyssa Rosemartin, 
# Hayley Limes, and Jeff Oliver. 
# See: https://github.com/alyssarosemartin/time-to-restore

# Erin Zylstra
# 2024-05-24
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
            n_fls = sum(!is.na(status_fl)),   # No. of flower status obs
            n_fos = sum(!is.na(status_fo)),   # No. of open flower status obs
            n_fli = sum(!is.na(midpoint_fl)), # No. of flower intensity values
            n_foi = sum(!is.na(midpoint_fo)), # No. of open flower intensity values
            n_ests_open = sum(!is.na(num_open)), # No. of observations with
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
summary(spp_ss)
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

# Characterize "peak" flowering -----------------------------------------------#
# Based on estimated number of open flowers

# Remove plant-year combinations where we have no estimates of the number of 
# open flowers
peak_df <- df %>%  
  group_by(individual_id, year) %>%
  filter(any(num_open != 0)) %>%
  ungroup() %>%
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
# peak_def <- "near max"
peak_def <- "near max or over 500"

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
# 50%: 9326
# 67%: 8039
# 75%: 7551

# Will use 67% threshold for now.
if (peak_def == "near max") {
  peak_df <- peak_df %>%
    mutate(peak = ifelse(near_max67 == 1, 1, 0))
} else {
  peak_df <- peak_df %>%
    mutate(peak = ifelse(near_max67 == 1 | over_500 == 1, 1, 0))
}
# checks:
# head(filter(peak_df, near_max67 == 0 & over_500 == 1))
# head(filter(peak_df, near_max67 == 1 & over_500 == 0))
# head(filter(peak_df, near_max67 == 1 & over_500 == 1))
  
# Will create a new dataframe (peaks) with summaries of data for each plant &
# year, including the number of daily observations and estimates of 
# peak onset/end/duration for each plant, year
peaks <- peak_df %>%
  filter(peak == 1) %>%
  group_by(site_id, common_name, individual_id, year) %>%
  summarize(peak_onset = min(day_of_year, na.rm = TRUE),
            peak_end = max(day_of_year, na.rm = TRUE),
            .groups = "keep") %>%
  mutate(peak_duration = peak_end - peak_onset + 1) %>%
  data.frame()

# Identify whether peak flowering phenophase was discontinuous or not.
# Note that we need a flower or open_flower status of 0, or an estimated number
# of flowers < 500 and < threshold*max between observations classifed as peak 
# to label the peak as discontinuous.
peaks$peak_discontinuous <- NA
for (i in 1:nrow(peaks)) {
  indiv <- peaks$individual_id[i]
  yr <- peaks$year[i]
  onset <- peaks$peak_onset[i]
  end <- peaks$peak_end[i]
  peak_df_sub <- peak_df %>%
    filter(individual_id == indiv & year == yr) %>%
    filter(day_of_year > onset & day_of_year < end)
  if (nrow(peak_df_sub) == 0) { 
    peaks$peak_discontinuous[i] <- 0
  } else {
    if (any(peak_df_sub$status_fl == 0, na.rm = TRUE) | 
        any(peak_df_sub$status_fo == 0, na.rm = TRUE) |
        any(peak_df_sub$peak == 0, na.rm = TRUE)) {
      peaks$peak_discontinuous[i] <- 1  
    } else {
      peaks$peak_discontinuous[i] <- 0
    }
  }
}

# Summaries of estimated peak duration/continuity:
sum(peaks$peak_duration == 1) / nrow(peaks) 
  # 51% of "peaks" have duration == 1
tail(sort(unique(peaks$peak_duration)), 30) 
  # 7 plant-years with duration > 200
sum(peaks$peak_discontinuous == 1) / nrow(peaks)
sum(peaks$peak_discontinuous[peaks$peak_duration > 1] == 1) / sum(peaks$peak_duration > 1)
  # About 12% of peaks discontinous (or 25% of peaks with duration > 1 day)
count(peaks, peak_duration, peak_discontinuous)
summary(glm(peak_discontinuous ~ peak_duration, 
            data = peaks, family = "binomial"))
  # As durations increase, so does probability of discontinous peak

# Attach relevant sample sizes
peaks <- peaks %>%
  left_join(select(samples, -c(common_name, site_id)), 
            by = c("individual_id", "year")) %>%
  relocate(starts_with("peak_"), .after = "TX")

# Visualize "peak" flowering --------------------------------------------------#
# Need to figure out how we're summarizing over years, individuals/locations

# Could plot number of open flowers by day (using different symbols to identify
# peak or not) for each plant, year. [Need to use peak_df for this]

peak_df %>% filter(common_name == "butterfly milkweed") %>% group_by(year) %>%
  summarize(n_plants = length(unique(individual_id)),
            n_obs = n(),
            mn_obs_per_yr = n_obs / n_plants)
peak_df %>% filter(common_name == "butterfly milkweed" & year == 2019) %>% 
  group_by(individual_id) %>%
  summarize(n_obs = n()) %>%
  data.frame()

# Here's an example of a plant that was observed 20 times in 2019 by the same
# person, consistently Jan-Apr and Sep-Dec. Lots of observations without 
# intensity values, which makes things difficult....
  bumi_14303_2019 <- peak_df %>% 
    filter(individual_id == 14303 & year == 2019) %>% 
    select(state, observation_date, day_of_year, status_fl, status_fo, 
           midpoint_fl, midpoint_fo, num_open, peak, person_id)
  # For visualizing, if status_fo == 0, then need num_open == 0
  bumi_14303_2019 <- bumi_14303_2019 %>%
    filter(!is.na(status_fl) & !is.na(status_fo)) %>%
    mutate(num_open = ifelse(is.na(num_open) & (status_fl == 0 | status_fo == 0), 
                             0, num_open)) %>%
    mutate(peak = ifelse(num_open == 0, 0, peak))
  bumi_14303_2019

# How should we accurately convey this information?
  # Is it fair to connect estimates of open flowers with a line if there's a lot
    # of time between estimates of observations?
  # Especially when there are observations without estimates of open flowers, 
    # it would like to also plot status information, but using solid "blocks" of 
    # flower/open flower time assuming that observations were often enough 
    # that we didn't miss status changes.
  
  # Could use stuff like this to create "flower" or "open flower" periods
  data.frame(unclass(rle(bumi_14303_2019$status_fl)))
  geom_rect(aes(xmin = 32, xmax = 263, ymin = 0, ymax = 50), 
            fill = "orange", alpha = 0.2, col = NA)
  
  ggplot() + 
    geom_line(data = bumi_14303_2019, aes(x = day_of_year, y = num_open)) +
    geom_point(data = filter(bumi_14303_2019, peak == 1),
               aes(x = day_of_year, y = num_open), 
               shape = 21, fill = "black", size = 2) +
    geom_point(data = filter(bumi_14303_2019, peak == 0),
               aes(x = day_of_year, y = num_open), 
               shape = 21, fill = "white", size = 2) +
    geom_hline(yintercept = 0.67 * max(bumi_14303_2019$num_open, na.rm = TRUE),
               linetype = 2) +
    labs(x = "Day of the year", y = "Estimated number of open flowers", 
         title = "Butterfly milkweed individual 14303, 2019") +
    theme_classic()
  

# Plant observed 46 times in 2019, all by the same observer
  test_peak <- peak_df %>% 
    filter(individual_id == 141477 & year == 2019) %>% 
    select(state, day_of_year, status_fl, status_fo, midpoint_fl, midpoint_fo, 
           num_open, peak, person_id)
  test_peak
  # For visualizing, if status_fo == 0, then need num_open == 0
  test_peak <- test_peak %>%
    filter(!is.na(status_fl) & !is.na(status_fo)) %>%
    mutate(num_open = ifelse(is.na(num_open) & (status_fl == 0 | status_fo == 0), 
                             0, num_open)) %>%
    filter(!is.na(num_open)) %>%
    mutate(peak = replace_na(peak, 0))
  test_peak
  ggplot() + 
    geom_line(data = test_peak, aes(x = day_of_year, y = num_open)) +
    geom_point(data = filter(test_peak, peak == 1),
               aes(x = day_of_year, y = num_open), 
               shape = 21, fill = "black", size = 2) +
    geom_point(data = filter(test_peak, peak == 0),
               aes(x = day_of_year, y = num_open), 
               shape = 21, fill = "white", size = 2) +
    labs(x = "Day of the year", y = "Estimated number of open flowers", 
         title = "Butterfly milkweed individual 141477, 2019")




# Could create a 7-day heatmap or "activity curve" with the proportion of
# individuals (over some area and years) that are in peak flowering [Need to
# use peak_df for that]

# Onset/end/duration good for linear models, where we can use the response for
# an individual and year and include explanatory variables like climate, year,
# and/or latitude



# Characterize flower phenophase ----------------------------------------------#
# Similar to peak, want to get an onset, duration, end, discontinuous flag for
# each plant/year. Unlike peak, this is likely easier done with the individual
# phenometrics dataset

