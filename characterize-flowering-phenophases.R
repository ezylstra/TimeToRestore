################################################################################
# Characterize flowering phenophases (onset, duration, end)

# Much of code is derived from that originally developed by Alyssa Rosemartin, 
# Hayley Limes, and Jeff Oliver. 
# See: https://github.com/alyssarosemartin/time-to-restore

# Erin Zylstra
# 2024-05-21
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
  group_by(common_name, individual_id, partner_group, site_id, site_name, 
           state, year) %>%
  summarize(n_obs = n(),                    # No. of daily observations
            n_obs_fls = sum(!is.na(fl_s)),  # No. of flower status obs
            n_obs_fos = sum(!is.na(fo_s)),  # No. of open flower status obs
            n_obs_fli = sum(!is.na(fl_i)),  # No. of flower intensity values
            n_obs_foi = sum(!is.na(fo_i)),  # No. of open flower intensity values
            n_estimates_open = sum(!is.na(num_open)), # No. of observations with
                                                      # an estiamted number of 
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
spp_ss <- count(filter(samples, n_estimates_open > 1), 
                common_name, priority, LA, NM, OK, TX) %>%
  arrange(desc(n))
spp_ss
summary(spp_ss)
  # Range = 2-1535 (2-295 if you exclude red maple); median = 20. 
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
count(filter(samples, n_estimates_open > 1 & state %in% states4), 
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
  #   select(common_name, observation_date, year, fl_i, fo_i, num_open, max_open)
  # filter(maxfo_annualvar, common_name == "eastern purple coneflower")
  # peak_df %>%
  #   filter(individual_id == 66579) %>%
  #   select(common_name, observation_date, year, fl_i, fo_i, num_open, max_open)

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
# 50%: 9277
# 67%: 8006
# 75%: 7520

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
    if (any(peak_df_sub$fl_s == 0, na.rm = TRUE) | 
        any(peak_df_sub$fo_s == 0, na.rm = TRUE) |
        any(peak_df_sub$peak == 0, na.rm = TRUE)) {
      peaks$peak_discontinuous[i] <- 1  
    } else {
      peaks$peak_discontinuous[i] <- 0
    }
  }
}
count(peaks, peak_discontinuous)
  # About 12% of peaks discontinous (or 25% of peaks with duration > 1 day)
count(peaks, peak_duration, peak_discontinuous)
  # As durations increase, so does probability of discontinous peak


# 122894, 2019 = 1 obs
test1_ind <- which(peaks$individual_id == 122894 & peaks$year == 2019)
i <- test1_ind
# 13583, 2019 = mult obs
test2_ind <- which(peaks$individual_id == 13583 & peaks$year == 2019)
i <- test2_ind


for (i in 1:nrow(peaks)) {
  indiv <- peaks$individual_id[i]
  yr <- peaks$year[i]
  onset <- peaks$peak_onset[i]
  end <- peaks$peak_end[i]
  peak_df_sub <- peak_df %>%
    filter(individual_id == indiv & year == yr) %>%
    filter(day_of_year > onset & day_of_year < end) %>%
    filter(!is.na(num_open))
  peaks$peak_discontinuous[i] <- ifelse(any(peak_df_sub$peak == 0), 1, 0)
}
# summaries of estimated peak duration:
sum(peaks$peak_duration == 1) / nrow(peaks) # 50% of "peaks" have duration == 1
tail(sort(unique(peaks$peak_duration)), 30) # 7 plant-years with duration > 200

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
peak_df %>% filter(individual_id == 13583 & year == 2019) %>% 
  select(state, observation_date, fl_s, fo_s, num_open, peak)


ggplot(data = bumi_73461) + 
  geom_point(aes(x = observation_date, y = num_open)) +
  geom_smooth(aes(x = observation_date, y = num_open))

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

