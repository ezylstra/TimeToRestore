################################################################################
# Characterize flowering phenophases (onset, duration, end)

# Much of code is derived from that originally developed by Alyssa Rosemartin, 
# Hayley Limes, and Jeff Oliver. 
# See: https://github.com/alyssarosemartin/time-to-restore

# Erin Zylstra
# 2024-05-06
################################################################################

require(dplyr)
require(lubridate)
require(tidyr)

rm(list = ls())

# Load processed NPN status/intensity data
df <- read.csv("data/flower-status-intensities-priorityspp.csv")

# Characterize peak flowering -------------------------------------------------#

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

# Select what criteria we'll use to delineate "peak" flowering:
# 1) estimated number of open flowers within x% of the maximum value for that
#    plant and year,
# 2) criteria above OR estimated number of open flowers > 500.
# peak_def <- "near max"
peak_def <- "near max or over 500"

# Identify the percentage threshold of open flowers, above which, will be
# considered peak flowering.
threshold <- 0.75

# Classify observations of plant on particular date as in peak or not
peak_df <- peak_df %>%
  group_by(individual_id, year) %>%
  mutate(near_max = if_else(num_open >= max_open * threshold, 1, 0)) %>%
  mutate(over_500 = ifelse(num_open > 500, 1, 0)) %>%
  ungroup() %>%
  data.frame()

if (peak_def == "near max") {
  peak_df <- peak_df %>%
    mutate(peak = ifelse(near_max == 1, 1, 0))
} else {
  peak_df <- peak_df %>%
    mutate(peak = ifelse(near_max == 1 | over_500 == 1, 1, 0))
}
# checks:
# head(filter(peak_df, near_max == 0 & over_500 == 1))
# head(filter(peak_df, near_max == 1 & over_500 == 0))
# head(filter(peak_df, near_max == 1 & over_500 == 1))
  
# Will create a new dataframe (peaks) with estimates of onset/end/duration for 
# each plant, year
peaks <- peak_df %>%
  filter(peak == 1) %>%
  group_by(site_id, common_name, individual_id, year) %>%
  summarize(peak_onset = min(day_of_year, na.rm = TRUE),
            peak_end = max(day_of_year, na.rm = TRUE),
            .groups = "keep") %>%
  mutate(peak_duration = peak_end - peak_onset + 1) %>%
  data.frame()

# Want to identify whether peak flowering phenophase was discontinuous or not
# Note that if there was no estimated number of open flowers (num_open = NA)
# between two observations that were considered part of the peak we will
# consider the peak continuous (ie,need a non-NA value of num_open between peak
# observations to consider the peak discontinuous). 
peaks$peak_discontinuous <- NA
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



# New column holds phenophase - Open Flowers
int_pmetric$phenophase_description <- "Open Flowers"
int_pmetric$phenophase_id <- 501

