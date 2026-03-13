################################################################################
# Create template for state info sheet components

# Erin Zylstra
# 2024-07-12
################################################################################

require(dplyr)
require(lubridate)
require(stringr)
require(tidyr)
require(ggplot2)
require(flextable)
require(cowplot)
require(terra)
require(tidyterra)
require(mgcv)
set_null_device(cairo_pdf) # To avoid some warnings about fonts

rm(list = ls())

# Set figure parameters -------------------------------------------------------#

  alphaline <- 0.6 # Transparency setting for lines, points
  alphapoly <- 0.4 # Transparency setting for polygons
  pdfw <- 8.5 # Width of pdf output
  pdfh <- 11  # Height of pdf output
  
  # Color palettes
  whitebox_palette <- "soft"
  # See: https://dieghernan.github.io/tidyterra/reference/scale_whitebox.html
  cols2 <- whitebox.colors(n = 2, palette = whitebox_palette)

# Load csv with sample sizes for weekly proportion of plants with open flowers #

  ofss <- read.csv("data/openflower_weeklyobs_samplesizes.csv")

# Identify priority species with sufficient data ------------------------------#
  
  # Selecting priority species for that average 6 or more observations in the
  # southcentral region per week (for all years combined)
  species <- ofss %>%
    filter(mn_wkobs_sc >= 6) %>%
    select(!contains("z789"))
  
  # List of common names
  common_names <- species$common_name

# Load and subset iNat data ---------------------------------------------------#
  
  inat_zip <- "data/iNat/inaturalist-41spp.zip"
  inat_csv <- str_replace(inat_zip, ".zip", ".csv")
  unzip(zipfile = inat_zip)
  inat_orig <- read.csv(inat_csv)
  file.remove(inat_csv)

  # Extract data for state priority species
  inat <- inat_orig %>%
    select(common_name, species, state, lat, lon, 
           coord_uncert, year, month, day) %>%
    filter(common_name %in% common_names)
  
  # Add state codes and create indicator for SC region
  inat <- inat %>%
    rename(state_name = state) %>%
    mutate(state = state.abb[match(state_name, state.name)],
           sc = 1 * state %in% c("LA", "NM", "OK", "TX"))
  
  # Reformat date
  inat <- inat %>%
    mutate(obsdate = make_date(year, month, day),
           doy = yday(obsdate))
  
# Load processed NPN status/intensity data and format -------------------------#

  df <- read.csv("data/flower-status-intensities-priorityspp.csv")
  
  # Extract data for species that had sufficient data in SC region
  df <- df %>%
    filter(common_name %in% common_names)
  
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
  
# For now, select a state to focus on -----------------------------------------#
  
  st <- "LA"
  
  # Set threshold for species inclusion based on the mean number of plants 
  # observed each week in the SC region (weeks 10-40, across years)
  weekly_nobs_min <- 8
  
  # Extract info just for those species that are priority in the state and have 
  # an average number of weekly open flower observations (in the SC region, 
  # across years) >= threshold listed above.
  species_st <- species %>%
    filter(get(st) == 1) %>%
    filter(mn_wkobs_sc >= weekly_nobs_min)
  common_names_st <- species_st$common_name
  inat_st <- filter(inat, common_name %in% common_names_st)
  df1_st <- filter(df1, common_name %in% common_names_st)
  
# For each species, create figures comparing the number of iNat observations --#
# inside/outside the SC region ------------------------------------------------#

  # Note: won't create iNat figures for species that are likely to be
  # photographed when NOT flowering (ie, red maple)
  
  # spp <- "common_buttonbush"
  spp <- "eastern baccharis"
  
  inat_spp <- inat_st %>%
    filter(common_name == spp)
  
  # If we look at daily totals, lines are really jagged, so exploring weekly 
  # sums instead
  inat_spp <- inat_spp %>%
    mutate(wk = week(obsdate)) %>%
    # Remove observations in week 53
    filter(wk < 53) %>%
    # Create wk_doy columns
    mutate(wk_date1 = parse_date_time(paste(2024, wk, 1, sep = "/"), "Y/W/w"),
           wk_date1 =  as.Date(wk_date1),
           wk_doy1 = yday(wk_date1),
           wk_date4 = parse_date_time(paste(2024, wk, 4, sep = "/"), "Y/W/w"),
           wk_date4 =  as.Date(wk_date4),
           wk_doy4 = yday(wk_date4))
  
  # Get weekly number of observations for all plants, regardless of location
  inat_spp_wk <- inat_spp %>%
    group_by(wk, wk_date4, wk_doy4) %>%
    summarize(n_obs = n(), .groups = "keep") %>%
    mutate(mn_obs = n_obs / length(unique(inat_spp$year))) %>%
    data.frame()
  ggplot(inat_spp_wk, aes(x = wk_date4, y = mn_obs)) +
    geom_line()
    
  # Get weekly number of observations, by region (SC, other)
  inat_spp_wkr <- inat_spp %>%
    group_by(sc, wk, wk_date4, wk_doy4) %>%
    summarize(n_obs = n(), .groups = "keep") %>%
    mutate(mn_obs = n_obs / length(unique(inat_spp$year)),
           sc_f = as.factor(ifelse(sc == 0, "Other", "SC"))) %>%
    data.frame()
  
  # To create figures with ticks between month labels on the x axis, need to do 
  # a little extra work. 
    # Dates where we want month labels (15th of month)
    x_lab <- as.Date(paste0("2024-", 1:12, "-15"))
    # Dates where we want ticks (1st of month)
    x_tick <- as.Date(c(paste0("2024-", 1:12, "-01"), "2025-01-01"))
    n_x_tick <- length(x_tick)
    # Will specify axis breaks & ticks at 1st and 15th of month. Make labels on
    # the 1st black and change color of tick marks on the 15th to NA.
    
  # Figure with a line for each region
  inat_r_plot <- ggplot(inat_spp_wkr, aes(x = wk_date4, y = mn_obs)) +
    geom_line(aes(group = sc_f, color = sc_f)) +
    scale_color_whitebox_d(palette = whitebox_palette, name = "Region") +
    scale_x_continuous(breaks = c(x_lab, x_tick),
                       labels = c(month.abb, rep("", n_x_tick))) +
    labs(y = paste0("Mean number of iNat observations per week")) +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8),
          axis.ticks.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                rep("black", n_x_tick))),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                      rep("white", n_x_tick))),
          panel.background = element_rect(fill = "gray95"),
          axis.title.x = element_blank())
  inat_r_plot
  
  # Just plotting mean number of observations per week in SC region
  inat_sc_plot <- ggplot(filter(inat_spp_wkr, sc == 1)) +
    geom_line(aes(x = wk_date4, y = mn_obs)) +
    scale_x_continuous(breaks = c(x_lab, x_tick),
                       labels = c(month.abb, rep("", n_x_tick))) +
    labs(y = paste0("Mean number of iNat observations per week")) +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8),
          axis.ticks.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                rep("black", n_x_tick))),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                      rep("white", n_x_tick))),
          panel.background = element_rect(fill = "gray95"),
          axis.title.x = element_blank())
  inat_sc_plot

# For each species, create figures comparing weekly open flower proportions ---# 
# nside/outside the SC region -------------------------------------------------#  
  
  # !! Can use a lot of code from weekly-open-flower-prop-species.R script here
  
# Create alternative visualizations of weekly open flower proportions for
# plants in the SC region
  
  # Raw portions with size of point relative to sample sizes
  # Heat map with raw proportions along with sample size bar chart below
  # GAM predictions
  # Overlaying raw proportions and GAM predictions
  
  # !! Can use a lot of code from weekly-open-flower-prop-species.R script here
  