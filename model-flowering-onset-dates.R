################################################################################
# Model variation in flowering onset dates

# Erin Zylstra
# 2024-07-31
################################################################################

require(dplyr)
require(lubridate)
require(stringr)
require(tidyr)
require(ggplot2)
require(ggpubr)
require(lme4)

rm(list = ls())

# Load relevant data ----------------------------------------------------------#

  # Flowering onset data, by plant-wateryr
  onset <- read.csv("data/flowering-onsets-2013-2023.csv")
  
  # Species-level summaries of flowering onset data (for species-phenophase
  # combinations that had a minimum of 30 plant-yrs of data)
  spp_onsets <- read.csv("output/flowering-onsets/flower-onset-data-summaries.csv")
 
# Select parameters for data inclusion/exclusion ------------------------------#
  
  # Set the maximum number of days prior to a "yes" that we'll need a "no" in 
  # order to use the observation in onset models
  daysprior <- 14    

# Reformat and subset onset data ----------------------------------------------#
  
  # Extract plant-wateryr data for flower onset
  onset_fl <- onset %>%
    select(-c(contains("fo"), first_obs, last_obs)) %>%
    rename_with(str_replace, pattern = "fl_", replacement = "") %>%
    rename_with(str_replace, pattern = "_fl", replacement = "") %>%
    mutate(phenophase = "flower") %>%
    relocate(phenophase, .after = "n_obs")
  
  # Extract plant-wateryr data for open flower onset
  onset_fo <- onset %>%
    select(-c(contains("fl"), first_obs, last_obs)) %>%
    rename_with(str_replace, pattern = "fo_", replacement = "") %>%
    rename_with(str_replace, pattern = "_fo", replacement = "") %>%
    mutate(phenophase = "open flower") %>%
    relocate(phenophase, .after = "n_obs") 
  
  # Combine the two dataframes. Now each row contains onset data for a unique 
  # combination of plant, wateryr, and phenophase
  onset <- rbind(onset_fl, onset_fo)
  
  # Remove plant-yrs when:
  # 1) flowers/open flowers weren't observed (firstyes = NA),
  # 2) "yes" to phenophase status wasn't preceded by a "no" (dayssince_lastno = NA)
  # 3) preceding "no" occurred more than XX days ago (dayssince_lastno > daysprior)
  onset <- onset %>%
    filter(!is.na(firstyes)) %>%
    filter(!is.na(dayssince_lastno)) %>%
    filter(dayssince_lastno <= daysprior)
  # count(onset, phenophase)

# Prioritize species to use in models -----------------------------------------#

  # Priority 1: 2 or more locations in SC region
  # Priority 2: 1 location in SC region
  # Priority 3: nearest plant location < 200 km from SC region
  # If a species-phenophase doesn't meet any of these criteria, probably won't 
    # create a model for onset date (and leaving priority = NA)

  spp_onsets <- spp_onsets %>%
    rowwise() %>%
    mutate(SC = sum(c_across(LA:TX))) %>%
    mutate(priority = case_when(
      SC > 1 ~ 1,
      SC == 1 ~ 2,
      dist_to_SC_km < 200 ~ 3,
      .default = NA
    )) %>%
    data.frame()

  # Remove species from onset data if they're not priority 1-3 (Will keep onset
  # data for both phenophases if the species is considered priority 1-3 for
  # at least one phenophase)
  spp_all <- spp_onsets %>%
    filter(!is.na(priority)) %>%
    select(common_name) %>%
    distinct() %>%
    unlist(use.names = FALSE)
  
  onset <- onset %>%
    filter(common_name %in% spp_all)
  
# Load weather data (download data if not done yet) ---------------------------#  
  
  sites <- onset %>%
    select(site_id, lon, lat) %>%
    distinct() 
  
  nsites <- nrow(sites)
  yrs <- (min(onset$wateryr) - 1):max(onset$wateryr)
  
  weather_filename <- paste0("data/climate/seasonal-weather-",
                             nsites, "sites-", 
                             first(yrs), "-", last(yrs), ".csv")

  # If not done already, download daymet data and summarize
  # Note: this will take quite a bit of time!
  if (!file.exists(weather_filename)) {
    source("download-daymet.R")
  }
  
  # Load data
  weather <- read.csv(weather_filename)
  
  # Seasons (as previously defined by NPN)
    # Spring = Mar-May
    # Summer = Jun-Aug
    # Fall = Sep-Nov
    # Winter = Dec-Feb (assigned to year for Jan-Feb)
  # Weather data
    # prcp = Accumulated precipitation (mm)
    # tmin = Mean of daily minimum temperatures (degC)
    # tmax = Mean of daily maximum temperatures (degC)
  
# Test workflow/run for one species -------------------------------------------#
  
  # Important to remember that all day numbers in the onset dataset are actually
  # day-of-WATER-year (dowy). Water year runs 1 Oct - 30 Sep (eg, wateryr 2022 =
  # 1 Oct 2021 - 30 Sep 2022). Creating simple table to look at important day 
  # numbers (no leap year).
  dowy <- data.frame(first_date = c(as.Date(paste0("2021-", 10:12, "-01")),
                                    as.Date(paste0("2022-", 1:12, "-01"))))
  dowy <- dowy %>%
    mutate(mon = month(first_date),
           n_days = days_in_month(first_date),
           dowy = as.numeric(first_date - make_date(2021, 9, 30)),
           doy = dowy - 92) %>%
    mutate(dowy = ifelse(dowy > 365, NA, dowy),
           doy = ifelse(doy < 1, NA, doy))
  dowy
    # Jan 1 = 93
    # Mar 1 = 152
    # Jun 1 = 244
    # Sep 1 = 336
  
  phenophases <- c("flower", "open flower")
  priority_levels <- 1:3
  
  # For now, picking common buttonbush, open flower (but structure is here
  # to do things in a loop)
  i = 2
  j = 1 
  pheno <- phenophases[i]
  spp_list <- spp_onsets %>%
    filter(phenophase == pheno) %>%
    filter(priority == priority_levels[j]) %>%
    select(common_name) %>%
    unlist(use.names = FALSE)
  
  spp <- spp_list[2]
  
  onsetsub <- onset %>%
    filter(phenophase == pheno) %>%
    filter(common_name == spp)

  # Identify when the plant flowers so we know what climate data to grab (to
  # ensure it's data preceding the event and not after). 
  median_onset <- median(onsetsub$firstyes)
  # If onset is before 1 Apr (dowy < 183), then use spring and summer weather 
    # from previous calendar year.
  # If onset is 1 Apr - 30 June (dowy = 183 - 273), then use spring weather in 
    # current year and summer weather in previous year
  # If onset is on or after 1 Jul (dowy = 274), then use spring and summer 
    # weather from current calendar year.
  spring <- ifelse(median_onset < 183, "previous", "current")
  summer <- ifelse(median_onset %in% 183:273, "previous", "current")
  
  # Put weather data in wide form
  weather <- weather %>%
    pivot_wider(names_from = season, 
                values_from = c(prcp, tmin, tmax),
                names_glue = "{season}_{.value}")
  
  # Extract weather data for sites, years in onset dataset
  weathersub <- expand.grid(
    site_id = sort(unique(onsetsub$site_id)),
    onset_yr = sort(unique(onsetsub$wateryr)),
    KEEP.OUT.ATTRS = FALSE
  )
  weathersub <- weathersub %>%
    rowwise() %>%
    mutate(spring_yr = ifelse(spring == "current", onset_yr, onset_yr - 1),
           summer_yr = ifelse(summer == "current", onset_yr, onset_yr - 1)) %>%
    mutate(fall_yr = onset_yr - 1,
           winter_yr = onset_yr)
  weathersub <- weathersub %>%
    left_join(select(weather, site_id, seasonyr, contains("spring")),
              by = c("site_id", "spring_yr" = "seasonyr")) %>%
    left_join(select(weather, site_id, seasonyr, contains("summer")),
              by = c("site_id", "summer_yr" = "seasonyr")) %>%
    left_join(select(weather, site_id, seasonyr, contains("fall")),
              by = c("site_id", "fall_yr" = "seasonyr")) %>%
    left_join(select(weather, site_id, seasonyr, contains("winter")),
              by = c("site_id", "winter_yr" = "seasonyr")) %>%
    select(-c(spring_yr, summer_yr, fall_yr, winter_yr))
  
  # Add weather data to onset data
  onsetsub <- onsetsub %>%
    left_join(weathersub, by = c("site_id", "wateryr" = "onset_yr"))

  # View scatterplots: matrix of onset date vs weather variable (rows = 
  # season, columns = precip, tmin, tmax)
  onsetsubl <- onsetsub %>%
    select(common_name, plantwateryr, firstyes, contains("spring"), 
           contains("summer"), contains("fall"), contains("winter")) %>%
    pivot_longer(cols = spring_prcp:winter_tmax, 
                 names_to = "weather_var",
                 values_to = "value") %>%
    separate(weather_var, sep = "_", into = c("season", "var")) %>%
    mutate(firstyes_doy = firstyes - 92) %>%
    mutate(firstyes_date = as.Date(firstyes_doy, origin = "2021-12-31")) %>%
    mutate(Var = case_when(
      var == "prcp" ~ "Precipitation (mm)",
      var == "tmax" ~ "Maximum temperature (C)",
      var == "tmin" ~ "Minimum temperature (C)")) %>%
    mutate(Season = factor(str_to_sentence(season),
                           levels = c("Winter", "Spring", "Summer", "Fall")))
  
  date_breaks <- dowy %>%
    filter(!is.na(doy)) %>%
    filter(mon %in% seq(1, 12, by = 3))
  
  ggplot(data = onsetsubl, aes(x = value, y = firstyes_doy)) +
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x) +
    facet_grid(rows = vars(Season), cols = vars(Var), scales = "free_x") +
    stat_cor(method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01) +
    scale_y_continuous(limits = c(1, 366), 
                       breaks = c(date_breaks$doy, 366),
                       labels = format(c(date_breaks$first_date, as.Date("2023-01-01")), 
                                       "%b")) +
    labs(x = "", y = "Onset date", 
         title = paste0(str_to_sentence(spp), ", ", pheno, " onset date"))
  
  # Get table with pearson corr coefs and associated p-values
  cors <- onsetsubl %>%
    group_by(season, var) %>%
    summarize(cor = cor.test(firstyes_doy, value)$estimate,
              p.cor = cor.test(firstyes_doy, value)$p.value,
              , .groups = "keep") %>%
    data.frame() %>%
    filter(p.cor < 0.10) %>%
    mutate(seas_var = paste0(season, "_", var))

  # Note that the r, P values are for simple linear correlations (Pearson), and
  # don't take things like repeated measures into account
  
  # Maybe we can use these figures to identify weather variables with the most
  # explanatory power (to include in linear models). Start with P <= 0.10?
  
  # After looking at results from just a few linear models, I'm not sure these
  # plots are helpful at all
  
  # Models that Alyssa et al. ran for this species/phenophase previously had
  # very different results. Looks like that's primarily because data from 2023
  # are very different. Majority of onset dates in 2023 are after doy 200, which 
  # rarely happened before 2023. Relationships between onset dates and weather
  # at least fall and spring tmins) seemed more significant before 2023 data 
  # were added in....
  
  # Correlation (firstyes ~ spring min temp) with all years data
  ggplot(onsetsub, aes(x = spring_tmin, y = firstyes_doy)) + 
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x) +
    stat_cor(method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01)
  
  # Correlation excluding 2023
  ggplot(filter(onsetsub, wateryr < 2023), 
         aes(x = spring_tmin, y = firstyes_doy)) + 
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x) +
    stat_cor(method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01)
  
  # Onset dates, all years
  select(onsetsub, site_id, plant_id, wateryr, firstyes_doy, fall_tmin)
  # Onset dates, 2023
  select(filter(onsetsub, wateryr == 2023), 
         site_id, lat, plant_id, wateryr, firstyes_doy, fall_tmin)


  # Run a few linear models with random effects
  seasons <- c("winter", "spring", "summer", "fall")
  vars <- c("prcp", "tmax", "tmin")
  
    seas <- seasons[2]
    var <- vars[3]
    seas_var <- paste0(seas, "_", var)
    m1formula <- as.formula(paste0("firstyes ~ ", seas_var, " + (1|site_id)"))
    m1 <- lmer(m1formula, data = onsetsub)
    summary(m1)
    m2formula <- as.formula(paste0("firstyes ~ ", seas_var, " + lat + (1|site_id)"))
    m2 <- lmer(m2formula, data = onsetsub)
    summary(m2)
    m3formula <- as.formula(paste0("firstyes ~ ", seas_var, " * lat + (1|site_id)"))
    m3 <- lmer(m3formula, data = onsetsub)
    summary(m3)
 
  anova(m2, m1)
  anova(m3, m1)

  summary(msimp1 <- lm(firstyes ~ spring_tmax, data = onsetsub))
  summary(msimp2 <- lm(firstyes ~ spring_tmax + lat, data = onsetsub))
  summary(msimp3 <- lm(firstyes ~ spring_tmax * lat, data = onsetsub))

  # Do we want to standardize variables before running models?