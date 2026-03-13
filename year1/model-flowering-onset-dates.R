################################################################################
# Model variation in flowering onset dates

# Erin Zylstra
# 2024-08-05
################################################################################

require(dplyr)
require(lubridate)
require(stringr)
require(tidyr)
require(ggplot2)
require(ggpubr)
require(nlme)
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

  # Priority 1: 2 or more plants located in SC region
  # Priority 2: 1 plant located in SC region
  # Priority 3: nearest plant location < 200 km from SC region
  # If a species-phenophase doesn't meet any of these criteria, probably won't 
    # create a model for onset date (leaving priority = NA)

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
  
# Load and format weather data (download data if not done yet) ----------------#  
  
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
  
  # Put weather data in wide form
  weather <- weather %>%
    pivot_wider(names_from = season, 
                values_from = c(prcp, tmin, tmax),
                names_glue = "{season}_{.value}")
  
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
  
  # Extract onset data for species, phenophase
  onsetsub <- onset %>%
    filter(phenophase == pheno) %>%
    filter(common_name == spp)
  
  # Look at distribution of measurements among sites and years
  count(onsetsub, site_id) # More than half of sites have 1 measurement
    count(count(onsetsub, site_id), n)
  count(onsetsub, wateryr) # 12 measurements in 2023, all other years 1-6
  count(onsetsub, site_id, wateryr) # 1-3 observations per site-year
  onsetsub %>%
    group_by(lat) %>%
    summarize(n_obs = n(),
              n_years = length(unique(wateryr))) %>%
    data.frame()
    # 2 sites ~ 37 deg lat that were observed in 5 years. Most other 
    # sites only have data from 1 year.
  
  # Identify when the plant flowers so we know what weather data to grab (to
  # ensure weather data precedes the event). 
  median_onset <- median(onsetsub$firstyes)
  # If onset is before 1 Apr (dowy < 183), then use spring and summer weather 
    # from previous calendar year.
  # If onset is 1 Apr - 30 June (dowy = 183 - 273), then use spring weather in 
    # current year and summer weather in previous year
  # If onset is on or after 1 Jul (dowy = 274), then use spring and summer 
    # weather from current calendar year.
  spring <- ifelse(median_onset < 183, "previous", "current")
  summer <- ifelse(median_onset < 274, "previous", "current")
  
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
  # Add firstyes doy (not dowy) 
  onsetsub <- onsetsub %>%
    mutate(firstyes_doy = firstyes - 92)
  
  # Put data in long form, for easy plotting
  onsetsubl <- onsetsub %>%
    select(common_name, plantwateryr, wateryr, firstyes, firstyes_doy, 
           contains("spring"), contains("summer"), contains("fall"), 
           contains("winter")) %>%
    pivot_longer(cols = spring_prcp:winter_tmax, 
                 names_to = "weather_var",
                 values_to = "value") %>%
    separate(weather_var, sep = "_", into = c("season", "var")) %>%
    mutate(Var = case_when(
      var == "prcp" ~ "Precipitation (mm)",
      var == "tmax" ~ "Maximum temperature (C)",
      var == "tmin" ~ "Minimum temperature (C)")) %>%
    mutate(Season = factor(str_to_sentence(season),
                           levels = c("Winter", "Spring", "Summer", "Fall")))

  # Create scatterplots: matrix of onset date vs weather variable (rows = 
  # season, columns = precip, tmin, tmax). NOTE: these scatterplots and r, P 
  # listed on each plot do not take things like repeated measures or latitude
  # into account.
  date_breaks <- dowy %>%
    filter(!is.na(doy)) %>%
    filter(mon %in% seq(1, 12, by = 3))
  
  # Scatterplot with all data
  ggplot(data = onsetsubl, aes(x = value, y = firstyes_doy)) +
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x) +
    facet_grid(rows = vars(Season), cols = vars(Var), scales = "free_x") +
    stat_cor(method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01) +
    scale_y_continuous(limits = c(1, 366), 
                       breaks = c(date_breaks$doy, 366),
                       labels = format(c(date_breaks$first_date, 
                                         as.Date("2023-01-01")), 
                                       "%b")) +
    labs(x = "", y = "Onset date", 
         title = paste0(str_to_sentence(spp), ", ", pheno, " onset date"))
  
  # Get table with pearson corr coefs and associated p-values
  # cors <- onsetsubl %>%
  #   group_by(season, var) %>%
  #   summarize(cor = cor.test(firstyes_doy, value)$estimate,
  #             p.cor = cor.test(firstyes_doy, value)$p.value,
  #             , .groups = "keep") %>%
  #   data.frame()
  
  # Same scatterplot, but highlight data from 2023 in red
  onsetsubl <- onsetsubl %>%
    mutate(y2023 = ifelse(wateryr == 2023, "2023", paste0(min(yrs), "-2022")),
           y2023 = factor(y2023, levels = c(paste0(min(yrs), "-2022"), "2023")))
  ggplot(data = onsetsubl, aes(x = value, y = firstyes_doy)) +
    geom_point(aes(group = y2023, color = y2023)) +
    scale_color_manual(values = c("black","salmon"), name = "Year") +
    stat_smooth(method = "lm", formula = y ~ x) +
    facet_grid(rows = vars(Season), cols = vars(Var), scales = "free_x") +
    stat_cor(method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01) +
    scale_y_continuous(limits = c(1, 366), 
                       breaks = c(date_breaks$doy, 366),
                       labels = format(c(date_breaks$first_date, 
                                         as.Date("2023-01-01")), 
                                       "%b")) +
    labs(x = "", y = "Onset date", 
         title = paste0(str_to_sentence(spp), ", ", pheno, " onset date"))  
  
  # Some onset dates in 2023 are much later than other years. When late onset
  # dates were associated with warmer tmin or tmax, it reduced significance of
  # correlations between temperature and onset date based on 2012-2022 data.
  
  # Look at relationship between spring min temps and onset date including/
  # excluding 2023 data
  cols <- c("2012-2023" = "darkorchid2", "2012-2022" = "gray20")
  ggplot(data = onsetsub, aes(x = spring_tmin, y = firstyes_doy)) +
    stat_smooth(data = onsetsub,
                method = "lm", formula = y ~ x, 
                aes(col = "2012-2023", fill = "2012-2023"), alpha = 0.2) +
    stat_smooth(data = filter(onsetsub, wateryr < 2023),
                method = "lm", formula = y ~ x, 
                aes(col = "2012-2022", fill = "2012-2022"), alpha = 0.2) +  
    geom_point(data = onsetsub, 
               aes(color = "2012-2023")) +
    geom_point(data = filter(onsetsub, wateryr < 2023),
               aes(color = "2012-2022")) +  
    scale_color_manual(values = cols, name = "Years") +
    scale_fill_manual(values = cols, guide = "none") +
    stat_cor(data = onsetsub,
             aes(color = "2012-2023"),      
             method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01, 
             label.x = 3, label.y = 75,
             show.legend = FALSE) +
    stat_cor(data = filter(onsetsub, wateryr < 2023),
             aes(color = "2012-2022"),      
             method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01, 
             label.x = 3, label.y = 82,
             show.legend = FALSE) +
    labs(x = "Minimum spring temperature (degC)",
         y = "Open flower onset DOY") +
    theme_bw()

  # Comparison of data: 2023 vs 2012-2022
  onsetsubl %>%
    filter(season == "spring", var == "tmin") %>%
    group_by(y2023) %>%
    summarize(n_onset = n(),
              # spr_tmin_min = round(min(value), 1),
              # spr_tmin_mean = round(mean(value), 1),
              # spr_tmin_max = round(max(value), 1),
              onset_min = min(firstyes_doy),
              onset_median = round(median(firstyes_doy)),
              onset_mean = round(mean(firstyes_doy)),
              onset_max = max(firstyes_doy)) %>%
    data.frame()
  
  # Exploring variation in onset dates to identify appropriate random effect
  # structure in linear models
  
  # Quick look to see how much onset dates vary by site (here, a lot)
  ggplot(data = onsetsub, aes(x = as.factor(site_id), y = firstyes_doy)) +
    geom_boxplot() +
    geom_jitter(color = "blue", width = 0.1, height = 0)
  
  # Quick look for correlations between onset dates and latitude (yes, positive)
  ggplot(data = onsetsub, aes(x = lat, y = firstyes_doy)) +
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x) +
    labs(x = "Latitude", y = "Onset DOY") +
    stat_cor(method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01)
  
  # Is latitude correlated with weather variables? (yes, highly negatively
  # correlated with all temperature variables [not with precip])
  round(cor(select(onsetsub, contains("winter"), contains("spring"),
                   contains("summer"), contains("fall")), onsetsub$lat), 2)

  # Quick look to see how much onset dates vary by year
  ggplot(data = onsetsub, aes(x = as.factor(wateryr), y = firstyes_doy)) +
    geom_boxplot() +
    geom_jitter(color = "blue", width = 0.1, height = 0)
  

  # Run a few linear models with random effects
  # (Note: using firstyes_doy and not firstyes since dates are later in year)
  
  seasons <- c("winter", "spring", "summer", "fall")
  vars <- c("prcp", "tmax", "tmin")
  
    # Pick one weather variable for now (spring_tmin)
    seas <- seasons[2]
    var <- vars[3]
    seas_var <- paste0(seas, "_", var)
    
    # Code I could use to run things in a loop over different weather vars...
    # m1formula <- as.formula(paste0("firstyes ~ ", seas_var, " + (1|site_id)"))
    # m1 <- lmer(m1formula, data = onsetsub)
    
    # Using lme4 package so we can create more complex random effect structures
    # (and hardcoding model formulas for now)
      
    # Model with year and site random effects (crossed)
    m.rsiteyr <- lmer(firstyes_doy ~ spring_tmin + (1|site_id) + (1|wateryr),
                      data = onsetsub, REML = TRUE)
    summary(m.rsiteyr)
    confint(m.rsiteyr)
      # Negative effect of spring min temps (though 95% CI does span 0)
      # Massive variation among sites, much less among years
    
      # On a side note: is a model with latitude instead of spring temps better?
      mlat.rsiteyr <- lmer(firstyes_doy ~ lat + (1|site_id) + (1|wateryr),
                           data = onsetsub, REML = TRUE)
      summary(mlat.rsiteyr)
      confint(mlat.rsiteyr)
      anova(mlat.rsiteyr, m.rsiteyr) # Model values almost identical...
      # Maybe not surprising given that so few sites that were observed
      # in multiple years. Can't separate temperature and latitudinal effects.
      # Could really use more data in the SC region and/or more sites observed 
      # in multiple years.
    
    # Compare model with spring temp + 2 RE with model that doesn't have year REs
    m.rsite <- lmer(firstyes_doy ~ spring_tmin + (1|site_id),
                    data = onsetsub, REML = TRUE)
    summary(m.rsite)
      # Smaller (and less signif) effect of spring min temps
    
    # For some reason, using anova refits models with ML (Though I think REML
    # is appropriate when the fixed effect structure is the same. ML is 
    # appropriate when comparing models with/without REs and models with 
    # different fixed effects strutures. See here for clear explanation:
    # https://ourcodingclub.github.io/tutorials/mixed-models/ which is based on 
    # Zuur et al. 2009. 
    # anova(m.rsiteyr, m.rsite)
    
    # logLik for model with both random effects about 1.5 points lower
    logLik(m.rsiteyr); logLik(m.rsite)
    AIC(m.rsiteyr); AIC(m.rsite)

    # Evaluate model with spring min temp and 2 random effects
    coef(m.rsiteyr) # look at site- and yr-specific intercepts
      # Interesting that 2021 and 2023 didn't have the highest intercepts
      # Is that because higher latitude sites were sampled in those years?
      ggplot(data = onsetsub, aes(x = as.factor(wateryr), y = lat)) +
        geom_boxplot() +
        geom_jitter(color = "blue", width = 0.1, height = 0)
    
    plot(m.rsiteyr) # resid vs fitted values. One much lower fitted
      # value stands out. Without that site, it's possible that there's a funnel
      # shape (less variance in residuals at higher fitted values). That seems
      # to track with plot below, where there's more variance in onset dates
      # at higher spring tmin values (lower fitted values)
      ggplot(onsetsub, aes(x = spring_tmin, y = firstyes_doy)) + 
        geom_point() +
        stat_smooth(method = "lm", formula = y ~ x) +
        stat_cor(method = "pearson", cor.coef.name = "r",
                 r.accuracy = 0.01, p.accuracy = 0.01, 
                 label.x = 15, label.y= 250)
    
        # Could try Levene's test to see if non-contant variance among groups
        # though many don't love it.
        car::leveneTest(onsetsub$firstyes_doy, onsetsub$site_id, center = "mean")
          # P = 0.20, so not convincing evidence that variance is non-constant
        # If this had been a problem, could theoretically use the nlme package
        # and assign weights to groups (varIdent(form = ~1|site_id)), but I 
        # think this is too complex for our data since many groups have just a 
        # single measurement.

    qqnorm(resid(m.rsiteyr))
    qqline(resid(m.rsiteyr))
      # This pattern (points on right above line, points on left below line)
      # indicates that residuals aren't normally distributed. Heavy-tailed 
      # distribution (too many extreme resids)
    
    # Other checks with residuals?
      # Do they vary with latitude?
      plot(resid(m.rsiteyr) ~ onsetsub$lat) 
        # More variation @ higher latitudes but not necessarily a pattern 
        # indicating that we need to add latitude to the model
      
      # Do residuals vary with year?
      ggplot(data = onsetsub, aes(x = as.factor(wateryr), 
                                  y = resid(m.rsiteyr))) +
        geom_boxplot() +
        geom_jitter(color = "blue", width = 0.1, height = 0)
        # Looks decent (no obvious pattern) with model that includes wateryr as 
        # random effect

# See what models look like for data-rich species (red maple) -----------------#  
  
  i = 2
  j = 1 
  pheno <- phenophases[i]
  spp_list <- spp_onsets %>%
    filter(phenophase == pheno) %>%
    filter(priority == priority_levels[j]) %>%
    select(common_name) %>%
    unlist(use.names = FALSE)
  spp_list
  
  spp <- spp_list[1]
  
  # Extract onset data for species, phenophase
  onsetsub <- onset %>%
    filter(phenophase == pheno) %>%
    filter(common_name == spp)
  
  # Using just data in/near SC (only SC state with plants is LA, but going to
  # include plants in MS that are very close to LA)
  onsetsc <- onsetsub %>%
    filter(state %in% c("LA", "MS"))
    count(onsetsc, site_id) # Two sites have lots of data (32, 54 measurements)
    count(onsetsc, wateryr) # 10-30 meas/yr in 2017-2023
    count(onsetsc, site_id, wateryr) # 1-9 observations per site-year

  # Identify when the plant flowers so we know what weather data to grab (to
  # ensure weather data precedes the event). 
  median_onset <- median(onsetsc$firstyes)
  spring <- ifelse(median_onset < 183, "previous", "current")
  summer <- ifelse(median_onset < 274, "previous", "current")

  # Extract weather data for sites, years in onset dataset
  weathersub <- expand.grid(
    site_id = sort(unique(onsetsc$site_id)),
    onset_yr = sort(unique(onsetsc$wateryr)),
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
  onsetsc <- onsetsc %>%
    left_join(weathersub, by = c("site_id", "wateryr" = "onset_yr"))
  # Add firstyes doy (not dowy) 
  onsetsc <- onsetsc %>%
    mutate(firstyes_doy = firstyes - 92)
  
  # Put data in long form, for easy plotting
  onsetscl <- onsetsc %>%
    select(common_name, plantwateryr, wateryr, firstyes, firstyes_doy, 
           contains("spring"), contains("summer"), contains("fall"), 
           contains("winter")) %>%
    pivot_longer(cols = spring_prcp:winter_tmax, 
                 names_to = "weather_var",
                 values_to = "value") %>%
    separate(weather_var, sep = "_", into = c("season", "var")) %>%
    mutate(Var = case_when(
      var == "prcp" ~ "Precipitation (mm)",
      var == "tmax" ~ "Maximum temperature (C)",
      var == "tmin" ~ "Minimum temperature (C)")) %>%
    mutate(Season = factor(str_to_sentence(season),
                           levels = c("Winter", "Spring", "Summer", "Fall")))
  
  # Create scatterplots: matrix of onset date vs weather variable (rows = 
  # season, columns = precip, tmin, tmax). NOTE: these scatterplots and r, P 
  # listed on each plot do not take things like repeated measures or latitude
  # into account.
  # For red maple, want to adjust breaks to capture individuals before Jan 1
  date_breaks <- dowy %>%
    filter(!is.na(dowy)) %>%
    filter(mon %in% seq(1, 12, by = 3))
  
  # Scatterplot with all data
  ggplot(data = onsetscl, aes(x = value, y = firstyes)) +
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x) +
    facet_grid(rows = vars(Season), cols = vars(Var), scales = "free_x") +
    stat_cor(method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01) +
    scale_y_continuous(limits = c(1, 366), 
                       breaks = c(date_breaks$dowy, 366),
                       labels = format(c(date_breaks$first_date, 
                                         as.Date("2022-10-01")), 
                                       "%b")) +
    labs(x = "", y = "Onset date", 
         title = paste0(str_to_sentence(spp), ", ", pheno, " onset date"))

  # Same scatterplot, but highlight data from 2023 in red
  onsetscl <- onsetscl %>%
    mutate(y2023 = ifelse(wateryr == 2023, "2023", paste0(min(yrs), "-2022")),
           y2023 = factor(y2023, levels = c(paste0(min(yrs), "-2022"), "2023")))
  ggplot(data = onsetscl, aes(x = value, y = firstyes)) +
    geom_point(aes(group = y2023, color = y2023)) +
    scale_color_manual(values = c("black","salmon"), name = "Year") +
    stat_smooth(method = "lm", formula = y ~ x) +
    facet_grid(rows = vars(Season), cols = vars(Var), scales = "free_x") +
    stat_cor(method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01) +
    scale_y_continuous(limits = c(1, 366), 
                       breaks = c(date_breaks$dowy, 366),
                       labels = format(c(date_breaks$first_date, 
                                         as.Date("2022-10-01")), 
                                       "%b")) +
    labs(x = "", y = "Onset date", 
         title = paste0(str_to_sentence(spp), ", ", pheno, " onset date"))  

  # Look at relationship between summer max temps and onset date including/
  # excluding 2023 data
  cols <- c("2012-2023" = "darkorchid2", "2012-2022" = "gray20")
  ggplot(data = onsetsc, aes(x = summer_tmax, y = firstyes)) +
    stat_smooth(data = onsetsc,
                method = "lm", formula = y ~ x, 
                aes(col = "2012-2023", fill = "2012-2023"), alpha = 0.2) +
    stat_smooth(data = filter(onsetsc, wateryr < 2023),
                method = "lm", formula = y ~ x, 
                aes(col = "2012-2022", fill = "2012-2022"), alpha = 0.2) +  
    geom_point(data = onsetsc, 
               aes(color = "2012-2023")) +
    geom_point(data = filter(onsetsc, wateryr < 2023),
               aes(color = "2012-2022")) +  
    scale_color_manual(values = cols, name = "Years") +
    scale_fill_manual(values = cols, guide = "none") +
    stat_cor(data = onsetsc,
             aes(color = "2012-2023"),      
             method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01, 
             label.x = 31.25, label.y = 76,
             show.legend = FALSE) +
    stat_cor(data = filter(onsetsc, wateryr < 2023),
             aes(color = "2012-2022"),      
             method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01, 
             label.x = 31.25, label.y = 80,
             show.legend = FALSE) +
    labs(x = "Minimum summer temperature (degC)",
         y = "Open flower onset DOWY") +
    theme_bw()
  
  # Exploring variation in onset dates to identify appropriate random effect
  # structure in linear models
  
  # Quick look to see how much onset dates vary by site
  ggplot(data = onsetsc, aes(x = as.factor(site_id), y = firstyes)) +
    geom_boxplot() +
    geom_jitter(color = "blue", width = 0.1, height = 0)
  
  # Quick look for correlations between onset dates and latitude
  ggplot(data = onsetsc, aes(x = lat, y = firstyes)) +
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x) +
    labs(x = "Latitude", y = "Onset DOWY") +
    stat_cor(method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01)
  
  # Is latitude correlated with weather variables?
  round(cor(select(onsetsc, contains("winter"), contains("spring"),
                   contains("summer"), contains("fall")), onsetsc$lat), 2)
  # Once we're looking at a much smaller area, less correlations
  # Spring, summer, fall tmins are negatively correlated with latitude (-0.62 to -0.75)
  
  # Winter tmin, tmax correlated (0.87)
  # seasonal tmins slightly correlated (~0.60)
  
  # Quick look to see how much onset dates vary by year
  ggplot(data = onsetsc, aes(x = as.factor(wateryr), y = firstyes)) +
    geom_boxplot() +
    geom_jitter(color = "blue", width = 0.1, height = 0)
  
  
  # Run a few linear models with random effects

  seasons <- c("winter", "spring", "summer", "fall")
  vars <- c("prcp", "tmax", "tmin")
  
  # Pick one weather variable for now (spring_tmin)
  seas <- seasons[3]
  var <- vars[2]
  seas_var <- paste0(seas, "_", var)
  
  # Code I could use to run things in a loop over different weather vars...
  # m1formula <- as.formula(paste0("firstyes ~ ", seas_var, " + (1|site_id)"))
  # m1 <- lmer(m1formula, data = onsetsub)
  
  # Using lme4 package so we can create more complex random effect structures
  # (and hardcoding model formulas for now)
  
  # Model with year and site random effects (crossed)
  m.rsiteyr <- lmer(firstyes ~ summer_tmax + (1|site_id) + (1|wateryr),
                    data = onsetsc, REML = TRUE)
  summary(m.rsiteyr)
  confint(m.rsiteyr)
  # Negative effect of spring min temps (though 95% CI does span 0)
  # Massive variation among sites, much less among years
  
  # On a side note: is a model with latitude instead of spring temps better?
  mlat.rsiteyr <- lmer(firstyes_doy ~ lat + (1|site_id) + (1|wateryr),
                       data = onsetsub, REML = TRUE)
  summary(mlat.rsiteyr)
  confint(mlat.rsiteyr)
  anova(mlat.rsiteyr, m.rsiteyr) # Model values almost identical...
  # Maybe not surprising given that so few sites that were observed
  # in multiple years. Can't separate temperature and latitudinal effects.
  # Could really use more data in the SC region and/or more sites observed 
  # in multiple years.
  
  # Model with just site RE
  m.rsite <- lmer(firstyes ~ summer_tmax + (1|site_id),
                  data = onsetsc, REML = TRUE)
  summary(m.rsite)
  anova(m.rsiteyr, m.rsite)
  
  # Maybe don't need year REs

  # Models that include latitude
  mlat.rsite <- lmer(firstyes ~ summer_tmax + lat + (1|site_id),
                     data = onsetsc, REML = TRUE)
  summary(mlat.rsite)
  # No evidence that we need to include latitude
  # Estimation problems with interaction model
  
  # Models with summer tmax, summer precip
  msump.rsite <- lmer(firstyes ~ summer_tmax + summer_prcp + (1|site_id),
                      data = onsetsc, REML = TRUE)
  summary(msump.rsite)
  anova(msump.rsite, m.rsite)
  # Adding summer prcp doesn't help
  
  # Models with summer tmax, winter tmin (cor = 0.57)
  mwinm.rsite <- lmer(firstyes ~ summer_tmax + winter_tmin + (1|site_id),
                      data = onsetsc, REML = TRUE)
  summary(mwinm.rsite)
  anova(mwinm.rsite, m.rsite)
  # Adding winter tmin doesn't help

  # Model diagnostics
  plot(m.rsite) # resid vs fitted values. Slight funnel shape, but not terrible.
  # Makes sense given plot below
  ggplot(onsetsc, aes(x = summer_tmax, y = firstyes)) + 
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x) +
    stat_cor(method = "pearson", cor.coef.name = "r",
             r.accuracy = 0.01, p.accuracy = 0.01, 
             label.x = 31.5, label.y= 75)

  qqnorm(resid(m.rsite))
  qqline(resid(m.rsite))
  # Not too bad
  
  # Other checks with residuals?
  plot(resid(m.rsite) ~ onsetsc$lat) 

  # Do residuals vary with year?
  ggplot(data = onsetsc, aes(x = as.factor(wateryr), 
                              y = resid(m.rsiteyr))) +
    geom_boxplot() +
    geom_jitter(color = "blue", width = 0.1, height = 0)
  