################################################################################
# Download daymet data and summarize by season
# Script is typically called from model-flowering-onset-dates.R

# Erin Zylstra
# 2024-07-23
################################################################################

require(daymetr)

# Create filenames ------------------------------------------------------------#

  daymet_csv <- paste0("data/climate/daymet-", nsites, "sites-", 
                       first(yrs), "-", last(yrs), ".csv")
  daymet_zip <- str_replace(daymet_csv, ".csv", ".zip")

# Download data ---------------------------------------------------------------#
  
  get_daymet <- function(i) {
    
    temp_lat <- sites[i, ] %>% pull(lat)
    temp_lon <- sites[i, ] %>% pull(lon)
    temp_site <- sites[i, ] %>% pull(site_id)
    
    temp_daymet <- download_daymet(
      lat = temp_lat,
      lon = temp_lon, 
      start = first(yrs),
      end = last(yrs)
    ) %>%
      #--- just get the data part ---#
      .$data %>% 
      #--- assign site_id so we can match with onset data ---#
      mutate(site_id = temp_site) %>%
      #--- get date from day of the year ---#
      mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"))
    
    return(temp_daymet)
  }
  
  daymet_df <- lapply(1:nsites, get_daymet) %>%
    bind_rows()
  
  daymet_df <- daymet_df %>%
    rename(prcp_mm = prcp..mm.day.,
           tmin_c = tmin..deg.c.,
           tmax_c = tmax..deg.c.,
           daylength_s = dayl..s.,
           shortwave_radiation_Wm2 = srad..W.m.2.,
           swe_kgm2 = swe..kg.m.2.,
           vapor_press_Pa = vp..Pa.)
  
  # Save original csv and put in a zip file. Remove csv because it's huge.
  write.csv(daymet_df, 
            file = daymet_csv,
            row.names = FALSE)
  zip(daymet_zip, files = daymet_csv)
  file.remove(daymet_csv)

# Create seasonal summaries ---------------------------------------------------#
  
  # Remove climate data we're not interested in
  daymet_d <- daymet_df %>%
    select(-c(daylength_s, shortwave_radiation_Wm2, swe_kgm2, vapor_press_Pa))

  # Will summarize data for each season (as defined by NPN)
    # Spring = Mar-May
    # Summer = Jun-Aug
    # Fall = Sep-Nov
    # Winter = Dec-Feb (assigned to year for Jan-Feb)
  
  # Assign "seasonyr" (so December is associated with following year)
  daymet_d <- daymet_d %>%
    mutate(date = ymd(date)) %>%
    mutate(mon = month(date)) %>%
    mutate(seasonyr = if_else(mon == 12, year + 1, year))
  
  # Remove winter seasons where we don't have data from all months (we won't 
  # need these seasons in models either), and add a season label
  daymet_d <- daymet_d %>%
    filter(seasonyr %in% yrs) %>%
    filter(!(seasonyr == min(seasonyr) & mon %in% 1:2)) %>%
    mutate(season = case_when(
      mon %in% 3:5 ~ "spring",
      mon %in% 6:8 ~ "summer",
      mon %in% 9:11 ~ "fall",
      .default = "winter"
    )) 
  
  # Summarize data by season (accumulated precip, mean tmin and tmax)
  weather <- daymet_d %>%
    group_by(site_id, seasonyr, season) %>%
    summarize(prcp = sum(prcp_mm),
              tmin = mean(tmin_c),
              tmax = mean(tmax_c),
              .groups = "keep") %>%
    mutate(season = factor(season, 
                           levels = c("winter", "spring", "summer", "fall"))) %>%
    data.frame()
  
  # Save seasonal data
  write.csv(weather, 
            file = weather_filename,
            row.names = FALSE)
  
  rm(daymet_d, daymet_df)

# View weather data (optional) ------------------------------------------------#  

  # Quick look at summaries, by season
  
  # ggplot(data = weather) +
  #   geom_histogram(aes(prcp), fill = "steelblue3") +
  #   facet_grid(rows = vars(season)) +
  #   labs(x = "Accumulated precipitation (mm)")
  # ggplot(data = weather) +
  #   geom_histogram(aes(tmin), fill = "steelblue3") +
  #   facet_grid(rows = vars(season)) +
  #   labs(x = "Average minimum temperature (degC)")
  # ggplot(data = weather) +
  #   geom_histogram(aes(tmax), fill = "steelblue3") +
  #   facet_grid(rows = vars(season)) +
  #   labs(x = "Average maximum temperature (degC)")
  
  # Could also look at summaries, by year
  
  # ggplot(data = filter(weather, season == "winter")) +
  #   geom_histogram(aes(tmin), fill = "steelblue3") +
  #   facet_wrap(vars(seasonyr)) +
  #   labs(x = "Average minimum temperature (degC)")
