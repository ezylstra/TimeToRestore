#' Create series from (previously obtained) status-intensity data (will create
#' data for all sites/species/phenophases included in SI dataset). Can specify
#' time period or default to time span covered by SI dataset by setting
#' start_date and end_date to NA
#' 
#' @param si_data dataframe with status-intensity data (and simplified column
#' names)
#' @param start_date date, provided as "YYYY-MM-DD", that indicates the earliest
#' date an output series ends. If start_date or end_date set to NA (default), 
#' function will generate all series within timeframe of status-intensity data.
#' @param end_date date, provided as "YYYY-MM-DD", that indicates the latest
#' date an output series begins (default = today's date)
#' @param max_yes_gap maximum number of days between consecutive yeses before
#' the 2nd yes is considered the start of a new series (default = 90).
#' 
#' @return A dataframe, where each row provides information about a series of
#' positive phenophase status observations ("yeses")

create_series <- function(si_data = df,
                          max_yes_gap = 90) {

  # Remove unnecssary columns
  si <- si_data %>%
    select(-c(intensity_value, midpoint))
  
  # Arrange by individual, date
  si <- si %>%
    arrange(spp, site, id, obsdate, php_id, status)
  
  # Discard observations with unknown status (-1), if present
  ser <- si %>%
    filter(status >= 0)
  
  # Extract information about each site (to append again later)
  sites <- ser %>%
    distinct(site, site_name, lpp, lat, lon, state, region, state4, state5)
  
  # Extract information about each species (to append again later)
  spps <- ser %>%
    distinct(spp_id, spp)
  
  # Extract information about phenophases (to append again later)
  phps <- ser %>%
    distinct(php_id, php)
  
  # Now, simplify main dataframe by removing this information
  ser <- ser %>%
    select(-c(site_name, lpp, lat, lon, state, region, state4, state5, 
              spp, php))
  
  # Summarize data for each individual, phenophase, and date; denote where
  # we have multiple observers and flag status conflicts
  ser <- ser %>%
    group_by(site, spp_id, id, php_id, obsdate, doy) %>%
    summarize(mult_observer = ifelse(n_distinct(person) > 1, 1, 0),
              person_id = paste(unique(person), collapse = ","),
              status_conflict = ifelse(n_distinct(status) == 2, 1, 0),
              status = max(status),
              .groups = "drop") %>%
    mutate(status_conflict_flag = case_when(
      status_conflict == 1 & mult_observer == 1 ~ "MultiObserver-StatusConflict",
      status_conflict == 1 & mult_observer == 0 ~ "OneObserver-StatusConflict",
      status_conflict == 0 ~ NA
    )) %>%
    data.frame()
  
  # # Add some fake data to test max gap thing...
  # fake <- data.frame(site = 24702,
  #                    spp_id = 210,
  #                    id = 1111111,
  #                    php_id = 500,
  #                    obsdate = ymd(
  #                      c("2021-02-01", "2021-02-02", "2021-02-05",
  #                        "2021-02-06", "2021-08-01", "2021-08-03",
  #                        "2021-08-05", "2021-10-01")),
  #                    doy = yday(c("2021-02-01", "2021-02-02", "2021-02-05",
  #                                 "2021-02-06", "2021-08-01", "2021-08-03",
  #                                 "2021-08-05", "2021-10-01")),
  #                    mult_observer = 0,
  #                    person_id = 99999,
  #                    status_conflict = 0,
  #                    status = c(0, rep(1, 5), 0, 1),
  #                    status_conflict_flag = 0)
  # ser <- rbind(fake, ser)
  
  # Identify unique combinations of individual plant and phenophase
  combos <- ser %>%
    distinct(site, spp_id, id, php_id)
  
  # Loop through each individual-phenophase combination
  for (i in 1:nrow(combos)) {
    ser1 <- ser %>%
      filter(id == combos$id[i] & php_id == combos$php_id[i])
    
    rles <- rle(ser1$status)
    
    # If there are no series of yeses, skip to next combination
    nseries <- sum(rles$value == 1)
    if (nseries == 0) {next}
    
    # Get row number (index) for the start of each 0 or 1 series
    starts <- cumsum(c(1, rles$lengths))
    starts <- starts[-length(starts)]
    
    # Extract table with information about each run of 0s and 1s
    series01 <- data.frame(value = rles$values,
                           length = rles$lengths,
                           startrow = starts) %>%
      mutate(endrow = length + startrow - 1)
    # Extract table with information about each run of 1s (series)
    yesseries <- series01 %>% filter(value == 1) %>%
      rename(series_yeses = length) %>%
      select(-value) %>%
      mutate(lastnorow = ifelse(startrow == 1, NA, startrow - 1)) %>%
      mutate(nextnorow = ifelse(endrow == nrow(ser1), NA, endrow + 1))
    
    # Identify any yesseries that have problems with max_yes_gap
    ys <- yesseries
    for (j in 1:nrow(ys)) {
      
      # If there's only one yes in a series, then there is no gap (and no
      # problem). But if there's more than one yes, extract the gaps between
      # consecutive yeses to see if they exceed the user-defined limit
      if (ys$series_yeses[j] == 1) {
        gaps <- 0
      } else {
        gaps <- as.numeric(ser1$obsdate[(ys$startrow[j] + 1):ys$endrow[j]] -
                             ser1$obsdate[ys$startrow[j]:(ys$endrow[j] - 1)]) 
      }
      # Identify new series that need to be created
      newseries <- which(gaps > max_yes_gap)
      
      # If there were consecutive yeses separated by more than the user-defined
      # maximum (max_yes_gap), then create a new series starting at the 2nd yes
      if (length(newseries) == 0) {
        yesseries_new <- ys[j,]
        yesseries_new$series_split_flag <- 0
      } else {
        yesseries_new <- data.frame(
          startrow = c(ys$startrow[j], ys$startrow[j] + newseries),
          endrow = c(ys$startrow[j] + newseries - 1, ys$endrow[j]),
          lastnorow = c(ys$lastnorow[j], rep(NA, length(newseries))),
          nextnorow = c(rep(NA, length(newseries)), ys$nextnorow[j])
        ) %>%
          mutate(series_yeses = endrow - startrow + 1, .before = startrow) %>%
          mutate(series_split_flag = 1)
      }
      
      if (j == 1) {
        yesseries <- yesseries_new
      } else {
        yesseries <- rbind(yesseries, yesseries_new)
      }
    }
    
    yesseries$first_yes_date <- ser1$obsdate[yesseries$startrow]
    yesseries$last_yes_date <- ser1$obsdate[yesseries$endrow] 
    
    # If there are no prior nos for any yes series, make all prior no dates NA.
    # Otherwise, add in date of prior no
    if (sum(is.na(yesseries$lastnorow)) == nrow(yesseries)) {
      yesseries$prior_no_date <- NA
    } else {
      yesseries$prior_no_date <- ser1$obsdate[yesseries$lastnorow]
    }
    # If there are no next nos for any yes series, make all next no dates NA.
    # Otherwise, add in date of next No
    if (sum(is.na(yesseries$nextnorow)) == nrow(yesseries)) {
      yesseries$next_no_date <- NA
    } else {
      yesseries$next_no_date <- ser1$obsdate[yesseries$nextnorow]
    }
    
    # Append information about the series
    yesseries$site <- combos$site[i]
    yesseries$spp_id <- combos$spp_id[i]
    yesseries$id <- combos$id[i]
    yesseries$php_id <- combos$php_id[i]
    yesseries <- yesseries %>%
      select(site, spp_id, id, php_id, first_yes_date, 
             prior_no_date, series_yeses, last_yes_date, next_no_date,
             series_split_flag)
    
    if (i == 1) {
      series <- yesseries
    } else {
      series <- rbind(series, yesseries)
    }
  }
  
  # Only proceed if there are one or more series....
  if (!exists("series")) {
    
    stop("No series for the selected species and phenophases")
    
  } else {
    
    # Append indicators/flags to the series dataset....
    # If >1 observer contributed yeses to the series, then multiple_observers = 1
    # If one or more yeses in a series had a status conflict, then flag
    
    flag1 <- ser %>%
      filter(status == 1)
    
    series <- series %>%
      mutate(multiple_observers = NA,
             person_id = NA,
             status_conflict_flag = NA)
    
    for (i in 1:nrow(series)) {
      # Find observation dates in flag1 that fall in series
      flag1sub <- flag1 %>%
        filter(id == series$id[i] & 
                 php_id == series$php_id[i]) %>%
        filter(obsdate >= series$first_yes_date[i] &
                 obsdate <= series$last_yes_date[i])
      
      observers <- paste0(flag1sub$person_id, collapse = ",")
      observers <- sort(unique(strsplit(observers, ",")[[1]]))
      series$multiple_observers[i] <- ifelse(length(observers) == 1, 0, 1)
      series$person_id[i] <- paste0(observers, collapse = ",")
      
      if (all(is.na(flag1sub$status_conflict_flag))) {
        series$status_conflict_flag[i] <- NA
      } else {
        flags <- sort(unique(flag1sub$status_conflict_flag))
        series$status_conflict_flag[i] <- paste0(flags, collapse = ",")
      }
    }
    
    # Calculate days since prior no, until next no, and series days
    series <- series %>%
      mutate(days_prior_no = as.numeric(first_yes_date - as.Date(prior_no_date)),
             days_next_no = as.numeric(as.Date(next_no_date) - last_yes_date),
             series_days = as.numeric(last_yes_date - first_yes_date) + 1)
    
    # Add in year, DOY, and Julian dates for first, last yes
    series <- series %>%
      mutate(first_yes_year = year(first_yes_date),
             first_yes_doy = yday(first_yes_date),
             first_yes_julian = julian_date(first_yes_date)) 
    
    # Add site, species, phenophase information back in and arrange columns
    series <- series %>%
      left_join(sites, by = "site") %>%
      left_join(spps, by = "spp_id") %>%
      left_join(phps, by = "php_id")
    
    series <- series %>%      
      select(site, site_name, lpp, lat, lon, state, region, state4, state5,
             spp_id, spp, id,
             php_id, php, first_yes_date, first_yes_doy,
             first_yes_year, first_yes_julian, days_prior_no,
             last_yes_date, series_yeses, series_days, multiple_observers,
             person_id, status_conflict_flag, series_split_flag)
  }
  
  return(series)
}
