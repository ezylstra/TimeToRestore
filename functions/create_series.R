## TODO: simplify this function because we don't need to download status data
# Input should be previously downloaded data that's been simplified (mostly
# in terms of variable names)


#' Create series dataset
#' 
#' @param start_date date, provided as "YYYY-MM-DD", that indicates the earliest
#' date an output series ends (default = "2009-01-01")
#' @param end_date date, provided as "YYYY-MM-DD", that indicates the latest
#' date an output series begins (default = today's date)
#' @param request_source character string indicating who is requesting the data
#' (required)
#' @param max_yes_gap maximum number of days between consecutive yeses before
#' the 2nd yes is considered the start of a new series (default = 90).
#' @param site_ids vector of NPN site/station IDs. If not specified, will return
#' series associated with all sites.
#' @param species_ids vector of NPN species IDs. If not specified, will return
#' series associated with all species.
#' @param phenophase_ids vector of NPN phenophse IDs. If not specified, will 
#' return series associated with all phenophases.
#'
#' @details 
#' 
#' @return A dataframe, where each row provides information about a series of
#' positive phenophase status observations ("yeses")

create_series <- function(start_date = "2009-01-01",
                          end_date = Sys.Date(),
                          request_source,
                          max_yes_gap = 90,
                          site_ids = NULL, 
                          species_ids = NULL, 
                          phenophase_ids = NULL) {
  
  # Extract the name of this function for reporting
  function_name <- as.character(match.call())[1]
  
  # Libraries required for this function to work
  dependencies <- c("dplyr", "lubridate", "rnpn")
  if (!all(unlist(lapply(X = dependencies, FUN = require, character.only = TRUE)))) {
    stop("At least one package required by ", function_name, 
         " could not be loaded: ", paste(dependencies, collapse = ", "),
         " are required.")
  }
  
  # Make sure dates are appropriate
  start_test <- suppressWarnings(parse_date_time(start_date, orders ="%Y-%m-%d"))
  if (is.na(start_test)) {
    stop(function_name, " requires a valid start date formatted as 'YYYY-MM-DD'")
  }
  end_test <- suppressWarnings(parse_date_time(end_date, orders ="%Y-%m-%d"))
  if (is.na(end_test)) {
    stop(function_name, " requires a valid end date formatted as 'YYYY-MM-DD'")
  }
  if (start_date >= end_date) {
    stop(function_name, " requires end date to be later than start date")
  }
  
  # Make sure request source is specified
  if (is.na(request_source)) {
    stop(function_name, " requires name of request source")
  }
  
  # Make sure max_yes_gap is numeric
  if (!is.numeric(max_yes_gap)) {
    stop(function_name, " requires numeric value for max_yes_gap")
  }
  
  # Extract data from period extending >= 1 year before and after desired period 
  # to get proper calculation of prior and subsequent nos ##### KEY STEP ######
  start_date <- ymd(start_date)
  start_yr <- year(start_date)
  start_extract <- start_yr - 1
  end_date <- ymd(end_date)
  end_yr <- year(end_date)
  end_extract <- end_yr + 1
  
  # Download status-intensity data 
  ##### Using rnpn for now to simplify things, but will need to access Cached 
  ##### Observation table in Database via other means eventually
  si_orig <- npn_download_status_data(
    request_source = request_source,
    years = start_extract:end_extract,
    station_ids = site_ids,
    species_ids = species_ids,
    phenophase_ids = phenophase_ids,
    additional_fields = "observedby_person_id"
  ) %>% data.frame()
  
  # Remove unnecssary columns for now and rename a few columns
  si <- si_orig %>%
    select(-c(update_datetime, abundance_value, intensity_category_id, 
              intensity_value)) %>%
    rename(person = observedby_person_id,
           elevation_m = elevation_in_meters,
           doy = day_of_year,
           status = phenophase_status)
  
  # Arrange by individual, date
  si <- si %>%
    arrange(common_name, species_id, site_id, individual_id, observation_date, 
            doy, phenophase_id, status)
  
  # Discard observations with unknown status (-1)
  ser <- si %>%
    filter(status >= 0)
  
  # Extract information about each site (to append again later)
  sites <- ser %>%
    distinct(site_id, latitude, longitude, elevation_m, state)
  
  # Extract information about each species (to append again later)
  spps <- ser %>%
    distinct(species_id, genus, species, common_name, kingdom)
  
  # Extract information about phenophases (to append again later)
  phps <- ser %>%
    distinct(phenophase_id, phenophase_description)
  
  # Now, simplify main dataframe by removing this information
  ser <- ser %>%
    select(-c(latitude, longitude, elevation_m, state, 
              genus, species, common_name, kingdom,
              phenophase_description))
  
  # Summarize data for each individual, phenophase, and date; denote where
  # we have multiple observers and flag status conflicts
  ser <- ser %>%
    group_by(site_id, species_id, individual_id, phenophase_id,
             observation_date, doy) %>%
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
  # fake <- data.frame(site_id = 24702,
  #                    species_id = 210,
  #                    individual_id = 1111111,
  #                    phenophase_id = 500,
  #                    observation_date = ymd(
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
    distinct(site_id, species_id, individual_id, phenophase_id)
  
  # Loop through each individual-phenophase combination
  for (i in 1:nrow(combos)) {
    ser1 <- ser %>%
      filter(individual_id == combos$individual_id[i] & 
               phenophase_id == combos$phenophase_id[i])
    
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
        gaps <- as.numeric(ser1$observation_date[(ys$startrow[j] + 1):ys$endrow[j]] -
                             ser1$observation_date[ys$startrow[j]:(ys$endrow[j] - 1)]) 
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
          lastnorow = c(ys$lastnorow[j], NA),
          nextnorow = c(NA, ys$nextnorow[j])
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
    
    yesseries$first_yes_date <- ser1$observation_date[yesseries$startrow]
    yesseries$last_yes_date <- ser1$observation_date[yesseries$endrow] 
    
    # If there are no prior nos for any yes series, make all prior no dates NA.
    # Otherwise, add in date of prior no
    if (sum(is.na(yesseries$lastnorow)) == nrow(yesseries)) {
      yesseries$prior_no_date <- NA
    } else {
      yesseries$prior_no_date <- ser1$observation_date[yesseries$lastnorow]
    }
    # If there are no next nos for any yes series, make all next no dates NA.
    # Otherwise, add in date of next No
    if (sum(is.na(yesseries$nextnorow)) == nrow(yesseries)) {
      yesseries$next_no_date <- NA
    } else {
      yesseries$next_no_date <- ser1$observation_date[yesseries$nextnorow]
    }
    
    # Append information about the series
    yesseries$site_id <- combos$site_id[i]
    yesseries$species_id <- combos$species_id[i]
    yesseries$individual_id <- combos$individual_id[i]
    yesseries$phenophase_id <- combos$phenophase_id[i]
    yesseries <- yesseries %>%
      select(site_id, species_id, individual_id, phenophase_id, first_yes_date, 
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
    
    # Only keep series that have first or last yes within user-defined 
    # start/end dates
    series <- series %>%
      filter((first_yes_date >= start_date & first_yes_date <= end_date) |
               (last_yes_date >= start_date & last_yes_date <= end_date))
    
    if (nrow(series) == 0) {
      
      stop("No series for the selected species and phenophases have a first ",
           "or last yes within user-selected start and end dates")
      
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
          filter(individual_id == series$individual_id[i] & 
                   phenophase_id == series$phenophase_id[i]) %>%
          filter(observation_date >= series$first_yes_date[i] &
                   observation_date <= series$last_yes_date[i])
        
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
        mutate(days_prior_no = as.numeric(first_yes_date - prior_no_date),
               days_next_no = as.numeric(next_no_date - last_yes_date),
               series_days = as.numeric(last_yes_date - first_yes_date) + 1)
      
      # Add in year, DOY, and Julian dates for first, last yes
      series <- series %>%
        mutate(first_yes_year = year(first_yes_date),
               first_yes_doy = yday(first_yes_date),
               first_yes_julian = julian_date(first_yes_date),
               last_yes_year = year(last_yes_date),
               last_yes_doy = yday(last_yes_date),
               last_yes_julian = julian_date(last_yes_date)) 
      
      # Add site, species, phenophase information back in and arrange columns
      series <- series %>%
        left_join(sites, by = "site_id") %>%
        left_join(spps, by = "species_id") %>%
        left_join(phps, by = "phenophase_id")
      
      series <- series %>%      
        select(site_id, latitude, longitude, elevation_m, state, species_id,
               genus, species, common_name, kingdom, individual_id, 
               phenophase_id, phenophase_description, first_yes_date,
               first_yes_year, first_yes_julian, prior_no_date, days_prior_no,
               last_yes_date, last_yes_year, last_yes_julian, next_no_date,
               days_next_no, series_yeses, series_days, multiple_observers,
               person_id, status_conflict_flag, series_split_flag)
    }
  }
  
  return(series)
}
