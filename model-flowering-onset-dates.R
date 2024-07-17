################################################################################
# Start developing linear models to explore climate effects on flowering
# phenology metrics

# Some of this code is derived from that originally developed by Alyssa 
# Rosemartin, Hayley Limes, and Jeff Oliver. 
# See: https://github.com/alyssarosemartin/time-to-restore

# Erin Zylstra
# 2024-07-16
################################################################################

require(dplyr)
require(lubridate)
require(stringr)
require(tidyr)
require(ggplot2)
require(cowplot)
require(terra)
require(tidyterra)
require(pdftools)

rm(list = ls())

# Note: throughout script fl = flower phenophase; fo = open flower phenophase

# Load shapefile with state boundaries ----------------------------------------#

  states <- vect("data/states/cb_2017_us_state_500k.shp")
  states <- subset(states, 
                   !states$STUSPS %in% c("HI", "AK", "VI", "MP", 
                                         "GU", "PR", "AS"))
  states$sc <- as.factor(ifelse(states$STUSPS %in% c("LA", "NM", "OK", "TX"), 
                                1, 0))  

# Load species info -----------------------------------------------------------#

  spp <- read.csv("data/priority-species.csv")
  spp <- spp %>%
    mutate(across(starts_with("votes_"), ~ replace_na(., 0))) %>%
    mutate(across(starts_with("revisit_"), ~ replace_na(., 0))) %>%
    mutate(LA = if_else(votes_LA > 0 | revisit_votes_LA > 0, 1, 0),
           NM = if_else(votes_NM > 0 | revisit_votes_NM > 0, 1, 0),
           OK = if_else(votes_OK > 0 | revisit_votes_OK > 0, 1, 0),
           TX = if_else(revisit_votes_TX > 0, 1, 0)) %>%
    select(common_name, genus, species, species_id, priority, LA, NM, OK, TX)

# Load processed NPN status/intensity data and format -------------------------#

  df <- read.csv("data/flower-status-intensities-priorityspp.csv")
  
  # Rename/remove/format columns
  df <- df %>%
    rename(lat = latitude,
           lon = longitude,
           plant_id = individual_id) %>%
    mutate(observation_date = ymd(observation_date)) %>%
    select(-c(person_id, n_observations, species_id, genus, species,
              midpoint_fl, midpoint_fo, num_open))

  # Identify the water year for each observation since some species (like red 
  # maple) can flower at the start of the calendar year. Then, calculate
  # day-of-wateryr (dowy)
  df <- df %>%
    mutate(mon = month(observation_date)) %>%
    mutate(wateryr = ifelse(mon %in% 10:12, year + 1, year)) %>%
    mutate(wateryrday1 = make_date(year = wateryr - 1, month = 9, day = 30)) %>%
    mutate(dowy = as.numeric(observation_date - wateryrday1)) %>%
    select(-c(mon, wateryrday1))
    
  # Create a unique ID for each plant-wateryr to make various matches easier
  df <- df %>%
    mutate(plantwateryr = paste0(plant_id, "_", wateryr))
  
# Summarize data for each plant-wateryr ---------------------------------------#
  
  # Just looking at onset of flowering and open-flower phenophases. 
  # Might see whether we have sufficient data to model phenophase duration
  # at a later date.
  
  # Create a new dataframe to hold plant-wateryr stats 
  py <- df %>%
    group_by(plantwateryr, plant_id, wateryr, common_name, state, lat, lon) %>%
    summarize(n_obs = n(),
              n_status_fl = sum(!is.na(status_fl)),
              n_status_fo = sum(!is.na(status_fo)),
              n_fl_yes = sum(!is.na(status_fl) & status_fl == 1),
              n_fo_yes = sum(!is.na(status_fo) & status_fo == 1),
              first_obs = min(dowy),
              last_obs = max(dowy),
              .groups = "keep") %>%
    data.frame()
  
  # Remove plant-wateryrs where the plant was observed only once or where
  # observers never confirmed flowers were present. (Note: this removes ~48%
  # of plant-wateryrs)
  py <- py %>%
    filter(n_obs > 1 & n_fl_yes > 0)
  
  # Add columns for summary stats
  py$fl_firstyes <- NA
  py$fl_dayssince_lastno <- NA
  py$fo_firstyes <- NA
  py$fo_dayssince_lastno <- NA
  
  for (i in 1:nrow(py)) {

    # Extract data for selected plant-wateryr and remove dates when flower
    # status wasn't recorded
    dffl <- df %>%
      select(plantwateryr, status_fl, dowy) %>%
      filter(plantwateryr == py$plantwateryr[i]) %>%
      filter(!is.na(status_fl))
    
    # Identify rows of dffl with first flowers = "yes" and first flowers = "no"
    first_yesfl_ind <- first(which(dffl$status_fl == 1))
    first_nofl_ind <- first(which(dffl$status_fl == 0))
    # Identify first day of wateryr with flowers confirmed present
    py$fl_firstyes[i] <- dffl$dowy[first_yesfl_ind]
    # Calculate number of days since previous observation without flowers. If 
    # there's no preceding observation with flowers = "no", then 
    # fl_dayssince_lastno = NA.
    if (first_yesfl_ind > first_nofl_ind & !is.na(first_nofl_ind)) {
      py$fl_dayssince_lastno[i] <- 
        py$fl_firstyes[i] - dffl$dowy[first_yesfl_ind - 1]
    } else {
      py$fl_dayssince_lastno[i] <- NA
    }
    
    # If there were no observations with open flowers = "yes" during
    # plant-wateryr, skip the rest of the loop
    if (py$n_fo_yes[i] == 0) {next}
    
    # Extract data for selected plant-wateryr and remove dates when open flower
    # status wasn't recorded
    dffo <- df %>%
      select(plantwateryr, status_fo, dowy) %>%
      filter(plantwateryr == py$plantwateryr[i]) %>%
      filter(!is.na(status_fo))
    
    # Identify rows of dffo with first open flowers = "yes" and first open 
    # flowers = "no"
    first_yesfo_ind <- first(which(dffo$status_fo == 1))
    first_nofo_ind <- first(which(dffo$status_fo == 0))
    # Identify first day of wateryr with open flowers confirmed present
    py$fo_firstyes[i] <- dffo$dowy[first_yesfo_ind]
    # Calculate number of days since previous observation without open flowers. 
    # If there's no preceding observation with poen flowers = "no", then 
    # fo_dayssince_lastno = NA.
    if (first_yesfo_ind > first_nofo_ind & !is.na(first_nofo_ind)) {
      py$fo_dayssince_lastno[i] <- 
        py$fo_firstyes[i] - dffo$dowy[first_yesfo_ind - 1]
    } else {
      py$fo_dayssince_lastno[i] <- NA
    }
  }

  # Look at distribution of values (across species):
  # par(mfrow = c(2, 1))
  # hist(py$fl_firstyes, breaks = 50, xlim = c(1, 365))
  # hist(py$fo_firstyes, breaks = 50, xlim = c(1, 365))

# Assess amount of onset data available for each species ----------------------#  

  # Set the maximum number of days prior to a "yes" that we'll need a "no" in 
  # order to use the observation in models for onset
  daysprior <- 14  
  
  # Flowers
  py_fl <- py %>%
    filter(!is.na(fl_dayssince_lastno)) %>%
    filter(fl_dayssince_lastno <= daysprior)
  dim(py_fl) # 2819 plantyrs
  
  fl_spp <- count(py_fl, common_name) %>%
    rename(n_plantyrs = n) %>%
    left_join(dplyr::select(spp, -c(genus, species, species_id)), 
              by = "common_name") %>%
    arrange(desc(n_plantyrs))
  fl_spp
  # 31 species: 11 that have >30 plantyrs

  # Open flowers
  py_fo <- py %>%
    filter(!is.na(fo_dayssince_lastno)) %>%
    filter(fo_dayssince_lastno <= daysprior)
  dim(py_fo) # 3829 plantyrs
  
  fo_spp <- count(py_fo, common_name) %>%
    rename(n_plantyrs = n) %>%
    left_join(select(spp, -c(genus, species, species_id)), 
              by = "common_name") %>%
    arrange(desc(n_plantyrs))
  fo_spp
  # 33 species: 12 that have >30 plantyrs
  
  # Note that we have more data to evaluate climate effects on the onset of 
  # the open flower phenophase than we do the flower phenophase. Surprising.
  
# Map locations of plants with onset data -------------------------------------#
  
  # Set minimum number of plant-yrs to move forward with species
  min_plantyrs <- 30
  
  # Get common names of species with sufficient data
  fl_common_names <- fl_spp$common_name[fl_spp$n_plantyrs >= min_plantyrs]
  fo_common_names <- fo_spp$common_name[fo_spp$n_plantyrs >= min_plantyrs]
  
  # Create map for each species with locations of plants that have one or more
  # years of flowering onset data
  fl_ggplots <- list()
    for (i in 1:length(fl_common_names)) {

    # Extract onset data
    py_fl_spp <- filter(py_fl, common_name == fl_common_names[i])
    
    # Isolate location information
    fl_sites <- py_fl_spp %>%
      group_by(lat, lon, state) %>%
      summarize(n_plantyrs = n(), .groups = "keep") %>%
      data.frame()
    fl_sitesv <- vect(fl_sites, geom = c("lon", "lat"), crs = "epsg:4269")
    
    # Create map with locations of plants with onset data
    spp_title <- paste0("Onset, flowering: ", 
                        str_to_sentence(fl_common_names[i]), 
                        " (", nrow(py_fl_spp), " plant-years)")
    spp_fl_map <- ggplot(data = states) +
      geom_spatvector(aes(fill = sc), color = "gray60") +
      scale_fill_discrete(type = c("white", "gray95"), guide = "none") +
      geom_spatvector(data = fl_sitesv, aes(size = n_plantyrs), 
                      col = "blue", pch = 1) +
      labs(size = "No. plant-years") +
      labs(title = spp_title) +
      theme(panel.background = element_rect(fill = rgb(t(col2rgb("lightblue")), 
                                                       alpha = 0.2 * 255, 
                                                       maxColorValue = 255)))
    fl_ggplots[[i]] <- spp_fl_map
  }
  
  # Combine flowering onset maps into a single pdf and save in output folder
  pages <- list()
  for (j in 1:ceiling(length(fl_ggplots) / 3)) {
  
    i <- (j - 1) * 3 + 1
    
    plot1 <- plot_grid(NULL, fl_ggplots[[i]], NULL, nrow = 1,
                       rel_widths = c(0.5, 7.5, 0.5))
    
    if (i + 1 > length(fl_ggplots)) {
      plot2 <- NULL
    } else {
      plot2 <- plot_grid(NULL, fl_ggplots[[i + 1]], NULL, nrow = 1,
                         rel_widths = c(0.5, 7.5, 0.5))
    }
    
    if (i + 2 > length(fl_ggplots)) {
      plot3 <- NULL
    } else {
      plot3 <- plot_grid(NULL, fl_ggplots[[i + 2]], NULL, nrow = 1,
                         rel_widths = c(0.5, 7.5, 0.5))
    }
    
    pages[[j]] <- plot_grid(NULL, plot1, NULL, plot2, NULL, plot3, NULL,
                            ncol = 1, 
                            rel_heights = c(0.6, 3.2, 0.1, 3.2, 0.1, 3.2, 0.6), 
                            align = "v")
  }
  
  for (j in 1:length(pages)) {
    pdf_name <- paste0("output/maps/flowering-onset-map-", j, ".pdf")
    ggsave(filename = pdf_name, 
           plot = pages[[j]],
           width = 8.5, 
           height = 11, 
           units = "in", 
           device = cairo_pdf)
  }
  
  pdf_name <- "output/maps/flowering-onset-maps.pdf"
  ind_map_pages <- list.files(path = "output/maps/",
                              pattern = "flowering-onset-map-",
                              full.names = TRUE)
  invisible(pdf_combine(ind_map_pages, output = pdf_name))
  invisible(file.remove(ind_map_pages))
  
  # Create map for each species with locations of plants that have one or more
  # years of open flower onset data
  fo_ggplots <- list()
  for (i in 1:length(fo_common_names)) {
    
    # Extract onset data
    py_fo_spp <- filter(py_fo, common_name == fo_common_names[i])
    
    # Isolate location information
    fo_sites <- py_fo_spp %>%
      group_by(lat, lon, state) %>%
      summarize(n_plantyrs = n(), .groups = "keep") %>%
      data.frame()
    fo_sitesv <- vect(fo_sites, geom = c("lon", "lat"), crs = "epsg:4269")
    
    # Create map with locations of plants with onset data
    spp_title <- paste0("Onset, open flowers: ", 
                        str_to_sentence(fo_common_names[i]), 
                        " (", nrow(py_fo_spp), " plant-years)")
    spp_fo_map <- ggplot(data = states) +
      geom_spatvector(aes(fill = sc), color = "gray60") +
      scale_fill_discrete(type = c("white", "gray95"), guide = "none") +
      geom_spatvector(data = fo_sitesv, aes(size = n_plantyrs), 
                      col = "blue", pch = 1) +
      labs(size = "No. plant-years") +
      labs(title = spp_title) +
      theme(panel.background = element_rect(fill = rgb(t(col2rgb("lightblue")), 
                                                       alpha = 0.2 * 255, 
                                                       maxColorValue = 255)))
    fo_ggplots[[i]] <- spp_fo_map
  }
  
  # Combine flowering onset maps into a single pdf and save in output folder
  pages <- list()
  for (j in 1:ceiling(length(fo_ggplots) / 3)) {
    
    i <- (j - 1) * 3 + 1
    
    plot1 <- plot_grid(NULL, fo_ggplots[[i]], NULL, nrow = 1,
                       rel_widths = c(0.5, 7.5, 0.5))
    
    if (i + 1 > length(fo_ggplots)) {
      plot2 <- NULL
    } else {
      plot2 <- plot_grid(NULL, fo_ggplots[[i + 1]], NULL, nrow = 1,
                         rel_widths = c(0.5, 7.5, 0.5))
    }
    
    if (i + 2 > length(fo_ggplots)) {
      plot3 <- NULL
    } else {
      plot3 <- plot_grid(NULL, fo_ggplots[[i + 2]], NULL, nrow = 1,
                         rel_widths = c(0.5, 7.5, 0.5))
    }
    
    pages[[j]] <- plot_grid(NULL, plot1, NULL, plot2, NULL, plot3, NULL,
                            ncol = 1, 
                            rel_heights = c(0.6, 3.2, 0.1, 3.2, 0.1, 3.2, 0.6), 
                            align = "v")
  }
  
  for (j in 1:length(pages)) {
    pdf_name <- paste0("output/maps/open-flower-onset-map-", j, ".pdf")
    ggsave(filename = pdf_name, 
           plot = pages[[j]],
           width = 8.5, 
           height = 11, 
           units = "in", 
           device = cairo_pdf)
  }
  
  pdf_name <- "output/maps/open-flower-onset-maps.pdf"
  ind_map_pages <- list.files(path = "output/maps/",
                              pattern = "open-flower-onset-map-",
                              full.names = TRUE)
  invisible(pdf_combine(ind_map_pages, output = pdf_name))
  invisible(file.remove(ind_map_pages))
