################################################################################
# Extract and summarize the amount of data available on flower/open flower onset

# Erin Zylstra
# 2024-07-23
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
require(pdftools)

# Note: throughout script fl = flower phenophase; fo = open flower phenophase

# Logical indicating whether to replace csv/pdf files if they already exist ---#

  replace <- TRUE

# Set parameters for pdf output -----------------------------------------------#

  pdfw <- 8.5 # Width of pdf output
  pdfh <- 11  # Height of pdf output

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
  
  # Just looking at onset of flower and open-flower phenophases. Might see 
  # whether we have sufficient data to model phenophase duration later on.

  # Putting all information available for a given plant-wateryr into a single
  # row because we might want to assess how often we're collecting both
  # measures. Plus, this better matches how we've stored status-intensity data.
  
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
  
  # Loop through each plant-wateryr (note that this takes a minute or so)
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

# Save flowering onset data to file -------------------------------------------#

  # First remove observations from wateryr 2024, since we only have data through
  # Dec 2023
  py <- filter(py, wateryr < 2024)
  
  # Extract year range
  wateryr_min <- min(py$wateryr)
  wateryr_max <- max(py$wateryr)
  
  data_filename <- paste0("data/flowering-onsets-", wateryr_min, "-",
                          wateryr_max, ".csv")
  if (!file.exists(data_filename) | 
      (file.exists(data_filename) & replace == TRUE)) {
    write.csv(py, 
              file = data_filename,
              row.names = FALSE)
  }

# Select parameters for data inclusion/exclusion ------------------------------#
  
  # Set the maximum number of days prior to a "yes" that we'll need a "no" in 
  # order to use the observation in onset models
  daysprior <- 14    

  # Set the minimum number of plant-yrs needed to include a species in 
  # phenophase summaries
  min_plantyrs <- 30
  
# Summarize flower/open flower onset data available for each species ----------#
  
  # Extract plant-wateryr data for flower onset
  py_fl <- py %>%
    select(-c(contains("fo"), first_obs, last_obs)) %>%
    rename_with(str_replace, pattern = "fl_", replacement = "") %>%
    rename_with(str_replace, pattern = "_fl", replacement = "") %>%
    mutate(phenophase = "flower") %>%
    relocate(phenophase, .after = "n_obs")
  
  # Extract plant-wateryr data for open flower onset
  py_fo <- py %>%
    select(-c(contains("fl"), first_obs, last_obs)) %>%
    rename_with(str_replace, pattern = "fo_", replacement = "") %>%
    rename_with(str_replace, pattern = "_fo", replacement = "") %>%
    mutate(phenophase = "open flower") %>%
    relocate(phenophase, .after = "n_obs") 
  
  # Combine the two dataframes. Now each row contains onset data for a unique 
  # combination of plant, wateryr, and phenophase
  py <- rbind(py_fl, py_fo)
  
  # Remove plant-yrs when:
  # 1) flowers/open flowers weren't observed (firstyes = NA),
  # 2) "yes" to phenophase status wasn't preceded by a "no" (dayssince_lastno = NA)
  # 3) preceding "no" occurred more than XX days ago (dayssince_lastno > daysprior)
  py <- py %>%
    filter(!is.na(firstyes)) %>%
    filter(!is.na(dayssince_lastno)) %>%
    filter(dayssince_lastno <= daysprior)
  count(py, phenophase)
  # 2803, 3822 plant-yrs for flowers, open flowers, respectively
  
  # Summarize data for each species and phenophase, keeping only those that 
  # have at least the minimum number of plant-yrs (min_plantyrs). 
  py_spp <- py %>%
    mutate(loc = paste0(lon, ", ", lat)) %>%
    group_by(phenophase, common_name) %>%
    summarize(n_plantyrs = n(),
              n_yrs = length(unique(wateryr)),
              yr_min = min(wateryr),
              yr_max = max(wateryr),
              n_locs = length(unique(loc)),
              n_locs_west = length(unique(loc[lon < (-105)])),
              n_locs_central = length(unique(loc[lon >= (-105) & lon < (-90)])),
              n_locs_east = length(unique(loc[lon >= (-90)])),
              n_locs_LA = length(unique(loc[state == "LA"])),
              n_locs_NM = length(unique(loc[state == "NM"])),
              n_locs_OK = length(unique(loc[state == "OK"])),
              n_locs_TX = length(unique(loc[state == "TX"])),
              lat_min = floor(min(lat)),
              lat_max = ceiling(max(lat)),
              lon_min = floor(min(lon)),
              lon_max = ceiling(max(lon)),
              .groups = "keep") %>%
    data.frame() %>%
    filter(n_plantyrs >= min_plantyrs) %>%
    mutate(prop_west = n_locs_west / n_locs,
           prop_central = n_locs_central / n_locs,
           prop_east = n_locs_east / n_locs) %>%
    mutate(region = case_when(
      prop_west > 0.5 ~ "west",
      prop_central > 0.5 ~ "central",
      prop_east > 0.5 ~ "east",
      prop_west == 0.5 & prop_central == 0.5 ~ "west-central",
      prop_west == 0.5 & prop_east == 0.5 ~ "west-east",
      prop_east == 0.5 & prop_central == 0.5 ~ "east-central",
      .default = NA)) %>%
    rename(LA = n_locs_LA,
           NM = n_locs_NM,
           OK = n_locs_OK,
           TX = n_locs_TX) %>%
    mutate(years = paste0(yr_min, "-", yr_max)) %>%
    select(phenophase, common_name, n_plantyrs, years, region, n_locs, LA, NM, OK,
           TX, lat_min, lat_max, lon_min, lon_max) %>%
    arrange(phenophase, desc(n_plantyrs))
  
  # If no plants for a given species-phenophase are in the SC region, calculate
  # the distance (in km) from the nearest plant to the SC region. (Will need to 
  # reproject polygons & points to CONUS Albers equal-area [ESRI:102003])
  
  sc_region <- states_sc <- subset(states, states$sc == 1)
  sc_region <- aggregate(sc_region)
  sc_region <- terra::project(sc_region, "ESRI:102003")
  
  # Extract plant locations
  py_locs <- py %>%
    select(phenophase, common_name, lon, lat) %>%
    distinct() %>%
    vect(geom = c("lon", "lat"), crs = "epsg:4326")
  py_locs <- terra::project(py_locs, "ESRI:102003")
  
  py_spp$dist_to_SC_km <- NA
  for (i in 1:nrow(py_spp)) {
    if (rowSums(py_spp[i, c("LA", "NM", "OK", "TX")]) > 0) {next}
    spp_ph_locs <- subset(py_locs, 
                          py_locs$common_name == py_spp$common_name[i] & 
                            py_locs$phenophase == py_spp$phenophase[i])
    dists <- terra::distance(spp_ph_locs, sc_region, unit = "km")
    py_spp$dist_to_SC_km[i] <- round(min(dists))
  }
  
  # Save this summary of flowering onset data to file
  csv_file <- "output/flowering-onsets/flower-onset-data-summaries.csv"
  if (!file.exists(csv_file) | (file.exists(csv_file) & replace == TRUE)) {
    write.csv(py_spp,
              file = csv_file,
              row.names = FALSE)
  }
  
# Create a pdf, summarizing data for flower phenophase ------------------------#

  pheno_name <- "flower"
  fl_pdf_name <- "output/flowering-onsets/flower-onset-data-summaries.pdf"
  
  # Create a flextable:
  py_fl_pdf <- py_spp %>%
    filter(phenophase == pheno_name) %>%
    mutate(lat_range = paste0(lat_min, "-", lat_max)) %>%
    select(-c(phenophase, lat_min, lat_max, lon_min, 
              lon_max, dist_to_SC_km)) %>%
    arrange(desc(n_plantyrs))
  
  fl_cap_start <- paste("Summary of data available to model variation in", 
                        pheno_name)
  cap_end <- paste("onset date in species with at least", min_plantyrs, 
                   "plant-years of data, where a plant-year is all observations",
                   "made of a plant in a particular water year (2023 = 1 Oct",
                   "2022 - 30 Sep 2023). A plant-year was considered to have a", 
                   "valid onset date if the first 'yes' was preceded by a 'no'",
                   "within", daysprior, "days. Region indicates where the",
                   "majority of plants were located in the US (e.g., east =",
                   "eastern US). Table lists the number of unique plant",
                   "locations in the US, the number of unique plant locations",
                   "in each of the four states in the southcentral region, and",
                   "the latitudinal range (in degrees) of observed plant",
                   "locations.")
  fl_cap <- paste(fl_cap_start, cap_end)
  
  fl_flext <- flextable(py_fl_pdf) %>%
    add_header_row(
      top = TRUE,
      values = c(rep("", 4), "No. locations", rep("", 5))) %>%
    set_header_labels(
      common_name = "Species",
      n_plantyrs = "No. plant-yrs",
      years = "Year range",
      region = "Region",
      n_locs = "Total",
      LA = "LA",
      NM = "NM",
      OK = "OK",
      TX = "TX",
      lat_range = "Latitudinal range") %>%
    merge_at(i = 1, j = 5:9, part = "header") %>%
    flextable::align(align = "center", j = 3:10, part = "all") %>%
    flextable::width(j = 1, width = 1.8) %>%
    flextable::width(j = 2, width = 0.8) %>%
    flextable::width(j = 3, width = 1) %>%
    flextable::width(j = 4, width = 0.7) %>%
    flextable::width(j = 5, width = 0.5) %>%
    flextable::width(j = 6:9, width = 0.3) %>%
    flextable::width(j = 10, width = 0.8) %>%
    add_header_lines(values = fl_cap) %>%
    hline_top(border = fp_border_default(width = 0), part = "header") %>%
    hline(i = 2, border = fp_border_default(width = 0), part = "header") %>%
    hline(i = 2, j = 5:9, border = fp_border_default(), part = "header")
  # fl_flext
  
  # Get common names of species with sufficient data
  fl_common_names <- py_fl_pdf$common_name
  
  # Create map for each species with locations of plants that have one or more
  # years of flowering onset data. Put ggplot objects in list (fl_ggplots)
  fl_ggplots <- list()
  for (i in 1:length(fl_common_names)) {
    
    # Extract onset data
    py_fl_spp <- py %>%
      filter(common_name == fl_common_names[i]) %>%
      filter(phenophase == "flower")
    
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
  
  # Create first page of pdf (title and table)
  fl_title <- ggdraw() + 
    draw_label("Onset of flower phenophase", 
               fontface = "bold", x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  fl_title_pdf <- plot_grid(NULL, fl_title, NULL, nrow = 1, 
                            rel_widths = c(0.5, pdfw - 1, 0.5))
  
  fl_flext_aspect <- flextable_dim(fl_flext)$aspect_ratio
  fl_flext_h <- fl_flext_aspect * (pdfw - 1)
  fl_flext_grob <- gen_grob(fl_flext, fit = "auto", just = "center")
  fl_flext_pdf <- plot_grid(NULL, fl_flext_grob, NULL, 
                            nrow = 1, 
                            rel_widths = c(0.5, pdfw - 1, 0.5))
  page1_fl <- plot_grid(NULL, fl_title_pdf, fl_flext_pdf, NULL,
                        ncol = 1, 
                        rel_heights = c(0.1, 0.5, fl_flext_h, 
                                        pdfh - (fl_flext_h + 0.6)))
  page1_fl_filename <- "output/flowering-onsets/flower-onset-page1.pdf"
  
  if (!file.exists(fl_pdf_name) | (file.exists(fl_pdf_name) & replace == TRUE)) {
    ggsave(filename = page1_fl_filename,
           plot = page1_fl,
           width = pdfw,
           height = pdfh,
           units = "in",
           device = cairo_pdf)
  }
    
  # Pages to follow with species maps
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
  
  if (!file.exists(fl_pdf_name) | (file.exists(fl_pdf_name) & replace == TRUE)) {
    
    for (j in 1:length(pages)) {
      pdf_name <- paste0("output/flowering-onsets/flower-onset-map-", j, ".pdf")
      ggsave(filename = pdf_name, 
             plot = pages[[j]],
             width = pdfw, 
             height = pdfh, 
             units = "in", 
             device = cairo_pdf)
    }

    # Combine all pages and remove the individual pdfs
    ind_map_pages <- list.files(path = "output/flowering-onsets/",
                                pattern = "flower-onset-map-",
                                full.names = TRUE)
    invisible(pdf_combine(c(page1_fl_filename, ind_map_pages), 
                          output = fl_pdf_name))
    invisible(file.remove(c(page1_fl_filename, ind_map_pages)))
  }
  
# Create a pdf, summarizing data for open flower phenophase -------------------#
  
  pheno_name <- "open flower"
  fo_pdf_name <- "output/flowering-onsets/open-flower-onset-data-summaries.pdf"
  
  # Create a flextable:
  py_fo_pdf <- py_spp %>%
    filter(phenophase == pheno_name) %>%
    mutate(lat_range = paste0(lat_min, "-", lat_max)) %>%
    select(-c(phenophase, lat_min, lat_max, lon_min, 
              lon_max, dist_to_SC_km)) %>%
    arrange(desc(n_plantyrs))
  
  fo_cap_start <- paste("Summary of data available to model variation in", 
                        pheno_name)
  fo_cap <- paste(fo_cap_start, cap_end)
  
  fo_flext <- flextable(py_fo_pdf) %>%
    add_header_row(
      top = TRUE,
      values = c(rep("", 4), "No. locations", rep("", 5))) %>%
    set_header_labels(
      common_name = "Species",
      n_plantyrs = "No. plant-yrs",
      years = "Year range",
      region = "Region",
      n_locs = "Total",
      LA = "LA",
      NM = "NM",
      OK = "OK",
      TX = "TX",
      lat_range = "Latitudinal range") %>%
    merge_at(i = 1, j = 5:9, part = "header") %>%
    flextable::align(align = "center", j = 3:10, part = "all") %>%
    flextable::width(j = 1, width = 1.8) %>%
    flextable::width(j = 2, width = 0.8) %>%
    flextable::width(j = 3, width = 1) %>%
    flextable::width(j = 4, width = 0.7) %>%
    flextable::width(j = 5, width = 0.5) %>%
    flextable::width(j = 6:9, width = 0.3) %>%
    flextable::width(j = 10, width = 0.8) %>%
    add_header_lines(values = fo_cap) %>%
    hline_top(border = fp_border_default(width = 0), part = "header") %>%
    hline(i = 2, border = fp_border_default(width = 0), part = "header") %>%
    hline(i = 2, j = 5:9, border = fp_border_default(), part = "header")
  # fo_flext 

  # Get common names of species with sufficient data
  fo_common_names <- py_fo_pdf$common_name
  
  # Create map for each species with locations of plants that have one or more
  # years of open flower onset data. Put ggplot objects in list (fo_ggplots)
  fo_ggplots <- list()
  for (i in 1:length(fo_common_names)) {
    
    # Extract onset data
    py_fo_spp <- py %>%
      filter(common_name == fo_common_names[i]) %>%
      filter(phenophase == "open flower")
    
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

  # Create first page of pdf (title and table)  
  fo_title <- ggdraw() + 
    draw_label("Onset of open flower phenophase", 
               fontface = "bold", x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  fo_title_pdf <- plot_grid(NULL, fo_title, NULL, nrow = 1, 
                            rel_widths = c(0.5, pdfw - 1, 0.5))
  
  fo_flext_aspect <- flextable_dim(fo_flext)$aspect_ratio
  fo_flext_h <- fo_flext_aspect * (pdfw - 1)
  fo_flext_grob <- gen_grob(fo_flext, fit = "auto", just = "center")
  fo_flext_pdf <- plot_grid(NULL, fo_flext_grob, NULL, 
                            nrow = 1, 
                            rel_widths = c(0.5, pdfw - 1, 0.5))
  page1_fo <- plot_grid(NULL, fo_title_pdf, fo_flext_pdf, NULL,
                        ncol = 1, 
                        rel_heights = c(0.1, 0.5, fo_flext_h, 
                                        pdfh - (fo_flext_h + 0.6)))
  page1_fo_filename <- "output/flowering-onsets/open-flower-onset-page1.pdf"
  
  if (!file.exists(fo_pdf_name) | (file.exists(fo_pdf_name) & replace == TRUE)) {
    ggsave(filename = page1_fo_filename,
           plot = page1_fo,
           width = pdfw,
           height = pdfh,
           units = "in",
           device = cairo_pdf)
  }
  
  # Pages to follow with species maps
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
  
  if (!file.exists(fo_pdf_name) | (file.exists(fo_pdf_name) & replace == TRUE)) {
  
    for (j in 1:length(pages)) {
      pdf_name <- paste0("output/flowering-onsets/open-flower-onset-map-", 
                         j, ".pdf")
      ggsave(filename = pdf_name, 
             plot = pages[[j]],
             width = pdfw, 
             height = pdfh, 
             units = "in", 
             device = cairo_pdf)
    }
    
    # Combine all pages and remove the individual pdfs
    ind_map_pages <- list.files(path = "output/flowering-onsets/",
                                pattern = "open-flower-onset-map-",
                                full.names = TRUE)
    invisible(pdf_combine(c(page1_fo_filename, ind_map_pages), 
                          output = fo_pdf_name))
    invisible(file.remove(c(page1_fo_filename, ind_map_pages)))
  }
  