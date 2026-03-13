################################################################################
# For priority species: Evaluate quantity and quality of data available to 
# visualize weekly proportion of plants with open flowers

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
require(pdftools)
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

# Load shapefile with state boundaries ----------------------------------------# 
  
  states <- vect("data/states/cb_2017_us_state_500k.shp")
  states <- subset(states, 
                   !states$STUSPS %in% c("HI", "AK", "VI", "MP", 
                                         "GU", "PR", "AS"))
  states$sc <- as.factor(ifelse(states$STUSPS %in% c("LA", "NM", "OK", "TX"), 
                                1, 0))  
  
# Load processed NPN status/intensity data and format -------------------------#
  
  df <- read.csv("data/flower-status-intensities-priorityspp.csv")
  
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
    group_by(common_name, plant_id, site_id, state, zone, year) %>%
    summarize(n_obs = n(),                      # No. of daily observations
              n_status_fl = sum(!is.na(status_fl)),  # No. of flower status obs
              n_status_fo = sum(!is.na(status_fo)),  # No. of open flower status obs
              n_value_fl = sum(!is.na(midpoint_fl)), # No. of flower intensity values
              n_value_fo = sum(!is.na(midpoint_fo)), # No. of open flower intensity values
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
  
  # Remove any plant-year with no open flower status recorded
  samples <- samples %>% 
    filter(n_status_fo > 0) %>%
    mutate(plant_yr = paste(plant_id, year, sep = "_"))
  df <- df %>%
    mutate(plant_yr = paste(plant_id, year, sep = "_")) %>%
    filter(plant_yr %in% samples$plant_yr)
  
# Add week to data, so we can calculate weekly proportions --------------------#
  
  # We'll also create wk_doy columns to assign each week with a day of the year:
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
  
  # Keep just one observation of each plant, each week. Sort so the most
  # advanced phenophase gets kept (if more than one value in a week)
  df1 <- df %>%
    arrange(common_name, plant_id, year, wk, 
            desc(status_fl), desc(status_fo)) %>%
    distinct(plant_id, year, wk, .keep_all = TRUE) %>%
    # Add an indicator for plants in southcentral states
    mutate(sc = 1 * state %in% c("LA", "NM", "OK", "TX"))
  
# Summarize weekly data available, by species ---------------------------------#  
  
  # Summarize quantity of data available across all plants, years
  propdat_wk <- df1 %>%
    group_by(common_name, wk) %>%
    summarize(n_obs = n(), .groups = "keep") %>%
    data.frame()
  # Calculate mean no. observations used to calculate proportion open in weeks 
  # 10-40 (~ March - September)
  propdat <- propdat_wk %>%
    filter(wk %in% 10:40) %>%
    group_by(common_name) %>%
    summarize(mn_wkobs = round(mean(n_obs), 1)) %>%
    data.frame()
  
  # Summarize quantity of data available across plants in SC region, all years
  propdat_wk_sc <- df1 %>%
    filter(sc == 1) %>%
    group_by(common_name, wk) %>%
    summarize(n_obs_sc = n(), .groups = "keep") %>%
    data.frame()
  propdat_sc <- propdat_wk_sc %>%
    filter(wk %in% 10:40) %>%
    group_by(common_name) %>%
    summarize(mn_wkobs_sc = round(mean(n_obs_sc), 1)) %>%
    data.frame()
  
  # Summarize number of plants, plant-years
  nplants <- df1 %>%
    group_by(common_name) %>%
    summarize(n_plants = length(unique(plant_id)),
              n_plants_sc = length(unique(plant_id[sc == 1])),
              n_plantyrs = length(unique(plant_yr)),
              n_plantyrs_sc = length(unique(plant_yr[sc == 1]))) %>%
    data.frame()
  
  # Combine everything
  propdata <- spp %>%
    right_join(nplants, by = "common_name") %>%
    left_join(propdat, by = "common_name") %>%
    left_join(propdat_sc, by = "common_name") %>%
    arrange(priority, desc(n_plantyrs_sc))
  
  # Look at species summary
  propdata
  
# Identify priority species with sufficient open-flower status data -----------#
  
  # Set a very liberal threshold for mean number of plants in the SC region 
  # observed each week in weeks 10-40. In other words, we'll only move 
  # forward with species that average X or more plants in the SC region 
  # observed each week (in weeks 10-40) when combining data across all years.
  threshold <- 6
  
  species <- propdata %>%
    filter(mn_wkobs_sc >= threshold)
  
  # List of common names
  common_names <- species$common_name
  
  # List of common names with underscores for filenames
  nice_names <- str_replace(common_names, " ", "_")
  
  species <- species %>%
    mutate(LA1 = ifelse(LA == 1, "LA", NA), 
           NM1 = ifelse(NM == 1, "NM", NA),
           OK1 = ifelse(OK == 1, "OK", NA),
           TX1 = ifelse(TX == 1, "TX", NA)) %>%
    unite("states", LA1:TX1, sep = ", ", na.rm = TRUE)
  
  # List of species, with priority info, for titles on species summary sheets
  species_titles <- paste0(str_to_sentence(common_names),
                           " (priority in ", species$states, ")")

# Loop through species, creating summary table and figures --------------------#

for (i in 1:length(common_names)) {

  # Extract data for species
  sppdat <- df1 %>% filter(common_name == common_names[i])
  
  # Isolate location information
  sppsites <- sppdat %>%
    select(site_id, lat, lon, state, sc) %>%
    distinct() %>%
    mutate(sc = as.factor(sc))

  # Create map with locations of monitored plants
  sppsitesv <- vect(sppsites, geom = c("lon", "lat"), crs = "epsg:4269")
  spp_map <- ggplot(data = states) +
    geom_spatvector(aes(fill = sc), color = "gray60") +
    scale_fill_discrete(type = c("white", "gray95"), guide = "none") +
    geom_spatvector(data = sppsitesv, aes(color = sc), size = 1.5) +
    scale_color_manual(values = c(cols2[1], "red"), guide = "none") +
    labs(title = paste0("Locations of monitored plants inside (red) and ",
                        "outside (green) the Southcentral region (shaded gray)")) +
    theme(panel.background = element_rect(fill = rgb(t(col2rgb("lightblue")), 
                                                     alpha = 0.2 * 255, 
                                                     maxColorValue = 255)),
          plot.title = element_text(family = "Arial", size = 11))
  # spp_map

  # Summarize how much data is available, by year and across years
  # (focus on weeks 10-40)
  spp_nobs <- sppdat %>%
    group_by(year) %>%
    summarize(nplants = length(unique(plant_id)),
              nplants_sc = length(unique(plant_id[sc == 1]))) %>%
    data.frame()
  spp_nobsperwk <- sppdat %>%
    filter(wk %in% 10:40) %>%
    group_by(year, wk) %>%
    summarize(nobs_open = n(), 
              nobs_open_sc = length(year[sc == 1]),
              .groups = "keep") %>%
    group_by(year) %>%
    summarize(nobs_weeklymn = mean(nobs_open),
              nobs_weeklymn_sc = mean(nobs_open_sc)) %>%
    data.frame()
  spp_nobs <- spp_nobs %>%
    left_join(spp_nobsperwk, by = "year")
  
  nobs_open <- sppdat %>%
    filter(wk %in% 10:40) %>%
    group_by(wk) %>%
    summarize(nobs_open = n(),
              nobs_open_sc = length(wk[sc == 1])) %>%
    data.frame()
  
  overall <- data.frame(
    year = "all",
    nplants = length(unique(sppdat$plant_id)),
    nplants_sc = length(unique(sppdat$plant_id[sppdat$sc == 1])),
    nobs_weeklymn = mean(nobs_open$nobs_open),
    nobs_weeklymn_sc = mean(nobs_open$nobs_open_sc))
  
  nobs_tab <- rbind(spp_nobs, overall) %>%
    mutate(across(contains("weekly"),  \(x) round(x, 1))) %>%
    relocate(nobs_weeklymn, .after = nplants) %>%
    mutate(year = ifelse(year == "all", "All years", year))
  
  # Create a warning when the mean number of plants observed each week in weeks 
  # 10-40 is less than the minimum number listed below
  weekly_nobs_min <- 8
  
  # Print message about data available for this species
  sppyears <- nobs_tab %>% 
    filter(nobs_weeklymn_sc >= weekly_nobs_min)
  if (nrow(sppyears) == 0) {
    sample_size_warning <- paste0("Data are very limited for ", common_names[i], 
                                  ". Mean number of observations per week ",
                                  "(combining data across years) < ", 
                                  weekly_nobs_min, ".")
    warning(sample_size_warning)
  } else{
    if (nrow(sppyears) == 1) {
      message("Can only characterize open flower phenophase for ", 
              common_names[i], " by combining data across all years.")
    } else {
      message("Data are sufficient to characterize open flower phenophase for ",
              common_names[i], " in ", paste(sppyears$year[-nrow(sppyears)], 
              collapse = ", "))
    }
  }

  # Put summary stats in a flextable so it can be combined with ggplot objects
  caption <- paste0("Number of plants and the mean number of open flower ",
                    "observations made each week (in weeks 10-40) for all ",
                    "plants in the Continental US and for plants in the ",
                    "Southcentral region.")
  if (exists("sample_size_warning")) {
    caption <- as_paragraph(
      as_chunk(caption), "\n", 
      as_chunk(paste0("Note: ", sample_size_warning), 
               fp_text_default(bold = TRUE, color = "red"))
    )
    rm(sample_size_warning)
  }
  
  flext <- flextable(nobs_tab) %>%
    add_header_row(
      top = TRUE,
      values = c("", 
                 "Continental US",
                 "",
                 "Southcentral region",
                 "")) %>%
    set_header_labels(
      year = "Year",
      nplants = "No. plants", 
      nobs_weeklymn = "Mean no. obs/wk",
      nplants_sc = "No. plants", 
      nobs_weeklymn_sc = "Mean no. obs/wk") %>%
    merge_at(i = 1, j = 2:3, part = "header") %>%
    merge_at(i = 1, j = 4:5, part = "header") %>%
    flextable::width(j = c(1, 2, 4), width = 1.2) %>%
    flextable::width(j = c(3, 5), width = 1.5) %>%
    flextable::align(align = "center", j = 2:5, part = "all") %>%
    hline(i = nrow(nobs_tab) - 1) %>%
    add_header_lines(values = caption) %>%
    hline_top(border = fp_border_default(width = 0), part = "header")
  # flext

  # Format data for figures with weekly proportion of open flowers, all years 
  # combined
  prop_allyrs <- sppdat %>%
    group_by(sc, wk, wk_date1, wk_doy1, wk_date4, wk_doy4) %>%
    summarize(n_obs = n(),
              n_open = sum(status_fo),
              .groups = "keep") %>%
    data.frame() %>%
    mutate(prop_open = n_open / n_obs) %>%
    mutate(sc_f = as.factor(ifelse(sc == 0, "Other", "SC")))
  
  # To create figures with ticks between month labels on the x axis, need to do 
  # a little extra work. 
    # Dates where we want month labels (15th of month)
    x_lab <- as.Date(paste0("2024-", 1:12, "-15"))
    # Dates where we want ticks (1st of month)
    x_tick <- as.Date(c(paste0("2024-", 1:12, "-01"), "2025-01-01"))
    n_x_tick <- length(x_tick)
    # Will specify axis breaks & ticks at 1st and 15th of month. Make labels on
    # the 1st black and change color of tick marks on the 15th to NA.

  # Create figure comparing weekly proportion of plants with open flowers inside
  # and outside of the SC region (to see whether it's reasonable to use data 
  # from plants across the US in state info sheets)
  prop_plot <- ggplot(data = prop_allyrs,
                      aes(x = wk_date4, y = prop_open, group = sc_f)) + 
    geom_point(aes(color = sc_f, size = n_obs), alpha = alphaline) +
    geom_line(aes(color = sc_f), alpha = alphaline) +
    scale_color_whitebox_d(palette = whitebox_palette) +
    scale_x_continuous(breaks = c(x_lab, x_tick),
                       labels = c(month.abb, rep("", n_x_tick))) +
    labs(y = paste0("Proportion of plants with open flowers")) +
    labs(color = "Region", size = "No. obs") +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8),
          axis.ticks.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                rep("black", n_x_tick))),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                      rep("white", n_x_tick))),
          panel.background = element_rect(fill = "gray95"),
          axis.title.x = element_blank())
  # prop_plot  

  # Extracting proportion based on plants in SC region only (combining
  # data across years)
  propSC_allyrs <- prop_allyrs %>%
    filter(sc == 1)
  
  # Use a GAM to model weekly proportions. (Note that we're using an iterative
  # process to select the smoothing parameter [k]. Seems that if k is set
  # pretty high, then sometimes you can get a very wide credible interval around
  # predictions during a portion of hte year when the curve is flat [ie, no 
  # plants flowering]. Seems reasonable to start by setting k somewhat low and
  # then checking to see if model fit is adequate. If not, can increase k)
  k_values <- seq(5, 50, by = 5)
  
  for (k_index in 1:length(k_values)) {
    gam1 <- gam(prop_open ~ s(wk_doy4, bs = "cc", k = k_values[k_index]), weights = n_obs,
                data = propSC_allyrs, method = "REML", family = "binomial")
    k_p <- k.check(gam1)[,"p-value"] 
    if(k_p > 0.1) break()
  }
  k <- k_values[k_index]

  k_message1 <- paste0("Smoothing parameter (k) = ", k, ".") 
  if (k_p <= 0.1 & k_index == length(k_values)) {
    k_message2 <- "Despite a large k value, still indications of poor model fit."
    message(k_message1, " ", k_message2)
  } else {
    message(k_message1)
  }

  # Make GAM predictions
  gam1_preds <- data.frame(
    wk_doy4 = min(propSC_allyrs$wk_doy4):max(propSC_allyrs$wk_doy4),
    wk_date4 = as.Date(min(propSC_allyrs$wk_date4):max(propSC_allyrs$wk_date4)))
  gam1_preds <- cbind(gam1_preds,
                      as.data.frame(predict(gam1,
                                            newdata = gam1_preds,
                                            type = "link",
                                            se.fit = TRUE))) %>%
    rename(fitl = fit) %>%
    mutate(lcll = fitl - 1.96 * se.fit,
           ucll = fitl + 1.96 * se.fit,
           fit = exp(fitl) / (1 + exp(fitl)),
           lcl = exp(lcll) / (1 + exp(lcll)),
           ucl = exp(ucll) / (1 + exp(ucll))) %>%
    select(-c(lcll, ucll))

  # Create figure with GAM predictions on top of raw weekly proportions
  gam1_plot <- ggplot(data = propSC_allyrs, 
                      aes(x = wk_date4, y = prop_open)) + 
    geom_point(aes(size = n_obs), alpha = alphaline, color = cols2[2]) +
    geom_line(alpha = alphaline, color = cols2[2]) +
    geom_ribbon(data = gam1_preds, 
                aes(x = wk_date4, y = fit, ymin = lcl, ymax = ucl), 
                color = "gray", alpha = alphapoly, linetype = 0) +
    geom_line(data = gam1_preds, aes(x = wk_date4, y = fit)) +
    scale_x_continuous(breaks = c(x_lab, x_tick),
                       labels = c(month.abb, rep("", n_x_tick))) +
    labs(y = paste0("Proportion of plants with open flowers"), 
         size = "No. obs", 
         title = k_message1) +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 8),
          axis.ticks.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                rep("black", n_x_tick))),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                      rep("white", n_x_tick))),
          panel.background = element_rect(fill = "gray95"),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 8))
  if (exists("k_message2")) {
    gam1_plot <- gam1_plot +
      labs(subtitle = k_message2) +
      theme(plot.subtitle = element_text(size = 8, face = "bold", color = "red"))
    rm(k_message2)
  }
  # gam1_plot
  
  # Create a "heat map" with raw weekly proportions. Adding a 2nd panel below
  # with weekly sample sizes (and the preset threshold identified with a dashed
  # line)
  propSC_allyrs <- propSC_allyrs %>%
    mutate(tile_ht = 0.5)

  g1 <- ggplot(propSC_allyrs) +
    geom_tile(aes(x = wk_date4, y = tile_ht, fill = prop_open)) +
    scale_fill_gradient(low = "gray95", high = cols2[2], 
                        na.value = "white",
                        name = "Proportion \nopen", 
                        guide = guide_colorbar(position = "inside")) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(breaks = c(x_lab, x_tick),
                       labels = c(month.abb, rep("", n_x_tick)),
                       expand = c(0.01, 0.01)) +
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y.right = element_line(color = "black"),
          axis.ticks.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                rep("black", n_x_tick))),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 10),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.position.inside = c(-0.05, 0.5),
          plot.margin = unit(c(5.5, 5.5, 5.5, 30), "pt"),
          legend.margin = margin(c(5, 0, 5, 0)))
  g2 <- ggplot(propSC_allyrs) +
    geom_col(aes(x = wk_date4, y = n_obs), width = 5, fill = "gray70") +
    geom_hline(yintercept = weekly_nobs_min, 
               linetype = "dashed", color = "black") +
    labs(y = "No. observations") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(breaks = c(x_lab, x_tick),
                       labels = c(month.abb, rep("", n_x_tick)),
                       expand = c(0.01, 0.01)) +
    theme(axis.line = element_line(color = "black"),
          axis.ticks.x = element_line(color = c(rep(NA, n_x_tick - 1), 
                                                rep("black", n_x_tick))),
          axis.title.x = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 10),
          plot.margin = unit(c(5.5, 5.5, 5.5, 30), "pt"))
  heatmap2 <- plot_grid(g1, g2, nrow = 2, rel_heights = c(2, 1), align = "v")
  # heatmap2

  # Create ggplot object with title for pdf summary pages
  title <- ggdraw() + 
    draw_label(species_titles[i], fontface = "bold", x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  title_pdf <- plot_grid(NULL, title, NULL, nrow = 1, 
                         rel_widths = c(0.5, 7.5, 0.5))
  
  # Combine table and map into a single pdf
  flext_pdf <- gen_grob(flext, fit = "fixed", just = "center")
  map_pdf <- plot_grid(NULL, spp_map, NULL, nrow = 1, 
                       rel_widths = c(0.5, 7.5, 0.5))
  titletablemap <- plot_grid(NULL, title_pdf, flext_pdf, map_pdf, NULL,
                             ncol = 1, 
                             rel_heights = c(0.1, 0.5, 4.5, 4.5, 0.5))
  pdf1_name <- paste0("output/weekly-open-flower-proportions/", 
                     nice_names[i], 
                     "-tablemap.pdf")
  ggsave(pdf1_name, 
         width = pdfw, 
         height = pdfh, 
         units = "in", 
         device = cairo_pdf)

  # Combine figures with weekly proportions into a single pdf
  prop_pdf <- plot_grid(NULL, prop_plot, NULL, nrow = 1, 
                        rel_widths = c(0.5, 7.5, 0.5))
  gam1_pdf <- plot_grid(NULL, gam1_plot, NULL, nrow = 1, 
                        rel_widths = c(0.5, 7.5, 0.5))
  heatmap2_pdf <- plot_grid(NULL, heatmap2, NULL, nrow = 1, 
                            rel_widths = c(0.5, 7.5, 0.5))
  all_figs <- plot_grid(prop_pdf, NULL, gam1_pdf, NULL, heatmap2_pdf,
                        ncol = 1, rel_heights = c(2.8, 0.2, 3, 0.2, 2.8), align = "v")
  titlefigures <- plot_grid(NULL, title_pdf, all_figs, NULL,
                             ncol = 1, 
                             rel_heights = c(0.1, 0.5, 9, 0.5))
  pdf2_name <- paste0("output/weekly-open-flower-proportions/", 
                      nice_names[i], 
                      "-figs.pdf")
  ggsave(pdf2_name, 
         width = pdfw, 
         height = pdfh, 
         units = "in", 
         device = cairo_pdf)

  # Combine two pages into a single pdf and remove the 1-pagers
  pdf_name <- paste0("output/weekly-open-flower-proportions/", 
                     nice_names[i], ".pdf")
  invisible(pdf_combine(c(pdf1_name, pdf2_name), output = pdf_name))
  file.remove(pdf1_name)
  file.remove(pdf2_name)
  message("PDF with summary of open flower data (weekly proportions) for ", 
          common_names[i], " saved to output folder")
}

# If an extraneous pdf was saved to the repo, delete it now:
if (file.exists("Rplot001.pdf")) {
  invisible(file.remove("Rplot001.pdf"))
}
