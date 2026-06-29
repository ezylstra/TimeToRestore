library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(shiny)
library(bslib)
library(sf)
library(ggplot2)
library(leaflet)
library(RColorBrewer)
library(sortable)
library(grid)
library(cowplot)

# Load/format observation data ------------------------------------------------#

df <- read.csv("data/status-intensity-flowers-May2026.csv")
df <- df %>%
  mutate(obsdate = ymd(obsdate)) %>%
  mutate(region = case_when(
    region == "Dallas" ~ "Dallas/Ft. Worth",
    region == "Austin" ~ "Austin/San Antonio",
    .default = region
  ))
# For now limit to 4 state area and 2025-2026
df <- df %>%
  filter(state4 == 1) %>%
  filter(yr >= 2025)

# If there's more than one observation of a plant on a given date, select one 
# with positive status and (higher) intensity value
df <- df %>%
  arrange(id, obsdate, php_id, desc(status), desc(midpoint)) %>%
  distinct(id, obsdate, php_id, .keep_all = TRUE)

# Put all data for a plant, date in the same row (wide form)
dfw <- df %>%
  select(-c(observation_id, php)) %>%
  pivot_wider(names_from = c(php_id),
              values_from = c(status, intensity_value, midpoint)) %>%
  data.frame()

# If flower = 0, then open flower must be 0. Remove observations with open
# flower = 1 and change any open flower = NA to 0.
dfw <- dfw %>%
  filter(!(!is.na(status_500) & !is.na(status_501) & status_500 == 0 & status_501 == 1)) %>%
  mutate(status_501 = ifelse(!is.na(status_500) & status_500 == 0, 0, status_501))

# If open flower = 1, then flower must be 1. Change any that are NA.
dfw <- dfw %>%
  mutate(status_500 = ifelse(!is.na(status_501) & status_501 == 1, 1, status_500))

# Finally, change any midpoint values to 0 when status is 0
dfw <- dfw %>%
  mutate(midpoint_500 = ifelse(!is.na(status_500) & status_500 == 0, 
                               0, midpoint_500)) %>%
  mutate(midpoint_501 = ifelse(!is.na(status_501) & status_501 == 0,
                               0, midpoint_501))

# Create table with location options
loc_options <- distinct(dfw, state, region) %>%
  mutate(region = ifelse(is.na(region), "Statewide", region)) %>%
  rbind(data.frame(state = "SC states (4)", 
                   region = "Entire SC region"))

# ui --------------------------------------------------------------------------#

ui <- page_navbar(
  title = "Time to Restore",
  sidebar = sidebar(
    width = "25%",
    div(style = "font-size:90%",
        navset_tab(
          nav_panel(
            title = "Info",
            br(),
            HTML("This app uses data submitted since January 2025 to <i>USA-NPN 
                 Nature's Notebook</i> for species identified as priorities as 
                 part of the <i>Time to Restore</i> project. Under the 
                 <i><b>Data summaries</i></b> tab, plots display the weekly or
                 biweekly proportion of plants with flowers or open flowers 
                 within the selected area (4-state region, state, or for Texas
                 only, metropolitan area). Under the <i><b>Map</i></b> tab, a 
                 map displays the locations of all oberved plants for the 
                 selected species and region. Users can click on a point to get  
                 information about the site and number of plants observed."),
            br(),
            br(),
            HTML("<b>Species</b>: Users can select one or more species to  
                 display, in the order of their choice. Available species for 
                 the selected region include only those with observations from 
                 10 or more plants, each with 5 or more observations per 
                 year."),
            br(), 
            br(),
            HTML("<b>Visualizations</b>: Users can select one of three options  
                 for displaying weekly/biweekly proportions of plants with  
                 flowers or open flowers. Each visualization type provides 
                 indications when sample sizes are problematic (e.g., proportion 
                 based on <5 observations). A <i><b>bar chart</b></i> is the 
                 easiest way to visualize variation in effort over time. A
                 <i><b>bubble plot</b></i> is most effective for visualizing
                 flowering 'peaks'. A <i><b>heat map</b></i> can be effective 
                 for comparing phenology among species, but users should use  
                 with caution given that it is difficult to identify when 
                 proportions are based on few observations."),
            br(),
            br(),
            HTML("<b>Funding</b>: The <i>Time to Restore</i> project is 
                 supported by the South Central Climate Adapatation Science  
                 Center, which is managed by the U.S. Geological Survey.")
          ),
          nav_panel(
            title = "Settings",
            br(),
            selectInput(inputId = "state",
                        label = "State",
                        choices = c(unique(loc_options$state), ""),
                        selected = ""),
            uiOutput("locChoices"),
            radioButtons(inputId = "php", 
                         label = "Phenophase",
                         choices = c("Flowers" = "flower", 
                                     "Open flowers" = "open"),
                         inline = TRUE),
            selectInput(inputId = "years",
                        label = "Years",
                        choices = c("2025-2026 (combined)", "2025", "2026"),
                        selected = "2025-2026 (combined)"),
            radioButtons(inputId = "period",
                         label = "Summarization period",
                         choices = c("Weekly" = "weekly",
                                     "Biweekly" = "biweekly"),
                         inline = TRUE),
            uiOutput("speciesChoices"),
            radioButtons(inputId = "vistype",
                         label = "Visualization type",
                         choices = c("Bar chart",
                                     "Bubble plot",
                                     "Heat map"))
          )
        )
    )
  ),
  fillable = "Map",
  nav_panel(title = "Data summaries",
            plotOutput(outputId = "plot")),
  nav_panel(title = "Map", 
            leafletOutput(outputId = "map"),
            em("Locations in map have been adjusted slightly to limit overlap"))
)

# server ----------------------------------------------------------------------#

server <- function(input, output, session) { # add session for observeEvent reset

  loc_table <- reactive({
    req(input$state != "")
    loc_options %>% filter(state == input$state)
  })

  output$locChoices <- renderUI({
    selectInput(inputId = "loc", 
                label = "Location",
                choices = c(unique(loc_table()$region), ""),
                selected = "")
  })

  # Filter by location
  df_filteredloc <- reactive({
    req(input$loc, input$loc != "")
    if (input$loc == "Entire SC region") {
      dfw 
    } else if (input$loc == "Statewide") {
      dfw %>% filter(state == input$state)
    } else {
      dfw %>% filter(region == input$loc)
    }
  })

  # Filter by year
  df_filteredyr <- reactive({
    if (input$years == "2025-2026 (combined)") {
      df_filteredloc()
    } else {
      df_filteredloc() %>%
        filter(yr == as.numeric(input$years))
    }
  })

  # Filter by phenophase. Rename the relevant status/midpoint columns to
  # generic names so all downstream code is phenophase-agnostic
  df_filteredphp <- reactive({
    if (input$php == "flower") {
      df_filteredyr() %>%
        rename(status = status_500,
               midpoint = midpoint_500)
    } else {
      df_filteredyr() %>%
        rename(status = status_501,
               midpoint = midpoint_501)
    }
  })

  # Find species that have sufficient data (species with >= 10 plants, each with
  # >= 5 observations over years of interest)
  speciessubset <- reactive({
    df_filteredphp() %>%
      group_by(id) %>%
      mutate(nobs = n()) %>%
      ungroup() %>%
      group_by(spp) %>%
      summarize(nplants = n_distinct(id),
                nplants5 = n_distinct(id[nobs >= 5])) %>%
      data.frame() %>%
      filter(nplants5 >= 10)
  })

  # Render a drag-and-drop rank_list instead of a plain selectInput.
  # Users drag items from the "Available" bucket into "Selected & ordered" to
  # both choose species and set their facet order in one widget.
  # TODO: figure out how to left justify this section ####
  output$speciesChoices <- renderUI({
    tagList(
      tags$style(HTML("
        /* Compress height of species names */
        .default-sortable .rank-list-item {
          padding: 3px 8px !important;
          margin-bottom: 2px !important;
        }
        /* Enable scrolling in Available bucket*/
        #rank-list-species_available {
          max-height: 200px !important;
          overflow-y: auto !important;
        }
      ")),
      bucket_list(
        header = "Species (drag to select & order)",
        group_name = "species_bucket",
        orientation = "vertical",
        add_rank_list(
          text = "Available",
          labels = unique(speciessubset()$spp),
          input_id = "species_available"
        ),
        add_rank_list(
          text = "Selected & ordered",
          labels = NULL,
          input_id = "species_order"
        )
      )
    )
  })

  # Reset species_order bucket whenever state, location or year changes
  # so stale selections from the previous filter context don't carry over.
  observeEvent(list(input$state, input$loc, input$years), {
    updateSelectInput(session, "species_order", selected = character(0))
  }, ignoreInit = TRUE)
  
  # Filter to selected species. Guard against NULL/empty selection,
  # and against stale selections left over from a previous state/location
  # (e.g. right after switching `state`, input$species may still hold
  # species names that don't exist in the newly-filtered data until the
  # speciesChoices bucket_list above gets rebuilt).
  df_filteredspp <- reactive({
    req(input$species_order)
    valid_species <- intersect(input$species_order, speciessubset()$spp)
    req(length(valid_species) > 0)
    df_filteredphp() %>%
      filter(spp %in% valid_species) %>%
      mutate(spp = factor(spp, levels = valid_species)) %>% 
      mutate(wk = week(obsdate)) %>%
      # Remove observations in week 53 (Dec 31 [and Dec 30 in leap years])
      filter(wk < 53) %>%
      mutate(wk2 = ifelse(wk%%2 == 0, wk - 1, wk))
  })

  # Collapse site/species to one row per site-species, jitter coords once
  map_frame <- reactive({
    obslocs <- df_filteredspp() %>%
      group_by(site, site_name, lat, lon, spp) %>%
      summarize(plants = n_distinct(id), .groups = "drop") %>%
      data.frame()

    obslocs %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      st_jitter(factor = 0.0001) %>%
      cbind(st_coordinates(.)) %>%   # extract jittered coords from THIS object
      rename(lon = X, lat = Y) %>%
      st_drop_geometry()
  })

  # Derive palette from the actual species present in map_frame
  spp_colors <- reactive({
    spp_list <- unique(map_frame()$spp)
    n <- length(spp_list)
    req(n >= 1) # nothing to draw if no species have sufficient data
    # brewer.pal needs n >= 3; pad/trim as needed
    pal <- brewer.pal(max(3, n), "Paired")[seq_len(n)]
    setNames(pal, spp_list) # named vector: spp -> color
  })

  # Create map
  output$map <- renderLeaflet({
    mf    <- map_frame()
    cols  <- spp_colors()
    spp_list <- names(cols)

    m <- leaflet(data = mf) %>%
      addTiles(options = tileOptions(opacity = 0.6))

    # seq_along() avoids the 1:0 trap when spp_list is empty
    for (i in seq_along(spp_list)) {
      m <- m %>% addCircleMarkers(
        lng = ~lon,
        lat = ~lat,
        data = filter(mf, spp == spp_list[i]),
        group = spp_list[i],
        radius = 5,
        color = "black",
        fillColor = cols[[i]],
        fillOpacity = 0.9,
        stroke = TRUE,
        weight = 2,
        popup = ~paste0(spp, "<br>",
                        "Site ID: ", site, "<br>",
                        "No. plants: ", plants))
    }

    m %>%
      addLayersControl(overlayGroups = spp_list,
                       options = layersControlOptions(collapse = FALSE)) %>%
      addLegend(position  = "bottomright",
                colors    = unname(cols),
                labels    = spp_list,
                opacity   = 0.9)
  })

  # Calculate proportions after keeping just one observation of each plant, each
  # week or 2-week period. Sort so the most advanced phenophase gets kept (if 
  # more than one value in period)
  props <- reactive({
    if (input$period == "weekly") {
      props_temp <- df_filteredspp() %>%
        mutate(period = wk) %>%
        mutate(wk_doy4 = (period * 7) - 3) %>%
        mutate(wk_date4 = parse_date_time(x = paste(2025, wk_doy4), 
                                          orders = "yj"))
    } else {
      props_temp <- df_filteredspp() %>%
        mutate(period = wk2) %>%
        mutate(wk_doy4 = period * 7) %>%
        mutate(wk_date4 = parse_date_time(x = paste(2025, wk_doy4), 
                                          orders = "yj"))
    }
    
    # Keeping year in so if we're combining info across years, we retain an 
    # observation of a plant in each year
    props_temp %>%
      filter(!is.na(status)) %>%
      arrange(id, yr, period, desc(status), desc(midpoint)) %>% 
      distinct(id, yr, period, .keep_all = TRUE) %>%
      group_by(spp, period, wk_date4) %>%
      summarize(nobs = n(),
                nyes = sum(status),
                prop = nyes/nobs,
                .groups = "drop") %>%
      data.frame() %>% 
      mutate(obs_group = ifelse(nobs < 5, "low", "sufficient"))
  })

  # Plot aesthetics
  text_size <- 14
  yaxis_bubble <- reactive({
    if (input$php == "flower") {
      "Proportion with flowers"
    } else {
      "Proportion with open flowers"
    }
  })
  fill_bar <- reactive({
    if (input$php == "flower") {
      "Flowers"
    } else {
      "Open flowers"
    }
  })
  bar_width <- reactive({
    if (input$period == "weekly") {3} else {7}
  })
  heatmap_width <- reactive({
    if (input$period == "weekly") {6} else {12}
  })

  # Create function for size breaks that include low values (for bubble plot)
  pretty_breaks <- function(x, n = 3) {
    pr <- pretty(x, n)
    pr <- pr[pr > 0]
    pr <- sort(unique(c(1, 4, pr)))
    return(pr)
  }
  
  # Bar chart
  output$plot <- renderPlot(
    height = function() max(300, 250 * length(input$species_order)),
    {
      req(input$species_order)
      p_data <- props()
      
      if (input$vistype == "Bar chart") {
        p_bar <- p_data %>%
          ggplot(aes(x = wk_date4, y = nobs)) +
          geom_bar(aes(alpha = obs_group, fill = "No", color = "No"),
                   stat = "identity", width = bar_width()) +
          geom_bar(aes(y = nyes, alpha = obs_group, fill = "Yes", color = "Yes"),
                   stat = "identity", width = bar_width()) +
          geom_hline(yintercept = 5, linetype = "dotted", color = "salmon3", linewidth = 0.8) +
          facet_wrap(~ spp, ncol = 1, scales = "free_x") +
          scale_alpha_manual(values = c("low" = 0.5, "sufficient" = 1),
                             guide = "none") +   # suppress separate alpha legend
          scale_fill_manual(values = c("No" = "gray", "Yes" = "steelblue3")) +
          scale_color_manual(values = c("No" = "gray", "Yes" = "steelblue3")) +
          scale_x_date(limits = c(as.Date("2025-01-01"), as.Date("2025-12-31")),
                       expand = 0.02,
                       date_breaks = "1 month",
                       date_labels = "%e %b") +
          labs(x = "Date", y = "No. observations", 
               fill = fill_bar(), color = fill_bar()) +
          theme_bw() +
          theme(legend.position = "top",
                axis.text = element_text(size = text_size),
                axis.title = element_text(size = text_size),
                legend.text = element_text(size = text_size),
                legend.title = element_text(size = text_size),
                strip.text = element_text(size = text_size + 1),
                panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "#b6d3b6"))
        
        # Pull the legend out of the plot as its own grob so it can be placed
        # next to the title instead of stacked above the plot
        legend_grob <- cowplot::get_legend(p_bar)
        
        title_text <- str_wrap(
          paste0("Red dotted line indicates minimum desired sample size to ",
                 "evaluate proportions. Periods with fewer samples have ",
                 "semi-transparent bars."),
          width = 80
        )
        title_grob <- grid::textGrob(
          title_text,
          x = 0, hjust = 0,
          gp = grid::gpar(fontsize = text_size)
        )
        
        # Anchor the legend to the right edge of its cell (instead of the
        # default center-justify, which can clip a wide legend at the
        # plot's right edge)
        legend_panel <- cowplot::ggdraw() +
          cowplot::draw_grob(legend_grob, x = 0.98, hjust = 1)
        
        # Title on the left, legend on the right, in one row
        top_row <- cowplot::plot_grid(title_grob, legend_grob,
                                      nrow = 1, rel_widths = c(4, 1))
        
        # That row stacked above the plot (with its own legend removed)
        total_height <- max(300, 250 * length(input$species_order))
        cowplot::plot_grid(top_row,
                           p_bar + theme(legend.position = "none"),
                           ncol = 1, 
                           rel_heights = c(50/total_height, 
                                           (total_height - 50)/total_height))

      } else if (input$vistype == "Bubble plot") {

        # Compute breaks the same way ggplot will: apply pretty_breaks to the
        # data range, then drop anything outside it (mimicking ggplot's behavior)
        data_range <- range(p_data$nobs, na.rm = TRUE)
        br <- pretty_breaks(p_data$nobs)
        br <- br[br >= data_range[1] & br <= data_range[2]]  # keep only in-range breaks
        min_sample_size <- 5
        
        aes_fill  <- ifelse(br < min_sample_size, "white",   "steelblue3")

        p_bubble <- p_data %>%
          ggplot(aes(x = wk_date4, y = prop)) +
          geom_line(color = "gray50") +
          geom_point(aes(size = nobs, fill = obs_group), color = "steelblue3",
                     shape = 21) +
          scale_fill_manual(values = c("low" = "white", "sufficient" = "steelblue3"),
                            guide = "none") +    # suppress separate fill legend
          facet_wrap(~ spp, ncol = 1, scales = "free_x") +
          scale_size_continuous(range = c(2, 10),
                                breaks = br,
                                guide = guide_legend(
                                  title = "No. observations",
                                  nrow = 1,
                                  override.aes = list(fill = aes_fill))
          ) +
          scale_x_date(limits = c(as.Date("2025-01-01"), as.Date("2025-12-31")),
                       expand = 0.02,
                       date_breaks = "1 month",
                       date_labels = "%e %b") +
          scale_y_continuous(limits = c(-0.05, 1.05)) +
          labs(x = "Date", y = yaxis_bubble()) +
          theme(legend.position = "top",
                axis.text = element_text(size = text_size),
                axis.title = element_text(size = text_size),
                legend.text = element_text(size = text_size),
                legend.title = element_text(size = text_size),
                strip.text = element_text(size = text_size + 1),
                panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "#b6d3b6"))
        
        # Pull the legend out of the plot as its own grob so it can be placed
        # next to the title instead of stacked above the plot
        legend_grob <- cowplot::get_legend(p_bubble)
        
        title_text <- str_wrap(
          paste0("Size of bubble varies with number of observations. ",
                 "White-filled bubbles indicate that the proportion is based on ",
                 "fewer than the minimum desired number of observations (",
                 min_sample_size, ")."),
          width = 100
        )
        title_grob <- grid::textGrob(
          title_text,
          x = 0, hjust = 0,
          gp = grid::gpar(fontsize = text_size)
        )
        
        # Anchor the legend to the right edge of its cell (instead of the
        # default center-justify, which can clip a wide legend at the
        # plot's right edge)
        legend_panel <- cowplot::ggdraw() +
          cowplot::draw_grob(legend_grob, x = 0.95, hjust = 1)
        
        # Title on the left, legend on the right, in one row
        top_row <- cowplot::plot_grid(title_grob, legend_grob, 
                                      nrow = 1, rel_widths = c(2, 1))
        
        # That row stacked above the plot (with its own legend removed)
        total_height <- max(300, 250 * length(input$species_order))
        cowplot::plot_grid(top_row,
                           p_bubble + theme(legend.position = "none"),
                           ncol = 1, 
                           rel_heights = c(50/total_height, 
                                           (total_height - 50)/total_height))
      
      } else {

        min_sample_size <- 5
        
        # Add a dummy column so we can inject an NA key into the legend
        p_data$na_label <- factor(ifelse(is.na(p_data$prop), 
                                         "No observations", NA),
                                  levels = "No observations")
        
        p_heat <- p_data %>%
          ggplot() +
          geom_bar(aes(x = wk_date4, y = 0.9, fill = prop), stat = "identity",
                   width = heatmap_width()) +
          # Invisible geom that exists only to create the legend key
          geom_tile(aes(x = wk_date4, y = -99, color = na_label),
                    fill = "white", linewidth = 0.5, width = heatmap_width(), 
                    height = 0) +
          scale_fill_gradient(low = "gray95", high = "steelblue3",
                              na.value = "white",
                              name = paste(yaxis_bubble(), "    "),
                              guide = guide_colorbar(order = 1)) +
          scale_color_manual(
            values = c("No observations" = "black"),
            labels = "",
            name = "No data",
            guide = guide_legend(
              order = 2,
              override.aes = list(
                fill  = "white",
                color = "black",
                size  = 6,
                shape = 22,
                linewidth = 0.5
              )
            )
          ) +
          geom_text(data = filter(p_data, obs_group == "low"),
                    aes(x = wk_date4, y = 0.95, label = nobs), 
                    color = "red", fontface = 2, size = 5) +
          geom_text(data = filter(p_data, obs_group == "sufficient"),
                    aes(x = wk_date4, y = 0.95, label = nobs), size = 4) +
          scale_y_continuous(expand = c(0, 0)) +
          coord_cartesian(ylim = c(0, 1)) +
          scale_x_date(limits = c(as.Date("2025-01-01"), as.Date("2026-01-01")),
                       expand = 0.02,
                       date_breaks = "1 month",
                       date_labels = "%e %b") +
          facet_wrap(~ spp, ncol = 1, scales = "free_x") +
          labs(x = "Date") +
          theme(legend.position = "top",
                legend.key.width = unit(1.5, "cm"),
                axis.text.x = element_text(size = text_size),
                axis.text.y = element_blank(),
                axis.title.x = element_text(size = text_size),
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank(),
                panel.background = element_rect(fill = "white"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.text = element_text(size = text_size),
                legend.title = element_text(size = text_size),
                strip.text = element_text(size = text_size + 1),
                strip.background = element_rect(fill = "#b6d3b6"))
        
        # Pull the legend out of the plot as its own grob so it can be placed
        # next to the title instead of stacked above the plot
        legend_grob <- cowplot::get_legend(p_heat)
        
        title_text <- str_wrap(
          paste0("Number of observations each proportion is based on is ",
                 "listed above the colored bar. Bolded red values indicate ",
                 "that there are fewer than the minimum desired number of ",
                 "samples (", min_sample_size, ")."),
          width = 100
        )
        title_grob <- grid::textGrob(
          title_text,
          x = 0, hjust = 0,
          gp = grid::gpar(fontsize = text_size)
        )
        
        # Anchor the legend to the right edge of its cell (instead of the
        # default center-justify, which can clip a wide legend at the
        # plot's right edge)
        legend_panel <- cowplot::ggdraw() +
          cowplot::draw_grob(legend_grob, x = 0.95, hjust = 1)
        
        # Title on the left, legend on the right, in one row
        top_row <- cowplot::plot_grid(title_grob, legend_grob,
                                      nrow = 1, rel_widths = c(1.5, 1))
        
        # That row stacked above the plot (with its own legend removed)
        total_height <- max(300, 250 * length(input$species_order))
        cowplot::plot_grid(top_row,
                           p_heat + theme(legend.position = "none"),
                           ncol = 1, 
                           rel_heights = c(50/total_height, 
                                           (total_height - 50)/total_height))
        
      }
    }
  )
}

# run app ---------------------------------------------------------------------#

shinyApp(ui = ui, server = server)
