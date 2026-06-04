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
# library(DT)
# library(glue)
# library(shinycssloaders)

# Load/format observation data ------------------------------------------------#

df <- read.csv("data/status-intensity-flowers-May2026.csv")
df <- df %>%
  mutate(obsdate = ymd(obsdate))
# For now limit to 4 state area
df <- filter(df, state4 == 1)

# Create table with location options
loc_options <- distinct(df, state, region) %>%
  mutate(region = ifelse(is.na(region), "Statewide", region))

# ui --------------------------------------------------------------------------#

ui <- page_navbar(
  title = "Time to Restore",
  sidebar = sidebar(
    width = "20%",
    div(style = "font-size:90%",
        navset_tab(
          nav_panel(
            title = "Info",
            br(),
            HTML("This app uses data submitted to <i>USA-NPN Nature's Notebook</i>
             to ...")
          ),
          nav_panel(
            title = "Settings",
            br(),
            selectInput(inputId = "state",
                        label = "State",
                        choices = unique(loc_options$state),
                        selected = "TX"),
            uiOutput("loc"),
            radioButtons(inputId = "php", 
                         label = "Phenophase",
                         choices = c("Flowers" = "flower", 
                                     "Open flowers" = "open flower"),
                         inline = TRUE),
            radioButtons(inputId = "vistype",
                         label = "Visualization type",
                         choices = c("Bar chart",
                                     "Bubble plot",
                                     "Heat map")),
            uiOutput("species"),
            selectInput(inputId = "years",
                        label = "Years",
                        choices = c("2025-2026 (combined)", "2025", "2026"),
                        selected = "2025-2026 (combined)"),
            radioButtons(inputId = "period",
                         label = "Summarization period",
                         choices = c("Biweekly" = "biweekly",
                                     "Weekly" = "weekly"))
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

server <- function(input, output) {

  loc_table <- reactive({
    loc_options %>% filter(state == input$state)
  })

  output$loc <- renderUI({
    selectizeInput("locInput", "Location",
                   choices = unique(loc_table()$region))
  })

}

# run app ---------------------------------------------------------------------#

shinyApp(ui = ui, server = server)
