################################################################################
# Exploring priority plant species

# Erin Zylstra
# 2024-04-29
################################################################################

library(rnpn)
library(dplyr)
# library(httr)
# library(readr)
# library(lubridate)
# library(ggplot2)

rm(list = ls())

# Load csv with priority species
spp <- read.csv("data/priority_species.csv", na.strings = c(NA, ""))
# 8 species have been designated as "top priority" with priority_level == 1
# Will make some rules to determine priority level for the rest
spp <- spp %>%
  mutate(rvote_LA = ifelse(revisit_votes_LA > 0, 1, 0),
         rvote_OK = ifelse(revisit_votes_OK > 0, 1, 0),
         rvote_NM = ifelse(revisit_votes_NM > 0, 1, 0),
         rvote_TX = ifelse(revisit_votes_TX > 0, 1, 0)) %>%
  mutate(revisit_nstates = ifelse(is.na(revisit_votes_LA), 0,
                                  rvote_LA + rvote_NM + rvote_OK + rvote_TX)) %>%
  mutate(v_LA = ifelse(is.na(votes_LA), 0, votes_LA),
         v_OK = ifelse(is.na(votes_OK), 0, votes_OK),
         v_NM = ifelse(is.na(votes_NM), 0, votes_NM)) %>%
  mutate(votes_high = ifelse(v_LA > 1 | v_OK > 1 | v_NM > 1, 1, 0)) %>%
  mutate(v_LA = ifelse(v_LA > 0, 1, 0),
         v_OK = ifelse(v_OK > 0, 1, 0),
         v_NM = ifelse(v_NM > 0, 1, 0)) %>%
  mutate(votes_nstates = v_LA + v_OK + v_NM) %>%
  select(-c(rvote_LA, rvote_OK, rvote_NM, rvote_TX, v_LA, v_OK, v_NM))
# checks:
# select(spp, c(common_name, genus, species, votes_LA, votes_OK, votes_NM,
#               votes_high, votes_nstates))
# select(spp, c(common_name, genus, species, revisit_votes_LA, revisit_votes_OK, 
#               revisit_votes_NM, revisit_votes_TX, revisit_nstates)) %>%
#   count(revisit_votes_LA, revisit_votes_OK, revisit_votes_NM, revisit_votes_TX, 
#         revisit_nstates)

count(filter(spp, is.na(priority_level)), 
             votes_high, votes_nstates, revisit_nstates)

spp <- spp %>%
  mutate(priority_level = ifelse(is.na(priority_level), 0, priority_level)) %>%
  mutate(priority = ifelse(priority_level == 1, 1,
                           ifelse(votes_high == 1 | votes_nstates > 1 | revisit_nstates > 1, 
                                  2, 3)))
# check:  
# count(spp, votes_high, votes_nstates, revisit_nstates, priority_level, priority)

# If we classify things as priority 2 if had multiple states voting species as
# a priority OR if 2 or more people in a state voted species as priority
count(spp, priority)
  # priority 1: 8
  # priority 2: 41
  # priority 3: 45

#------------------------------------------------------------------------------#
# Pick up here...

# Next step: See if plant names are correct and/or they match up with entries
# in NPN database
# Dataframe with all plant species in NPN database
species_list <- npn_species(kingdom = "Plantae") %>% data.frame()
