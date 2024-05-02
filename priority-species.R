################################################################################
# Exploring priority plant species

# Erin Zylstra
# 2024-05-02
################################################################################

library(rnpn)
library(dplyr)

# Load csv with plant species that received votes from one or more states
spp <- read.csv("data/species-with-votes.csv", na.strings = c(NA, ""))

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

# Adding information indicating which species were included in previous work 
# that was done with "16 priority species"
old16 <- read.csv("data/old/previous-16-priority-species.csv")
spp <- spp %>%
  mutate(prev16 = ifelse(common_name %in% old16$common_name, 1, 0))

# For now, classifying things as priority 2 if had multiple states voted species 
# a priority; priority 3 if it was considered one of top 16 priority species in 
# previous iterations of the work; priority 4 if two or more people in a state 
# voted species a priority, and priority 5 are species where only a single 
# person in one state considered the species a priority

spp <- spp %>%
  mutate(priority_level = ifelse(is.na(priority_level), 0, priority_level)) %>%
  mutate(priority = ifelse(priority_level == 1, 1,
                           ifelse(votes_nstates > 1 | revisit_nstates > 1, 2,
                                  ifelse(prev16 == 1, 3,
                                         ifelse(votes_high == 1, 4, 5)))))
# check:  
# count(spp, priority, priority_level, votes_nstates, revisit_nstates, 
#       prev16, votes_high)

count(spp, priority)
  # priority 1: 8
  # priority 2: 15
  # priority 3: 5
  # priority 4: 21
  # priority 5: 45

# Change entries with species = NA to species = "spp."
spp <- spp %>%
  mutate(species = ifelse(is.na(species), "spp.", species)) %>%
  select(-priority_level)

# See how many of the entries with/without a species name are high priority:
count(filter(spp, species != "spp."), priority)
  # If we exclude the spp's:
  # priority 1: 8
  # priority 2: 13
  # priority 3: 5
  # priority 4: 17
  # priority 5: 40
filter(spp, species == "spp.", priority == 2)
  # Excluding Liatris spp (note that L. aspera is in priority 1 category)
  # Excluding Solidago spp (note that S. rugosa is in priority 2 category)

# Create high-priority species list (levels 1:3)
spp_hp <- spp %>%
  filter(species != "spp." & priority %in% 1:3)

# Dataframe with species in NPN database
species_list <- npn_species() %>% 
  rename(NPN_common_name = common_name) %>%
  select(species_id, NPN_common_name, genus, species) %>%
  data.frame()

# Match priority species with NPN species
spp_hp <- spp_hp %>%
  left_join(species_list, by = c("genus", "species"))
filter(spp_hp, is.na(species_id))
# Just one high priority species that doesn't appear in NPN database: 
# Coreopsis tinctoria (golden tickseed). Priority 2: 1 vote in LA and NM

# Simplify state information
spp_hp <- spp_hp %>%
  mutate(LA = ifelse((!is.na(votes_LA) & votes_LA > 0) |
                       (!is.na(revisit_votes_LA) & revisit_votes_LA > 0), 1, 0)) %>%
  mutate(OK = ifelse((!is.na(votes_OK) & votes_OK > 0) |
                       (!is.na(revisit_votes_OK) & revisit_votes_OK > 0), 1, 0)) %>%
  mutate(NM = ifelse((!is.na(votes_NM) & votes_NM > 0) |
                       (!is.na(revisit_votes_NM) & revisit_votes_NM > 0), 1, 0)) %>%
  mutate(TX = ifelse(!is.na(revisit_votes_TX) & revisit_votes_TX > 0, 1, 0)) %>%
  select(common_name, genus, species, species_id, LA, OK, NM, TX, priority)

# Write to file
# write.csv(spp_hp, "data/priority-species.csv", row.names = FALSE)
