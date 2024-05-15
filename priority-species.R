################################################################################
# Exploring priority plant species

# Erin Zylstra
# 2024-05-15
################################################################################

library(rnpn)
library(dplyr)
library(tidyr)

# Load csv with plant species that received votes from one or more states in
# original or "revisited" spreadsheets.
spp <- read.csv("data/species-with-votes.csv", na.strings = c(NA, ""))

# Replace missing species epithets with "spp." and replace NAs with 0s in
# priority_level column
spp <- spp %>%
  mutate(species = ifelse(is.na(species), "spp.", species)) %>%
  mutate(priority_level = replace_na(priority_level, 0))

# Create rvote column indicating if species got any votes in revisited spreadsheet
spp <- spp %>%
  rowwise() %>% 
  mutate(rvote_sum = sum(c_across(starts_with("revisit")))) %>%
  ungroup() %>%
  mutate(rvote = ifelse(!is.na(rvote_sum) & rvote_sum > 0, 1, 0)) %>%
  select(-rvote_sum) %>%
  data.frame()

# Add up the total number of votes (across states) that each species received in 
# original spreadsheets for LA, OK, NM
spp <- spp %>%
  rowwise() %>%
  mutate(nvotes_orig = sum(c_across(starts_with("votes")), na.rm = TRUE)) %>%
  ungroup() %>%
  data.frame()

# Append indicator of which species were included in previous work with 
# "16 priority species"
old16 <- read.csv("data/old/previous-16-priority-species.csv")
old16 <- old16 %>%
  select(genus, species) %>%
  mutate(prev16 = 1)
spp <- spp %>%
  left_join(old16, by = c("genus", "species")) %>%
  mutate(prev16 = replace_na(prev16, 0))

# 8 species have been designated as "top priority" with priority_level == 1.
# Making rules to classify the rest of the species:
# priority 2: received one or more vote in revisited spreadsheet (rvote == 1)
# priority 3: was in previous group of 16 priority species
# priority 4: received multiple votes in original spreadsheets (nvotes_orig > 1)
# priority 5: none of the above.
spp <- spp %>%
  mutate(
    priority = case_when(
      priority_level == 1 ~ 1,
      rvote == 1 ~ 2,
      prev16 == 1 ~ 3,
      nvotes_orig > 1 ~ 4,
      .default = 5
    )
  )
# checks:
count(spp, priority, priority_level, rvote, prev16, nvotes_orig)
count(spp, priority)
  # priority 1: 8 spp
  # priority 2: 19 spp
  # priority 3: 6 spp
  # priority 4: 26 spp
  # priority 5: 36 spp

# Create dataframe with species in NPN database
species_list <- npn_species() %>% 
  select(species_id, common_name, genus, species) %>%
  data.frame()

# Just work with priority 1-4 for now.
hp <- filter(spp, priority < 5)

# Match up priority species with NPN database entries. Make the NPN common name
# the primary one and the common name in state spreadsheets secondary.
hp <- hp %>%
  rename(common_name_states = common_name) %>%
  left_join(species_list, by = c("genus", "species")) %>%
  relocate(common_name_states, .after = "common_name")

# Create separate dataframe for taxa that don't have specific epithets
hp_genus <- filter(hp, species == "spp.")

# For now, just retain genus-level taxa if they were priority 1 or 2 and 
# replace each row with all NPN species in that genus
hp_genus <- hp_genus %>%
  filter(priority < 3) %>%
  select(-c(species, species_id, common_name)) %>%
  left_join(species_list, by = "genus") %>%
  relocate(species, .after = "genus") %>%
  relocate(common_name_states, .after = "common_name")

# Remove spp's from hp dataframe and remove unnecessary column
hp <- hp %>%
  filter(species != "spp.") %>%
  select(-priority_level)

# Write dataframes to file, keeping taxa without specific epithets separate
# write.csv(hp, "data/priority-species.csv", row.names = FALSE)
# write.csv(hp_genus, "data/priority-genera.csv", row.names = FALSE)
