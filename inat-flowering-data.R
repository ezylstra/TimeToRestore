library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(ggplot2)
library(grid)
library(cowplot)
library(terra)
library(tidyterra)
library(rnpn)

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

# Desired minimum sample size for calculating proportions
min_sample_size <- 5

dfw <- dfw %>%
  mutate(wk = week(obsdate)) %>%
  # Remove observations in week 53 (Dec 31 [and Dec 30 in leap years])
  filter(wk < 53) %>%
  mutate(wk2 = ifelse(wk%%2 == 0, wk - 1, wk))

# Load iNat (via Phenobase) data ----------------------------------------------#
# iNat observations with flower annotation (don't have separate open flower
# annotation)

phenofolder <- "data/iNat/phenobase"
phenozips <- list.files(phenofolder, full.names = TRUE)

for (i in seq_along(phenozips)) {
  filename <- unz(phenozips[i], "data.csv")
  pheno1 <- read.csv(filename)
  if (i == 1) {
    inatfull <- pheno1
  } else {
    inatfull <- bind_rows(inatfull, pheno1)
  }
}
rm(pheno1)

# Simplify data
inat <- inatfull %>%
  mutate(date = ymd(eventDate)) %>%
  rename(coordU = coordinateUncertaintyInMeters) %>%
  select(scientificName, year, latitude, longitude, date, coordU)

# Look at what points have tons of geographic uncertainty
# inat %>%
#   mutate(uncertain = ifelse(is.na(coordU) | coordU > 10000, 1, 0)) %>%
#   ggplot() +
#   geom_point(aes(x = longitude, y = latitude, color = factor(uncertain)))
# No obvious patterns

# What does it mean if the coordU field is NA?
# For now will delete these, along with records that have very high values
# (location uncertainty > 10 km)
inat <- inat %>%
  filter(!is.na(latitude)) %>%
  filter(!is.na(coordU) & coordU <= 10000) %>%
  select(-coordU)

# Simplify species names (iNat/phenobase provides subspecies for some)
inat <- inat %>%
  separate_wider_delim(scientificName, 
                       delim = " ", 
                       names = c("genus", "species", "subspecies"),
                       too_few = "align_start")

# Get NPN species_ids and common names
npn_spp <- npn_species() %>%
  select(species_id, common_name, genus, species)
inat <- inat %>%
  left_join(npn_spp, by = c("genus", "species")) %>%
  data.frame()

# Assign states for each plant location
states <- vect("data/states/cb_2017_us_state_500k.shp")
states <- project(states, "epsg:4326")
states <- select(states, STUSPS)
inatlocs <- inat %>%
  select(longitude, latitude) %>%
  vect(crs = "epsg:4326")
statefill <- terra::extract(states, inatlocs)
inat$state <- statefill$STUSPS
# Limit to US states
inat <- inat %>%
  filter(!is.na(state))

# Create month, week, bi-weekly assignments
inat <- inat %>%
  mutate(month = month(date),
         wk = week(date)) %>%
  mutate(wk2 = ifelse(wk%%2 == 0, wk - 1, wk))

# For now, limiting both datasets to 5 species:
spp5 <- c("American beautyberry", "eastern redbud", "wax mallow",
          "blackeyed Susan", "firewheel")
dfw5 <- dfw %>%
  filter(spp %in% spp5)
inat5 <- inat %>%
  filter(common_name %in% spp5)

# Add Texas region to iNat data
regions <- vect("data/tx-regions/tx-regions.shp")
inat5locs <- inat5 %>%
  select(longitude, latitude) %>%
  vect(crs = "epsg:4326")
regionfill <- terra::extract(regions, inat5locs)
inat5$region <- regionfill$region
inat5 <- inat5 %>%
  mutate(scregion = ifelse(state %in% c("LA", "OK", "NM", "TX"), 1, 0)) %>%
  filter(wk2 < 53) %>%
  mutate(wk_doy4 = wk2 * 7) %>%
  mutate(wk_date4 = parse_date_time(x = paste(2025, wk_doy4), 
                                    orders = "yj"))

# Data summaries by region/state
table(dfw5$spp, dfw5$region, useNA = "ifany")
table(dfw5$spp, dfw5$state, useNA = "ifany")
table(inat5$common_name[inat5$scregion == 1], 
      inat5$region[inat5$scregion == 1], 
      useNA = "ifany")
table(inat5$common_name[inat5$scregion == 1], 
      inat5$state[inat5$scregion == 1], 
      useNA = "ifany")

# Create function for size breaks that include low values (for bubble plot)
pretty_breaks <- function(x, n = 3) {
  pr <- pretty(x, n)
  pr <- pr[pr > 0]
  pr <- sort(unique(c(1, min_sample_size - 1, pr)))
  return(pr)
}

# American beautyberry: flowers, bi-weekly data, Texas ------------------------#
ambe_nn <- dfw5 %>%
  filter(spp == "American beautyberry") %>%
  filter(state == "TX") %>%
  rename(status = status_500,
         midpoint = midpoint_500)

ambe_props_temp <- ambe_nn %>%
  mutate(period = wk2) %>%
  mutate(wk_doy4 = period * 7) %>%
  mutate(wk_date4 = parse_date_time(x = paste(2025, wk_doy4), 
                                    orders = "yj")) %>%
  filter(!is.na(status)) %>%
  arrange(id, yr, period, desc(status), desc(midpoint)) %>% 
  distinct(id, yr, period, .keep_all = TRUE) %>%
  group_by(spp, period, wk_date4) %>%
  summarize(nobs = n(),
            nyes = sum(status),
            prop = nyes/nobs,
            .groups = "drop") %>%
  data.frame() %>% 
  mutate(obs_group = ifelse(nobs < min_sample_size, "low", "sufficient"))

p_data <- ambe_props_temp
  # Compute breaks the same way ggplot will: apply pretty_breaks to the
  # data range, then drop anything outside it (mimicking ggplot's behavior)
  data_range <- range(p_data$nobs, na.rm = TRUE)
  br <- pretty_breaks(p_data$nobs)
  br <- br[br >= data_range[1] & br <= data_range[2]]

  aes_fill  <- ifelse(br < min_sample_size, "white", "steelblue3")
  # text_size <- 14

  p_bubble <- p_data %>%
    ggplot(aes(x = period, y = prop)) +
    geom_line(color = "gray50") +
    geom_point(aes(size = nobs, fill = obs_group), color = "steelblue3",
               shape = 21) +
    facet_wrap(~ spp) +
    scale_fill_manual(values = c("low" = "white", "sufficient" = "steelblue3"),
                      guide = "none") +    # suppress separate fill legend
    scale_size_continuous(range = c(2, 7),
                          breaks = br,
                          guide = guide_legend(
                            title = "NPN observations",
                            nrow = 1,
                            override.aes = list(fill = aes_fill))
    ) +
    scale_x_continuous(limits = c(1, 52), expand = 0.05) +
    scale_y_continuous(limits = c(-0.05, 1.05)) +
    labs(y = "Proportion with flowers") +
    theme(legend.position = "top",
          axis.title.x = element_blank())
  
ambe_inat_all <- inat5 %>%
  filter(common_name == "American beautyberry") %>%
  group_by(wk2) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  data.frame() %>%
  mutate(area = "iNaturalist: Rangewide in US")
  
fit_all <- loess(log(nobs) ~ wk2, span = 0.25, data = ambe_inat_all)
smoothed_all <- data.frame(
  wk2 = seq(min(ambe_inat_all$wk2), max(ambe_inat_all$wk2), length = 100)
  ) %>%
  mutate(nobs = exp(predict(fit_all, newdata = .)),
         area = "iNaturalist: Rangewide in US")

ambe_inat_tx <- inat5 %>%
  filter(common_name == "American beautyberry") %>%
  filter(state == "TX") %>%
  group_by(wk2) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  data.frame() %>%
  mutate(area = "iNaturalist: Texas")

fit_tx <- loess(log(nobs) ~ wk2, span = 0.25, data = ambe_inat_tx)
smoothed_tx <- data.frame(
  wk2 = seq(min(ambe_inat_tx$wk2), max(ambe_inat_tx$wk2), length = 100)
  ) %>%
  mutate(nobs = exp(predict(fit_tx, newdata = .)),
         area = "iNaturalist: Texas")

ambe_inat <- rbind(ambe_inat_all,
                   ambe_inat_tx)
smoothed_inat <- rbind(smoothed_all,
                       smoothed_tx)
p_inat <- ambe_inat %>%
  ggplot(aes(x = wk2, y = nobs)) +
  geom_point(size = 2) + 
  geom_line(data = smoothed_inat, color = "steelblue") +
  facet_wrap(~area, ncol = 1, scales = "free_y") +
  labs(x = "Week of year", y = "Number of flower photos")

ambe_plots <- plot_grid(p_bubble, p_inat, ncol = 1, align = "v", 
                        rel_heights = c(0.4, 0.6))
ggsave("output/nn-inat-comparisons/american-beautyberry-tx.png",
       ambe_plots, 
       width = 8, 
       height = 8, 
       units = "in", 
       dpi = 600)

# blackeyed Susan: flowers, bi-weekly data, Texas -----------------------------#
besu_nn <- dfw5 %>%
  filter(spp == "blackeyed Susan") %>%
  filter(state == "TX") %>%
  rename(status = status_500,
         midpoint = midpoint_500)

besu_props_temp <- besu_nn %>%
  mutate(period = wk2) %>%
  mutate(wk_doy4 = period * 7) %>%
  mutate(wk_date4 = parse_date_time(x = paste(2025, wk_doy4), 
                                    orders = "yj")) %>%
  filter(!is.na(status)) %>%
  arrange(id, yr, period, desc(status), desc(midpoint)) %>% 
  distinct(id, yr, period, .keep_all = TRUE) %>%
  group_by(spp, period, wk_date4) %>%
  summarize(nobs = n(),
            nyes = sum(status),
            prop = nyes/nobs,
            .groups = "drop") %>%
  data.frame() %>% 
  mutate(obs_group = ifelse(nobs < min_sample_size, "low", "sufficient"))

p_data <- besu_props_temp
# Compute breaks the same way ggplot will: apply pretty_breaks to the
# data range, then drop anything outside it (mimicking ggplot's behavior)
data_range <- range(p_data$nobs, na.rm = TRUE)
br <- pretty_breaks(p_data$nobs)
br <- br[br >= data_range[1] & br <= data_range[2]]

aes_fill  <- ifelse(br < min_sample_size, "white", "steelblue3")

p_bubble <- p_data %>%
  ggplot(aes(x = period, y = prop)) +
  geom_line(color = "gray50") +
  geom_point(aes(size = nobs, fill = obs_group), color = "steelblue3",
             shape = 21) +
  facet_wrap(~ spp) +
  scale_fill_manual(values = c("low" = "white", "sufficient" = "steelblue3"),
                    guide = "none") +    # suppress separate fill legend
  scale_size_continuous(range = c(2, 7),
                        breaks = br,
                        guide = guide_legend(
                          title = "NPN observations",
                          nrow = 1,
                          override.aes = list(fill = aes_fill))
  ) +
  scale_x_continuous(limits = c(1, 52), expand = 0.05) +
  scale_y_continuous(limits = c(-0.05, 1.05)) +
  labs(y = "Proportion with flowers") +
  theme(legend.position = "top",
        axis.title.x = element_blank())

besu_inat_all <- inat5 %>%
  filter(common_name == "blackeyed Susan") %>%
  group_by(wk2) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  data.frame() %>%
  mutate(area = "iNaturalist: Rangewide in US")

fit_all <- loess(log(nobs) ~ wk2, span = 0.25, data = besu_inat_all)
smoothed_all <- data.frame(
  wk2 = seq(min(besu_inat_all$wk2), max(besu_inat_all$wk2), length = 100)
) %>%
  mutate(nobs = exp(predict(fit_all, newdata = .)),
         area = "iNaturalist: Rangewide in US")

besu_inat_tx <- inat5 %>%
  filter(common_name == "blackeyed Susan") %>%
  filter(state == "TX") %>%
  group_by(wk2) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  data.frame() %>%
  mutate(area = "iNaturalist: Texas")

fit_tx <- loess(log(nobs) ~ wk2, span = 0.25, data = besu_inat_tx)
smoothed_tx <- data.frame(
  wk2 = seq(min(besu_inat_tx$wk2), max(besu_inat_tx$wk2), length = 100)
) %>%
  mutate(nobs = exp(predict(fit_tx, newdata = .)),
         area = "iNaturalist: Texas")

besu_inat <- rbind(besu_inat_all,
                   besu_inat_tx)
smoothed_inat <- rbind(smoothed_all,
                       smoothed_tx)
p_inat <- besu_inat %>%
  ggplot(aes(x = wk2, y = nobs)) +
  geom_point(size = 2) + 
  geom_line(data = smoothed_inat, color = "steelblue") +
  facet_wrap(~area, ncol = 1, scales = "free_y") +
  labs(x = "Week of year", y = "Number of flower photos")

besu_plots <- plot_grid(p_bubble, p_inat, ncol = 1, align = "v", 
                        rel_heights = c(0.4, 0.6))
ggsave("output/nn-inat-comparisons/blackeyed-Susan-tx.png",
       besu_plots, 
       width = 8, 
       height = 8, 
       units = "in", 
       dpi = 600)

# firewheel: flowers, bi-weekly data, Texas -----------------------------#
fire_nn <- dfw5 %>%
  filter(spp == "firewheel") %>%
  filter(state == "TX") %>%
  rename(status = status_500,
         midpoint = midpoint_500)

fire_props_temp <- fire_nn %>%
  mutate(period = wk2) %>%
  mutate(wk_doy4 = period * 7) %>%
  mutate(wk_date4 = parse_date_time(x = paste(2025, wk_doy4), 
                                    orders = "yj")) %>%
  filter(!is.na(status)) %>%
  arrange(id, yr, period, desc(status), desc(midpoint)) %>% 
  distinct(id, yr, period, .keep_all = TRUE) %>%
  group_by(spp, period, wk_date4) %>%
  summarize(nobs = n(),
            nyes = sum(status),
            prop = nyes/nobs,
            .groups = "drop") %>%
  data.frame() %>% 
  mutate(obs_group = ifelse(nobs < min_sample_size, "low", "sufficient"))

p_data <- fire_props_temp
# Compute breaks the same way ggplot will: apply pretty_breaks to the
# data range, then drop anything outside it (mimicking ggplot's behavior)
data_range <- range(p_data$nobs, na.rm = TRUE)
br <- pretty_breaks(p_data$nobs)
br <- br[br >= data_range[1] & br <= data_range[2]]

aes_fill  <- ifelse(br < min_sample_size, "white", "steelblue3")

p_bubble <- p_data %>%
  ggplot(aes(x = period, y = prop)) +
  geom_line(color = "gray50") +
  geom_point(aes(size = nobs, fill = obs_group), color = "steelblue3",
             shape = 21) +
  facet_wrap(~ spp) +
  scale_fill_manual(values = c("low" = "white", "sufficient" = "steelblue3"),
                    guide = "none") +    # suppress separate fill legend
  scale_size_continuous(range = c(2, 7),
                        breaks = br,
                        guide = guide_legend(
                          title = "NPN observations",
                          nrow = 1,
                          override.aes = list(fill = aes_fill))
  ) +
  scale_x_continuous(limits = c(1, 52), expand = 0.05) +
  scale_y_continuous(limits = c(-0.05, 1.05)) +
  labs(y = "Proportion with flowers") +
  theme(legend.position = "top",
        axis.title.x = element_blank())

fire_inat_all <- inat5 %>%
  filter(common_name == "firewheel") %>%
  group_by(wk2) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  data.frame() %>%
  mutate(area = "iNaturalist: Rangewide in US")

fit_all <- loess(log(nobs) ~ wk2, span = 0.25, data = fire_inat_all)
smoothed_all <- data.frame(
  wk2 = seq(min(fire_inat_all$wk2), max(fire_inat_all$wk2), length = 100)
) %>%
  mutate(nobs = exp(predict(fit_all, newdata = .)),
         area = "iNaturalist: Rangewide in US")

fire_inat_tx <- inat5 %>%
  filter(common_name == "firewheel") %>%
  filter(state == "TX") %>%
  group_by(wk2) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  data.frame() %>%
  mutate(area = "iNaturalist: Texas")

fit_tx <- loess(log(nobs) ~ wk2, span = 0.25, data = fire_inat_tx)
smoothed_tx <- data.frame(
  wk2 = seq(min(fire_inat_tx$wk2), max(fire_inat_tx$wk2), length = 100)
) %>%
  mutate(nobs = exp(predict(fit_tx, newdata = .)),
         area = "iNaturalist: Texas")

fire_inat <- rbind(fire_inat_all,
                   fire_inat_tx)
smoothed_inat <- rbind(smoothed_all,
                       smoothed_tx)
p_inat <- fire_inat %>%
  ggplot(aes(x = wk2, y = nobs)) +
  geom_point(size = 2) + 
  geom_line(data = smoothed_inat, color = "steelblue") +
  facet_wrap(~area, ncol = 1, scales = "free_y") +
  labs(x = "Week of year", y = "Number of flower photos")

fire_plots <- plot_grid(p_bubble, p_inat, ncol = 1, align = "v", 
                        rel_heights = c(0.4, 0.6))
ggsave("output/nn-inat-comparisons/firewheel-tx.png",
       fire_plots, 
       width = 8, 
       height = 8, 
       units = "in", 
       dpi = 600)

# eastern redbud: flowers, bi-weekly data, Texas -----------------------------#
eare_nn <- dfw5 %>%
  filter(spp == "eastern redbud") %>%
  filter(state == "TX") %>%
  rename(status = status_500,
         midpoint = midpoint_500)

eare_props_temp <- eare_nn %>%
  mutate(period = wk2) %>%
  mutate(wk_doy4 = period * 7) %>%
  mutate(wk_date4 = parse_date_time(x = paste(2025, wk_doy4), 
                                    orders = "yj")) %>%
  filter(!is.na(status)) %>%
  arrange(id, yr, period, desc(status), desc(midpoint)) %>% 
  distinct(id, yr, period, .keep_all = TRUE) %>%
  group_by(spp, period, wk_date4) %>%
  summarize(nobs = n(),
            nyes = sum(status),
            prop = nyes/nobs,
            .groups = "drop") %>%
  data.frame() %>% 
  mutate(obs_group = ifelse(nobs < min_sample_size, "low", "sufficient"))

p_data <- eare_props_temp
# Compute breaks the same way ggplot will: apply pretty_breaks to the
# data range, then drop anything outside it (mimicking ggplot's behavior)
data_range <- range(p_data$nobs, na.rm = TRUE)
br <- pretty_breaks(p_data$nobs)
br <- br[br >= data_range[1] & br <= data_range[2]]

aes_fill  <- ifelse(br < min_sample_size, "white", "steelblue3")

p_bubble <- p_data %>%
  ggplot(aes(x = period, y = prop)) +
  geom_line(color = "gray50") +
  geom_point(aes(size = nobs, fill = obs_group), color = "steelblue3",
             shape = 21) +
  facet_wrap(~ spp) +
  scale_fill_manual(values = c("low" = "white", "sufficient" = "steelblue3"),
                    guide = "none") +    # suppress separate fill legend
  scale_size_continuous(range = c(2, 7),
                        breaks = br,
                        guide = guide_legend(
                          title = "NPN observations",
                          nrow = 1,
                          override.aes = list(fill = aes_fill))
  ) +
  scale_x_continuous(limits = c(1, 52), expand = 0.05) +
  scale_y_continuous(limits = c(-0.05, 1.05)) +
  labs(y = "Proportion with flowers") +
  theme(legend.position = "top",
        axis.title.x = element_blank())

eare_inat_all <- inat5 %>%
  filter(common_name == "eastern redbud") %>%
  group_by(wk2) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  data.frame() %>%
  mutate(area = "iNaturalist: Rangewide in US")

fit_all <- loess(log(nobs) ~ wk2, span = 0.25, data = eare_inat_all)
smoothed_all <- data.frame(
  wk2 = seq(min(eare_inat_all$wk2), max(eare_inat_all$wk2), length = 100)
) %>%
  mutate(nobs = exp(predict(fit_all, newdata = .)),
         area = "iNaturalist: Rangewide in US")

eare_inat_tx <- inat5 %>%
  filter(common_name == "eastern redbud") %>%
  filter(state == "TX") %>%
  group_by(wk2) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  data.frame() %>%
  mutate(area = "iNaturalist: Texas")

fit_tx <- loess(log(nobs) ~ wk2, span = 0.25, data = eare_inat_tx)
smoothed_tx <- data.frame(
  wk2 = seq(min(eare_inat_tx$wk2), max(eare_inat_tx$wk2), length = 100)
) %>%
  mutate(nobs = exp(predict(fit_tx, newdata = .)),
         area = "iNaturalist: Texas")

eare_inat <- rbind(eare_inat_all,
                   eare_inat_tx)
smoothed_inat <- rbind(smoothed_all,
                       smoothed_tx)
p_inat <- eare_inat %>%
  ggplot(aes(x = wk2, y = nobs)) +
  geom_point(size = 2) + 
  geom_line(data = smoothed_inat, color = "steelblue") +
  facet_wrap(~area, ncol = 1, scales = "free_y") +
  labs(x = "Week of year", y = "Number of flower photos")

eare_plots <- plot_grid(p_bubble, p_inat, ncol = 1, align = "v", 
                        rel_heights = c(0.4, 0.6))
ggsave("output/nn-inat-comparisons/eastern-redbud-tx.png",
       eare_plots, 
       width = 8, 
       height = 8, 
       units = "in", 
       dpi = 600)

# wax mallow: flowers, bi-weekly data, Texas -----------------------------#
wama_nn <- dfw5 %>%
  filter(spp == "wax mallow") %>%
  filter(state == "TX") %>%
  rename(status = status_500,
         midpoint = midpoint_500)

wama_props_temp <- wama_nn %>%
  mutate(period = wk2) %>%
  mutate(wk_doy4 = period * 7) %>%
  mutate(wk_date4 = parse_date_time(x = paste(2025, wk_doy4), 
                                    orders = "yj")) %>%
  filter(!is.na(status)) %>%
  arrange(id, yr, period, desc(status), desc(midpoint)) %>% 
  distinct(id, yr, period, .keep_all = TRUE) %>%
  group_by(spp, period, wk_date4) %>%
  summarize(nobs = n(),
            nyes = sum(status),
            prop = nyes/nobs,
            .groups = "drop") %>%
  data.frame() %>% 
  mutate(obs_group = ifelse(nobs < min_sample_size, "low", "sufficient"))

p_data <- wama_props_temp
# Compute breaks the same way ggplot will: apply pretty_breaks to the
# data range, then drop anything outside it (mimicking ggplot's behavior)
data_range <- range(p_data$nobs, na.rm = TRUE)
br <- pretty_breaks(p_data$nobs)
br <- br[br >= data_range[1] & br <= data_range[2]]

aes_fill  <- ifelse(br < min_sample_size, "white", "steelblue3")

p_bubble <- p_data %>%
  ggplot(aes(x = period, y = prop)) +
  geom_line(color = "gray50") +
  geom_point(aes(size = nobs, fill = obs_group), color = "steelblue3",
             shape = 21) +
  facet_wrap(~ spp) +
  scale_fill_manual(values = c("low" = "white", "sufficient" = "steelblue3"),
                    guide = "none") +    # suppress separate fill legend
  scale_size_continuous(range = c(2, 7),
                        breaks = br,
                        guide = guide_legend(
                          title = "NPN observations",
                          nrow = 1,
                          override.aes = list(fill = aes_fill))
  ) +
  scale_x_continuous(limits = c(1, 52), expand = 0.05) +
  scale_y_continuous(limits = c(-0.05, 1.05)) +
  labs(y = "Proportion with flowers") +
  theme(legend.position = "top",
        axis.title.x = element_blank())

wama_inat_all <- inat5 %>%
  filter(common_name == "wax mallow") %>%
  group_by(wk2) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  data.frame() %>%
  mutate(area = "iNaturalist: Rangewide in US")

fit_all <- loess(log(nobs) ~ wk2, span = 0.25, data = wama_inat_all)
smoothed_all <- data.frame(
  wk2 = seq(min(wama_inat_all$wk2), max(wama_inat_all$wk2), length = 100)
) %>%
  mutate(nobs = exp(predict(fit_all, newdata = .)),
         area = "iNaturalist: Rangewide in US")

wama_inat_tx <- inat5 %>%
  filter(common_name == "wax mallow") %>%
  filter(state == "TX") %>%
  group_by(wk2) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  data.frame() %>%
  mutate(area = "iNaturalist: Texas")

fit_tx <- loess(log(nobs) ~ wk2, span = 0.25, data = wama_inat_tx)
smoothed_tx <- data.frame(
  wk2 = seq(min(wama_inat_tx$wk2), max(wama_inat_tx$wk2), length = 100)
) %>%
  mutate(nobs = exp(predict(fit_tx, newdata = .)),
         area = "iNaturalist: Texas")

wama_inat <- rbind(wama_inat_all,
                   wama_inat_tx)
smoothed_inat <- rbind(smoothed_all,
                       smoothed_tx)
p_inat <- wama_inat %>%
  ggplot(aes(x = wk2, y = nobs)) +
  geom_point(size = 2) + 
  geom_line(data = smoothed_inat, color = "steelblue") +
  facet_wrap(~area, ncol = 1, scales = "free_y") +
  labs(x = "Week of year", y = "Number of flower photos")

wama_plots <- plot_grid(p_bubble, p_inat, ncol = 1, align = "v", 
                        rel_heights = c(0.4, 0.6))
ggsave("output/nn-inat-comparisons/wax-mallow-tx.png",
       wama_plots, 
       width = 8, 
       height = 8, 
       units = "in", 
       dpi = 600)



  





# Extract observations for 4-state area
inat4 <- inat %>%
  filter(state %in% c("NM", "OK", "TX", "LA"))

inat_month <- inat %>%
  group_by(common_name, month) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  group_by(common_name) %>%
  mutate(nmonths = n_distinct(month),
         totalobs = sum(nobs)) %>%
  ungroup() %>%
  filter(totalobs > 100) %>%
  data.frame()

# ggplot(inat_month, aes(x = month, y = nobs)) +
#   geom_smooth(method = "loess", span = 0.35, se = FALSE) +
#   geom_point() +
#   scale_x_continuous(breaks = seq(1, 11, by = 2)) +
#   scale_y_log10() +
#   facet_wrap(~common_name, scales = "free_y")

smoothed_month <- inat_month %>%
  group_by(common_name) %>%
  group_modify(function(.x, .y) {
    fit <- loess(log(nobs) ~ month, data = .x, span = 0.35)
    tibble(month = seq(min(.x$month), max(.x$month), length.out = 100)) %>%
      mutate(nobs = exp(predict(fit, newdata = .)))
  })

ggplot(inat_month, aes(x = month, y = nobs)) +
  geom_point() +
  geom_line(data = smoothed_month, color = "steelblue") +
  facet_wrap(~common_name, scales = "free_y")

# Biweekly
inat_biweek <- inat4 %>%
  group_by(common_name, wk2) %>%
  summarize(nobs = n(), .groups = "drop") %>%
  group_by(common_name) %>%
  mutate(nperiods = n_distinct(wk2),
         totalobs = sum(nobs)) %>%
  ungroup() %>%
  filter(totalobs > 100) %>%
  data.frame()

smoothed_biweek <- inat_biweek %>%
  group_by(common_name) %>%
  group_modify(function(.x, .y) {
    fit <- loess(log(nobs) ~ wk2, data = .x, span = 0.25)
    tibble(wk2 = seq(min(.x$wk2), max(.x$wk2), length.out = 100)) %>%
      mutate(nobs = exp(predict(fit, newdata = .)))
  })

ggplot(inat_biweek, aes(x = wk2, y = nobs)) +
  geom_point() +
  geom_line(data = smoothed_biweek, color = "steelblue") +
  facet_wrap(~common_name, scales = "free_y")

# Weekly (seems like there's too much noise at weekly timescale)
# inat_week <- inat4 %>%
#   group_by(common_name, wk) %>%
#   summarize(nobs = n(), .groups = "drop") %>%
#   group_by(common_name) %>%
#   mutate(nweeks = n_distinct(wk),
#          totalobs = sum(nobs)) %>%
#   ungroup() %>%
#   filter(totalobs > 100) %>%
#   data.frame()
# 
# smoothed_week <- inat_week %>%
#   group_by(common_name) %>%
#   group_modify(function(.x, .y) {
#     fit <- loess(log(nobs) ~ wk, data = .x, span = 0.35)
#     tibble(wk = seq(min(.x$wk), max(.x$wk), length.out = 100)) %>%
#       mutate(nobs = exp(predict(fit, newdata = .)))
#   })
# 
# ggplot(inat_week, aes(x = wk, y = nobs)) +
#   geom_point() +
#   geom_line(data = smoothed_week, color = "steelblue") +
#   facet_wrap(~common_name, scales = "free_y")
