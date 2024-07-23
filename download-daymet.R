################################################################################
# Download daymet data.
# Script is typically called from model-flowering-onset-dates.R

# Erin Zylstra
# 2024-07-23
################################################################################

require(daymetr)

get_daymet <- function(i) {
  
  temp_lat <- sites[i, ] %>% pull(lat)
  temp_lon <- sites[i, ] %>% pull(lon)
  temp_site <- sites[i, ] %>% pull(site_id)
  
  temp_daymet <- download_daymet(
    lat = temp_lat,
    lon = temp_lon, 
    start = first(yrs),
    end = last(yrs)
  ) %>%
    #--- just get the data part ---#
    .$data %>% 
    #--- assign site_id so we can match with onset data ---#
    mutate(site_id = temp_site) %>%
    #--- get date from day of the year ---#
    mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"))
  
  return(temp_daymet)
}

daymet_df <- lapply(1:nsites, get_daymet) %>%
  bind_rows()

daymet_df <- daymet_df %>%
  rename(prcp_mm = prcp..mm.day.,
         tmin_c = tmin..deg.c.,
         tmax_c = tmax..deg.c.,
         daylength_s = dayl..s.,
         shortwave_radiation_Wm2 = srad..W.m.2.,
         swe_kgm2 = swe..kg.m.2.,
         vapor_press_Pa = vp..Pa.)

# Save csv and put in a zip file. Remove csv and original zip file (it's huge)
write.csv(daymet_df, 
          file = daymet_file,
          row.names = FALSE)
zip_file <- str_replace(daymet_file, ".csv", ".zip")
zip(zip_file, files = daymet_file)
file.remove(daymet_file)
