#' Create julian date
#' 
#' @param dates vector of dates, provided as "YYYY-MM-DD"
#'
#' @details 
#' 
#' @return Julian Date (numeric)

julian_date <- function(dates) {
  
  # Extract the name of this function for reporting
  function_name <- as.character(match.call())[1]
  
  # Libraries required for this function to work
  dependencies <- "lubridate"
  if (!all(unlist(lapply(X = dependencies, FUN = require, character.only = TRUE)))) {
    stop("At least one package required by ", function_name, 
         " could not be loaded: ", paste(dependencies, collapse = ", "),
         " are required.")
  }

  jd <- as.numeric(as.POSIXct(dates)) # get Unix timestamp
  jd <- jd / 86400 # convert seconds to days
  jd <-jd + 2440587.5 # Adds the Julian Day Number for the Unix epoch (Jan 1, 1970 at midnight UTC)
  jd <- floor(jd)
  
  return(jd)
}
  
  