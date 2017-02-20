# This script reads the very large eBird file
#===============================================================================*
# ---- SET-UP ----
#===============================================================================*

# Smart installer will check list of packages that are installed, install any
# necessary package that is missing, and load the library:

smartInstallAndLoad <- function(packageVector){
  for(i in 1:length(packageVector)){
    package <- packageVector[i]
    if(!package %in% rownames(installed.packages())){
      install.packages(packageVector[i],repos="http://cran.rstudio.com/", dependencies=TRUE)
    }
  }
  lapply(packageVector, library, character.only = TRUE)
}

# Load and potentially install libraries:

smartInstallAndLoad(c('stringr', 'dplyr','lubridate', 'readr',
                      'R.utils', 'plyr'))

# Override some libraries for tidyverse functions:

filter <- dplyr::filter
select <- dplyr::select
round_any <- plyr::round_any
rename <- dplyr::rename

#-------------------------------------------------------------------------------*
# ---- FUNCTION TO READ EBIRD DATA ----
#-------------------------------------------------------------------------------*

# This function reads eBird list data in a sets of a given number of rows. Users
# may wish to select different columns.

splitReader <- function(pathToFile, rowSkip, maxRows){
  data <- read.delim(
    pathToFile, header = FALSE, sep = '\t', quote = '',
    col.names = colNames,  skip = rowSkip, nrows = maxRows,
    stringsAsFactors  = FALSE) %>%
    # Modify the below if you want different columns:
    rename(observationID = SAMPLING.EVENT.IDENTIFIER,
                  observer = OBSERVER.ID, 
                  lat = LATITUDE,  lon = LONGITUDE,
                  localityType = LOCALITY.TYPE,
                  protocol = PROTOCOL.TYPE,
                  date = OBSERVATION.DATE, time =  TIME.OBSERVATIONS.STARTED,
                  sp = SCIENTIFIC.NAME,
                  count = OBSERVATION.COUNT,
                  durMinutes = DURATION.MINUTES,
                  effortDist = EFFORT.DISTANCE.KM,
                  nObservers = NUMBER.OBSERVERS) %>%
    select(observationID, observer, lat, lon,
           localityType, protocol,
           date, time, sp, count,
           durMinutes, effortDist, nObservers) %>%
    filter(
      localityType %in% c('H','P'), # points selected from list (H) or on map (P)
      (str_detect(protocol, 'Traveling')|str_detect(protocol, 'Rusty'))
    ) %>%
    select(-c(localityType)) %>%
    tbl_df
  return(data)
}

#-------------------------------------------------------------------------------*
# ---- GET RUSTY OBSERVATIONS ----
#-------------------------------------------------------------------------------*

# Set path to file (examples, do not run):

pathToFile <- 'C:/Users/Brian/Desktop/ebd_rusbla_relMay-2016/ebd_rusbla_relMay-2016.txt'

pathToFile <- 'C:/Users/Brian/Desktop/rubl_summer_2016/ebd_US_200601_201605_relMay-2016/ebd_US_200601_201605_relMay-2016.txt'

# Get vector of column names:

colNames <- read.delim(
  pathToFile, skip = 0, nrows = 1,
  header = TRUE) %>% names 

# Get the total number of rows in the file (this will take a while to run):

tRows <- countLines(pathToFile)[1]

# Get the number of records per set:

rowsPerRead <- 1E7

nSets <- plyr::round_any(tRows/rowsPerRead, 1, f = ceiling)

# Read in tRows at a time 

outList <- vector('list', length = nSets)

for(i in 1:nSets){
  if(i == 1){
    rowsToSkip <- 1
  } else {
    rowsToSkip <- (i-1) * rowsPerRead
  }
  outList[[i]] <- splitReader(pathToFile, rowsToSkip, rowsPerRead)
}

# Convert the list to a data frame:

rublObservations <- bind_rows(outList) %>%
  distinct

# Save file:

saveRDS(rublObservations, file = 'rublObservations.RDS')

#-------------------------------------------------------------------------------*
#---- GET SAMPLING EVENTS ---- 
#-------------------------------------------------------------------------------*

# Get the unique eBird lists (unique user and date):

samplingEvensts <- rublObservations %>%
  select(-c(sp, count)) %>%
  distinct %>%
  mutate(date = as.Date(date)) %>%
  filter(month(date) %in% 1:5) %>%
  distinct

# Save file:

saveRDS(samplingEvents, 'samplingEvents.RDS')
