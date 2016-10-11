# Make samples with data file for spring blitz
#===================================================================================================*
# ---- SET-UP ----
#===================================================================================================*

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

smartInstallAndLoad(c('dplyr', 'tidyr','stringi','stringr', 'sp',
                      'lubridate', 'raster', 'dismo', 'ggplot2'))

# note: 

# setwd('C:/Users/Brian/Desktop/gits/RUBL/rubl_winter/') # Helm
# setwd('C:/Users/Default.Default-THINK/Desktop/gits/RUBL/rubl_winter') # Thinkpad

options(stringsAsFactors = FALSE)

# Set path to land cover data:

# pathToRasterData <- 'C:/Users/Brian/Dropbox/rubl_12_15/'  # Helm
# pathToRasterData <- '/Users/bsevans/Dropbox/rubl_12_15/'   # MacBook Air
# pathToRasterData <- 'C:/Users/Default.Default-THINK/Dropbox/rubl_12_15/' # Thinkpad

# Set path to eBird list data:

# pathToEbirdListData <- 'C:/Users/Brian/Dropbox/eBirdListData.csv'   # Helm
# pathToEbirdListData <- 'C:/Users/Default.Default-THINK/Dropbox/eBirdListData.csv' # ThinkPad
# pathToEbirdListData <- '/Users/bsevans/Dropbox/eBirdListData.csv'   # MacBook Air

# Set path to climate rasters (if available):
# 
# rasterDirTmin <- 'C:/Users/Brian/Desktop/rubl_summer_2016/minTempRasters/tmin_' # Helm
# rasterDirPPt <- 'C:/Users/Brian/Desktop/rubl_summer_2016/pptRasters/ppt_' # Helm

#===================================================================================================*
# ---- BIOCLIMATIC DATA ----
#===================================================================================================*

#---------------------------------------------------------------------------------------------------*
# ---- GET LAND COVER DATA  ----
#---------------------------------------------------------------------------------------------------*

loadEnv <- function(rasterDirectory){
  require(raster)
  # Find the raster data (searches for and ID's all files that end in ".asc":
  
  raster_data <- list.files(rasterDirectory,
                            pattern='\\.asc$', full=T)    
  
  # Create a raster stack of raster layers:
  
  env.stack <- stack(raster_data)
  
  # Add raster stack values to memory:
  
  values(env.stack) <- getValues(env.stack)
  
  # Add projection information to the raster stack:
  
  projection(env.stack) <- CRS('+proj=longlat +datum=WGS84')
  
  names(env.stack) <- c('dev_hi','dev_li','flood','forh', 'form', 'grass',
                        'pasture','ppt','rowcrop', 'shrub','tmin', 'upfor',
                        'weth', 'wetw', 'woodland')
  return(env.stack)
}

#***************************************************************************************************
# RUN SECTION!
#***************************************************************************************************

# Get raster data:

rStack <- loadEnv(paste0(pathToRasterData, 'lc_asc'))

# Remove minimum temperature and precipitation:

rStack <- rStack[[-c(8,11)]]

# Raster data to be used for addresses:

r <- rStack[[1]]

projInfo <- raster::projection(r)

#***************************************************************************************************
#***************************************************************************************************

#---------------------------------------------------------------------------------------------------*
# ----  FUNCTIONS TO DOWNLOAD CLIMATE DATA, IF NECESSARY ----
#---------------------------------------------------------------------------------------------------*

# Function to download tmin raster for a given date, if necessary:

downloadTminRaster <- function(year, month, day){
  yearMonth <- paste0(year, '0', month, '00') %>% as.numeric
  date <- yearMonth + day
  if(year != 2016){
    ftpAddress <- paste('ftp://prism.oregonstate.edu/daily/tmin',
                        year,
                        'PRISM_tmin_stable_4kmD1',
                        sep = '/')
  } else {
    ftpAddress <- paste('ftp://prism.oregonstate.edu/daily/tmin',
                        year,
                        'PRISM_tmin_provisional_4kmD1',
                        sep = '/')
  }
  tmin <- tempfile()
  ftp <- paste(ftpAddress, date, 'bil.zip', sep = '_')
  download.file(ftp, tmin)
  # Unzip:
  tmin <- unzip(tmin)
  # Convert to raster format:
  tmin <- raster(tmin[1])
  # Provide projection information:
  crs(tmin) <- "+proj=longlat +datum=WGS84"
  return(crop(tmin, rStack[[1]]))
}


# Function to download climate data for a given data, if necessary:

downloadPPTRaster <- function(year, month, day){
  yearMonth <- paste0(year, '0', month, '00') %>% as.numeric
  date <- yearMonth + day
  if(year != 2016){
    ftpAddress <- paste('ftp://prism.oregonstate.edu/daily/ppt',
                        year,
                        'PRISM_ppt_stable_4kmD2',
                        sep = '/')
  } else {
    ftpAddress <- paste('ftp://prism.oregonstate.edu/daily/ppt',
                        year,
                        'PRISM_ppt_provisional_4kmD2',
                        sep = '/')
  }
  ppt <- tempfile()
  ftp <- paste(ftpAddress, date, 'bil.zip', sep = '_')
  download.file(ftp, ppt)
  # Unzip:
  ppt <- unzip(ppt)
  # Convert to raster format:
  ppt <- raster(ppt[1])
  # Provide projection information:
  crs(ppt) <- "+proj=longlat +datum=WGS84"
  return(crop(ppt, rStack[[1]]))
}

# Example usage (not run):
# 
# downloadTminRaster(year = '2016', month = '03', day = '15')
# downloadPPTRaster(year = '2016', month = '03', day = '15')

#===================================================================================================*
# ---- PREPARE BIRD OBSERVATIONS ----
#===================================================================================================*

#---------------------------------------------------------------------------------------------------*
# ---- GET LIST DATA FOR RUSTY AND BACKGROUND POINTS ----
#---------------------------------------------------------------------------------------------------*

# Rusty observations:

rustyListsSpring <- read.csv('rublEbird.csv') %>%
  tbl_df %>%
  # Filter to study extent:
    filter(lon > extent(r)[1] & lon < extent(r)[2],
         lat > extent(r)[3] & lat < extent(r)[4]) %>%
  # Subset to dates associated with the spring blitz:
  mutate(date = as.Date(date)) %>%
  mutate(month = lubridate::month(date),
         day = lubridate::day(date)) %>%
  dplyr::filter(
    year(date) > 2013,
    month >= 3 & month < 6,
    !(month == 5 & day > 10)
  ) %>%
  dplyr::select(-c(month, day)) %>%
  # Remove counts recorded as 'X':
  filter(count != 'X') %>%
  mutate(count = as.numeric(count))

# Make a vector of rusty lists where count is recorded as 'X':

rustyXobservations <- read.csv('rublEbird.csv') %>%
  filter(count == 'X') %>%
  .$observationID

# Get eBird list data:

eBirdListsSpring <- read.csv(pathToEbirdListData) %>%
  tbl_df %>%
  # Filter to the study extent:
  filter(lon > extent(r)[1] & lon < extent(r)[2],
         lat > extent(r)[3] & lat < extent(r)[4]) %>%
  # Subset lists to only dates associated with the rusty observations:
  filter(date %in%
           (rustyListsSpring$date %>% 
              unique %>% 
              as.character)) %>%
  # Remove observationIDs where rusty count was reported as X:
  filter(!observationID %in% rustyXobservations) %>%
  # Add count data:
  left_join(rustyListsSpring %>%
              dplyr::select(observationID, count),
            by = 'observationID') %>%
  # Change na counts (no observation match) to 0:
  mutate(count = ifelse(is.na(count), 0, count))

# Add cell addresses:

eBirdListsSpring$cellAddress <- cellFromXY(
  r,
  eBirdListsSpring %>%
    dplyr::select(lon, lat) %>%
    data.frame %>%
    SpatialPoints(proj4string = CRS(projInfo)) 
)

# Extract raster land cover data by cell ID:

envByCell <- data.frame(
  cellAddress = (eBirdListsSpring$cellAddress %>% unique),
  raster::extract(
    x = rStack,
    y = (eBirdListsSpring$cellAddress %>% unique),
    df = TRUE)
) %>%
  tbl_df %>%
  dplyr::select(-c(ID, tmin, ppt))

# Join sampling data to environment data and add year field:

eBirdSamplingEnvSpring <- left_join(
  eBirdListsSpring,
  envByCell,
  by = 'cellAddress'
) %>%
  mutate(date = as.Date(date),
         year = lubridate::year(date)) %>%
  filter(year %in% 2014:2016)

#---------------------------------------------------------------------------------------------------*
# ---- ADD CLIMATE DATA TO POINTS ----
#---------------------------------------------------------------------------------------------------*

# Function to load and prepare climate data:

rasterPrepTminPpt <- function(rasterDirTmin, rasterDirPPt, date){
  # Get ppt and tmin rasters for given date:
  tminR <- raster(paste0(rasterDirTmin, date))
  pptR <- raster(paste0(rasterDirPPt, date))
  
  # Create a raster stack of raster layers:
  
  tminPptStack = stack(tminR, pptR)
  
  # Add raster stack values to memory:
  
  values(tminPptStack) = getValues(tminPptStack)
  
  # Add projection information to the raster stack:
  
  newproj = CRS('+proj=longlat +datum=WGS84')
  
  names(tminPptStack) = c('tmin','ppt')
  return(tminPptStack)
}

# Function to extract tmin, ppt to points and add the columns to the swd file:

extractClimateByDate <- function(rasterDirTmin, rasterDirPPt, dateValue, ptFile){
  # Prepare point file:
  pts <- ptFile %>%
    mutate(date = as.character(date)) %>%
    filter(date == dateValue) %>%
    as.data.frame
  spPts <- SpatialPointsDataFrame(coords = pts[,c('lon','lat')],
                                  data = pts,
                                  proj4string = CRS(projInfo))
  # Get raster file:
  rasterStack <- rasterPrepTminPpt(rasterDirTmin, rasterDirPPt, dateValue)
  # Extract To pts
  ptsEnv <- raster::extract(rasterStack, spPts)
  ptsEnv <- cbind(spPts@data, ptsEnv) %>%
    filter(!is.na(tmin))
  return(ptsEnv)
}

# Get unique vector of dates:

dateValues <- eBirdSamplingEnvSpring %>% .$date %>% unique 

# Make an empty list to store the data:

outList <- vector('list', length = length(dateValues))

# For loop to extract the rasters for a given date (may take a while to run):

for(i in 1:length(dateValues)){
  outList[[i]] <- extractClimateByDate(
    rasterDirTmin,
    rasterDirPPt,
    dateValues[i],
    eBirdSamplingEnvSpring
  )
}

# Make samples with data for spring, across sampling periods:

springSWDfull <- bind_rows(outList) %>%
  tbl_df %>%
  dplyr::mutate(
    month = lubridate::month(date),
    day = lubridate::day(date),
    year = lubridate::year(date),
    # Set sampling periods:
    samplingPeriod = ifelse(
      # Sampling period 1: March 1 - 15
      month == 3 & day <= 15, 1,
      ifelse(
        # Sampling period 2: March 16 - 29
        month == 3 & day <= 29, 2,
        # Sampling period 3: March 30 - April 12
        ifelse(
          month == 3|(month == 4 & day <= 12), 3,
               ifelse(
                 # Sampling period 4: April 13 - April 26
                 month == 4 & day <= 26, 4,
                 # Sampling period 5
                 ifelse(month == 4|month == 5 & day <= 10,
                 5, NA)))))
  ) %>%
  group_by(cellAddress, year, samplingPeriod) %>%
  mutate(count = max(count),
            tmin = mean(tmin),
            ppt = mean(ppt)) %>%
  ungroup %>%
  dplyr::select(cellAddress, count, year, 
                samplingPeriod, dev_hi:woodland, ppt, tmin) %>%
  distinct %>%
  dplyr::select(-cellAddress)

# Function to prepare samples with data for model running:

prepSWD <- function(swdIn, minFlockSize, maxFlockSize, samplingPeriodValue){
  swdIn %>%
    dplyr::filter(!(count > 0 & count < minFlockSize),
                  samplingPeriod == samplingPeriodValue) %>%
    mutate(
      pa = ifelse(count >= minFlockSize & count <= maxFlockSize,
                     1, 0)
    ) %>%
    dplyr::select(pa, dev_hi:tmin) %>%
    mutate(tmin2 = tmin^2)
}
  
  
  




