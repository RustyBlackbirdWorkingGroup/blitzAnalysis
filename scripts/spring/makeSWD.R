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

# Override some libraries for tidyverse functions:

filter <- dplyr::filter
select <- dplyr::select

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

# Set path to land cover data:

# pathToRasterData <- 'C:/Users/Brian/Dropbox/rubl_12_15/'  # Helm
# pathToRasterData <- '/Users/bsevans/Dropbox/rubl_12_15/'   # MacBook Air
# pathToRasterData <- 'C:/Users/Default.Default-THINK/Dropbox/rubl_12_15/' # Thinkpad

# Get land cover data:

rStack <- loadEnv(paste0(pathToRasterData, 'lc_asc'))

# Remove minimum temperature and precipitation:

rStack <- rStack[[-c(8,11)]]

# Raster data to be used for addresses:

r <- rStack[[1]]

# Get projection info from raster stack:

projInfo <- raster::projection(r)

#***************************************************************************************************
#***************************************************************************************************

#---------------------------------------------------------------------------------------------------*
# ----  FUNCTIONS TO DOWNLOAD CLIMATE DATA, IF NECESSARY ----
#---------------------------------------------------------------------------------------------------*

# Function to download the minimum temperature raster for a given date:

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

# Function to download the precipitation raster a given data:

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

#---------------------------------------------------------------------------------------------------*
# ----  FUNCTIONS TO GENERATE PPT AND TMIN RASTERS FOR A GIVEN SAMPLING PERIOD ----
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

# Get a data frame of dates associated with each period and year:

periodStarts <- c('03-01', '03-15', '03-29', '04-12', '04-26')
years <- 2014:2016

outList <- vector('list', length = 3)
for(i in 1:length(years)){
  periodStartDates <- paste(years[i], periodStarts, sep = '-')
  periodList <- vector('list', length = 5)
  for(j in 1:length(periodStarts)){
    periodStart <- as.Date(periodStartDates[j])
    periodList[[j]] <- data.frame(
      samplingPeriod = j,
      date = c(periodStart, periodStart + days(1:13))
    )
  }
  outList[[i]] <- bind_rows(periodList) %>%
    mutate(year = years[i])
}

periodDates <- bind_rows(outList)
  
# Function to get summary climate data for a given period and year:

summarizeClimate <- function(rasterDirTmin, rasterDirPPt, periodValue, yearValue){
  dateValues <- periodDates %>%
    filter(year == yearValue,
           samplingPeriod == periodValue) %>%
    .$date
  # Get list of tmin and ppt rasters for each date:
  tminList <- vector('list', length = length(dateValues))
  pptList <- vector('list', length = length(dateValues))
  for(i in 1:length(dateValues)){
    tminList[[i]] <- raster(paste0(rasterDirTmin, dateValues[i])) %>%
      crop(r)
    pptList[[i]] <-  raster(paste0(rasterDirPPt, dateValues[i])) %>%
      crop(r)
  }
  # Calculate the minimum (mintemp) and summed precipitation across the period:
  tminR <- min(stack(tminList))
  pptR <- sum(stack(pptList))
  # Return stack, with tmin2
  outStack <- stack(list(ppt = pptR, tmin = tminR, tmin2 = tminR^2))
  return(outStack)
}

#===================================================================================================*
# ---- PREPARE BIRD OBSERVATIONS ----
#===================================================================================================*

# Function to match date with period:

periodDateMatch <- function(periodDateValue){
  periodDates %>%
    filter(dates == periodDateValue) %>%
    .$period
}

#---------------------------------------------------------------------------------------------------*
# ---- PREPARE SAMPLING DATA ----
#---------------------------------------------------------------------------------------------------*

# Rusty observations:

rustyListsSpring <- read.csv('rublEbird.csv') %>%
  tbl_df %>%
  #filter to study extent:
   filter(lon > extent(r)[1] & lon < extent(r)[2],
         lat > extent(r)[3] & lat < extent(r)[4]) %>%
  # Subset to dates associated with the spring blitz:
  mutate(date = as.Date(date)) %>%
  filter(date %in% periodDates$date) %>%
  # Remove counts recorded as 'X':
  filter(count != 'X') %>%
  mutate(count = as.numeric(count)) %>%
  select(observationID, count)

# Make a vector of rusty lists where count is recorded as 'X':

rustyXobservations <- read.csv('rublEbird.csv') %>%
 filter(count == 'X') %>%
  .$observationID

# Get eBird list data:

eBirdListsSpring <- read.csv(pathToEbirdListData) %>%
  select(observationID, lat, lon, date) %>%
  tbl_df %>%
  #filter to the study extent and study dates:
  filter(lon > extent(r)[1] & lon < extent(r)[2],
                lat > extent(r)[3] & lat < extent(r)[4]) %>% #,
  mutate(date = as.Date(date)) %>%
  filter(as.Date(date) %in% periodDates$date) %>%
  mutate(year = lubridate::year(date)) %>%
  # Remove observationIDs where rusty count was reported as X:
  filter(!observationID %in% rustyXobservations) %>%
  # Add count data:
  left_join(rustyListsSpring, by = 'observationID') %>%
  # Change na counts (no observation match) to 0:
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  # Add period data:
  inner_join(periodDates %>%
               select(-year), by = 'date')

# Add cell addresses:

eBirdListsSpring$cellAddress <- cellFromXY(
  r,
  eBirdListsSpring %>%
    select(lon, lat) %>%
    data.frame %>%
    SpatialPoints(proj4string = CRS(projInfo)) 
)

# Get max count for each cell address, year, and sampling period:

springSampling <- eBirdListsSpring %>%
  group_by(cellAddress, year, samplingPeriod) %>%
  summarize(count = max(count)) %>%
  ungroup

#---------------------------------------------------------------------------------------------------*
# ---- ADD LAND COVER AND CLIMATE DATA (MAKE SWD) ----
#---------------------------------------------------------------------------------------------------*

# Extract raster land cover data by cell ID:

envByCell <- data.frame(
  cellAddress = (springSampling$cellAddress %>% unique),
  raster::extract(
    x = rStack,
    y = (springSampling$cellAddress %>% unique),
    df = TRUE)
) %>%
  tbl_df %>%
  select(-ID)

# Join sampling data to environment data and add year field:

samplingEnvSpring <- left_join(
  springSampling,
  envByCell,
  by = 'cellAddress') %>%
  filter(!is.na(dev_hi))

# Function to extract precip and tmin data for a given cell:

samplingWithClimateFun <- function(rasterDirTmin, rasterDirPPt, periodValue, yearValue){
  # Filter sampling frame to period and year:
  samplingSubset <- samplingEnvSpring %>%
    filter(year == yearValue & samplingPeriod == periodValue)
  # Get raster summary data:
  climateStack <- summarizeClimate(rasterDirTmin, rasterDirPPt, periodValue, yearValue)
  # Extract to cell addresses:
  climateByCell <- data.frame(
    cellAddress = (samplingSubset$cellAddress %>% unique),
    raster::extract(
      x = climateStack,
      y = (samplingSubset$cellAddress %>% unique),
      df = TRUE)) %>%
    tbl_df %>%
    select(-ID)
  # Join to sampling subset associated with period and year value:
  samplingSubsetEnv <- left_join(
    samplingSubset,
    climateByCell,
    by = 'cellAddress')
  return(samplingSubsetEnv)
}

# For loop to add climate data to points:

periods <- 1:5
years <- 2014:2016

outListPeriod <- vector('list', length = 5)

for(i in 1:5){
  outListYear <- vector('list', length = 3)
  for(j in 1:length(years)){
    outListYear[[j]] <- samplingWithClimateFun(rasterDirTmin, rasterDirPPt, periods[i], years[j]) 
  }
  outListPeriod[[i]] <- bind_rows(outListYear)
}

# Make samples with data:

swdSpring <- bind_rows(outListPeriod) %>%
  select(-year) %>%
  filter(!is.na(tmin))

#---------------------------------------------------------------------------------------------------*
# ---- PREPARE SAMPLES ----
#---------------------------------------------------------------------------------------------------*

# Function to prepare samples with data for model running:

prepSWD <- function(swdIn, minFlockSize, maxFlockSize, samplingPeriodValue){
  swdSamplingPeriod <- swdIn %>%
    filter(samplingPeriod == samplingPeriodValue)
  bind_rows(
    swdSamplingPeriod %>%
      filter(count == 0) %>%
      mutate(pa = 0,
             k = NA),
    swdSamplingPeriod %>%
      filter(count >= minFlockSize & count <= maxFlockSize) %>%
      mutate(pa = 1,
             k = dismo::kfold(x = pa, k = 5, by = pa))
  ) %>%
    select(cellAddress, pa, dev_hi:tmin2, k)
}
