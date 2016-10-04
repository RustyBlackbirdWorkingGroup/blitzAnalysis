# Summarize point data by raster cell


# note: 
# setwd('C:/Users/Brian/Desktop/gits/RUBL/rubl_winter/') # Helm
# setwd('C:/Users/Default.Default-THINK/Desktop/gits/RUBL/rubl_winter') # Thinkpad

library(sp)
library(raster)
library(dplyr)
library(stringr)
library(dismo)

options(stringsAsFactors = FALSE)

# pathToRasterData <- 'C:/Users/Brian/Dropbox/rubl_12_15/'  # Helm
# pathToRasterData <- '/Users/bsevans/Dropbox/rubl_12_15/'   # MacBook Air
# pathToRasterData <- 'C:/Users/Default.Default-THINK/Dropbox/rubl_12_15/' # Thinkpad


# pathToEbirdListData <- 'C:/Users/Brian/Dropbox/eBirdListData.csv'   # Helm
# pathToEbirdListData <- 'C:/Users/Default.Default-THINK/Dropbox/eBirdListData.csv' # ThinkPad
# pathToEbirdListData <- '/Users/bsevans/Dropbox/eBirdListData.csv'   # MacBook Air

#---------------------------------------------------------------------------------------------------*
# ---- FUNCTIONS ----
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

# Load tmin raster for a given date:

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

# Load tmin raster for a given date:

downloadTminRasterMonthly <- function(year){
  month = c('01', '02')
  outList <- vector('list', length = 2)
  for(i in 1:2){
    ftp <- paste0('ftp://prism.oregonstate.edu/monthly/tmin/',
                  year, '/',
                  'PRISM_tmin_stable_4kmM2_',
                  paste0(year, month[i]),
                  '_bil.zip')
    tmin <- tempfile()
    download.file(ftp, tmin)
    # Unzip:
    tmin <- unzip(tmin)
    # Convert to raster format:
    tmin <- raster(tmin[1])
    # Provide projection information:
    crs(tmin) <- "+proj=longlat +datum=WGS84"
    outList[[i]] <- crop(tmin, rStack[[1]])
  }
  writeRaster(overlay(stack(outList), fun = mean), 
              paste0('tmin', year),
              overwrite = TRUE)
  }


downloadPPTRasterMonthly <- function(year){
  month = c('01', '02')
  outList <- vector('list', length = 2)
  for(i in 1:2){
    ftp <- paste0('ftp://prism.oregonstate.edu/monthly/ppt/',
                  year, '/',
                  'PRISM_ppt_stable_4kmM3_',
                  paste0(year, month[i]),
                  '_bil.zip')
    ppt <- tempfile()
    download.file(ftp, ppt)
    # Unzip:
    ppt <- unzip(ppt)
    # Convert to raster format:
    ppt <- raster(ppt[1])
    # Provide projection information:
    crs(ppt) <- "+proj=longlat +datum=WGS84"
    outList[[i]] <- crop(ppt, rStack[[1]])
  }
  writeRaster(overlay(stack(outList), fun = mean), 
              paste0('ppt', year),
              overwrite = TRUE)
}
# 
# years <- 2006:2016
# 
# for(i in 1:length(years)){
#   downloadTminRasterMonthly(years[i])
#   downloadPPTRasterMonthly(years[i])
# }

#---------------------------------------------------------------------------------------------------*
# ---- GET CELL DATA AND SUMMARIZE SAMPLING EFFORT TO DATE AND CELL ----
#---------------------------------------------------------------------------------------------------*

# Get raster data:

rStack <- loadEnv(paste0(pathToRasterData, 'lc_asc'))

# Raster data to be used for addresses:

r <- rStack[[1]]

projInfo = raster::projection(r)

# Get rusty blackbird data:

rustyLists <- read.csv('rublEbird.csv') %>%
  tbl_df %>%
  # Subset to dates associated with the winter blitz:
  mutate(date = as.Date(date)) %>%
  mutate(month = lubridate::month(date),
         day = lubridate::day(date)) %>%
  filter(month == 1|
           (month == 2 & day <15),
         lubridate::year(date) < 2015) %>%
  dplyr::select(-c(month, day)) %>%
  # Remove counts recorded as 'X':
  filter(count != 'X') %>%
  mutate(count = as.numeric(count))
  
# Make a vector of rusty lists where count is recorded as 'X':

rustyXobservations <- read.csv('rublEbird.csv') %>%
  filter(count == 'X') %>%
  .$observationID
  
# Get eBird list data:

eBirdLists <- read.csv(pathToEbirdListData) %>%
  tbl_df %>%
  # Filter to the study extent:
  filter(lon > extent(r)[1] & lon < extent(r)[2],
         lat > extent(r)[3] & lat < extent(r)[4]) %>%
  # Subset lists to only dates associated with the rusty observations:
  filter(date %in%
           (rustyLists$date %>% 
              unique %>% 
              as.character)) %>%
  # Remove observationIDs where rusty count was reported as X:
  filter(!observationID %in% rustyXobservations) %>%
  # Add count data:
  left_join(rustyLists %>%
              dplyr::select(observationID, count),
            by = 'observationID') %>%
  # Change na counts (no observation match) to 0:
  mutate(count = ifelse(is.na(count), 0, count))

# Add cell addresses

eBirdLists$cellAddress <- cellFromXY(
  r,
  eBirdLists %>%
    dplyr::select(lon, lat) %>%
    data.frame %>%
    SpatialPoints(proj4string = CRS(projInfo)) 
)

# Extract raster land cover data by cell ID:

envByCell <- data.frame(
  cellAddress = (eBirdLists$cellAddress %>% unique),
  raster::extract(
  x = rStack,
  y = (eBirdLists$cellAddress %>% unique),
  df = TRUE)
) %>%
  tbl_df %>%
  dplyr::select(-c(ID, tmin, ppt))

# Join sampling data to environment data and add year field:

eBirdSamplingEnv <- left_join(
  eBirdLists,
  envByCell,
  by = 'cellAddress'
) %>%
  mutate(date = as.Date(date),
         year = lubridate::year(date))

# Add precipitation and temperature data:

years <- 2006:2014

samplingByYearList <- vector('list', length = length(years))

for(i in 1:length(years)){
  samplingSubset <- eBirdSamplingEnv %>%
    dplyr::filter(year == years[i]) %>%
    dplyr::select(cellAddress, year, protocol, count, dev_hi:woodland) %>%
    group_by(cellAddress, protocol) %>%
    dplyr::mutate(count = max(count)) %>%
    ungroup %>%
    distinct
  # Get precipitation and min temperature rasters:
  pptR <- raster(paste0(pathToRasterData, 'climateRasters/ppt',years[i]))
  tminR <- raster(paste0(pathToRasterData, 'climateRasters/tmin',years[i]))
  # Extract precip and temperature to points (by cell address):
  samplingSubset$ppt <- raster::extract(pptR, samplingSubset$cellAddress)
  samplingSubset$tmin <- raster::extract(tminR, samplingSubset$cellAddress)
  # Output list item, removing NA tmin and ppt (outside of extent)
  samplingByYearList[[i]] <- samplingSubset %>%
    dplyr::filter(!is.na(tmin), !is.na(ppt)) #%>%
}

samplingAcrossYears <- bind_rows(samplingByYearList)

#---------------------------------------------------------------------------------------------------*
# ---- ATTACH RUBL SAMPLES FOR A GIVEN CELL AND DATE ----
#---------------------------------------------------------------------------------------------------*

prepSWD <- function(minFlockSize,maxFlockSize, years, protocolChoice = 'all'){
  swdBG <- samplingAcrossYears %>%
    group_by(cellAddress, year) %>%
    dplyr::mutate(count = max(count)) %>%
    ungroup %>%
    dplyr::select(cellAddress, year, count, dev_hi:tmin) %>%
    dplyr::filter(count < 1, year %in% years) %>%
    distinct %>%
    dplyr::mutate(pa = 0) %>%
    dplyr::select(pa, dev_hi:tmin)
  if(protocolChoice == 'all'){
    paFrame <- samplingAcrossYears %>%
      dplyr::filter(year %in% years) %>%
      group_by(cellAddress, year) %>%
      dplyr::mutate(count = max(count)) %>%
      ungroup %>%
      dplyr::filter(count >= minFlockSize & count <= maxFlockSize) %>%
      mutate(pa = 1) %>%
      dplyr::select(pa, dev_hi:tmin) %>%
      bind_rows(swdBG) %>%
      mutate(tmin2 = tmin^2) %>%
      distinct
  }
  if(protocolChoice == 'eb') {
    paFrame <- samplingAcrossYears %>%
      dplyr::filter(year %in% years,
                    str_detect(protocol, 'Traveling Count')) %>%
      group_by(cellAddress, year) %>%
      dplyr::mutate(count = max(count)) %>%
      ungroup %>%
      dplyr::filter(count >= minFlockSize & count <= maxFlockSize) %>%
      mutate(pa = 1) %>%
      dplyr::select(pa, dev_hi:tmin) %>%
      bind_rows(swdBG) %>%
      mutate(tmin2 = tmin^2) %>%
      distinct
  }
  if(protocolChoice == 'blitz'){
    paFrame <- samplingAcrossYears %>%
      dplyr::filter(year %in% years,
                    str_detect(protocol, 'Blitz')) %>%
      group_by(cellAddress, year) %>%
      dplyr::mutate(count = max(count)) %>%
      ungroup %>%
      dplyr::filter(count >= minFlockSize & count <= maxFlockSize) %>%
      mutate(pa = 1) %>%
      dplyr::select(pa, dev_hi:tmin) %>%
      bind_rows(swdBG) %>%
      mutate(tmin2 = tmin^2)
  }
  # Add k to paFrame and return output
  return(
    paFrame %>%
      mutate(k = dismo::kfold(x = pa, k = 5, by = pa)) %>%
      data.frame
  )
}



