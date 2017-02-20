# maxent script for spring blitz
#===============================================================================*
# ---- SET-UP ----
#===============================================================================*

# Smart installer will check list of packages that are installed, install any
# necessary package that is missing, and load the library:

smartInstallAndLoad <- function(packageVector){
  for(i in 1:length(packageVector)){
    package <- packageVector[i]
    if(!package %in% rownames(installed.packages())){
      install.packages(packageVector[i],
                       repos="http://cran.rstudio.com/", 
                       dependencies=TRUE)
    }
  }
  lapply(packageVector, library, character.only = TRUE)
}

# Load and potentially install libraries:

smartInstallAndLoad(c('dplyr', 'tidyr','stringi','stringr', 'sp',
                      'lubridate','raster', 'dismo', 'ggplot2', 'rJava',
                      'grid', 'gridExtra', 'maps', 'maptools', 'rgeos'))

# Override some libraries for tidyverse functions:

filter <- dplyr::filter
select <- dplyr::select

# Set outPlots directory (example):

outPlotsDir <- 'C:/Users/Brian/Desktop/gits/blitzAnalysis/outPlots/'

# Because several steps in the process takes a while to run, I save the disk
# image at several pts:

diskImage <- 'spring_10-13.RData'

# Likewise, the location to save RDS files due to long run times:

saveRDSlocation <- 'C:/Users/Brian/Dropbox/randomRUBLmodels/'

# IMPORTANT! Prior to running the below, run the file "makeSWD.R" Also, in order
# to run this script, you need to place the maxent jar file in the folder: 
# win-library/[version]/dismo/java

#-------------------------------------------------------------------------------*
# ---- Basic functions----
#-------------------------------------------------------------------------------*

se <- function(x) {sd(x)/sqrt(length(x))}

conf95 <- function(x) se(x)*1.96 

cropRasterByMinLat <- function(rLayer, minLat){
  extentR <- extent(rLayer)
  extentR@ymin <- minLat
  crop(rLayer, extentR)
}

#-------------------------------------------------------------------------------*
# Get swd (across observation types):
#-------------------------------------------------------------------------------*

swdList <- vector('list', length = 5)

for(i in 1:length(swdList)){
  flockList <- vector('list', length = 3)
  names(flockList) <- c('Small', 'Medium', 'Large')
  flockList$Small <- prepSWD(swdSpring, 1,19, i)
  flockList$Medium <- prepSWD(swdSpring, 20,99, i)
  flockList$Large <- prepSWD(swdSpring, 100,Inf, i)
  swdList[[i]] <- flockList
}

#===============================================================================*
# ---- MODEL RUNNING AND CALIBRATION ----
#===============================================================================*

# Basic function to run a model:

maxentRun <- function(swd, betaMultiplier,
                      kFold = 'noCrossValidate', excludeVariables = NULL){
  swd <- swd %>%
    select(-cellAddress)
  if(!is.null(excludeVariables)){
    swd <- swd %>%
      select(-one_of(excludeVariables))
  }
  swdAbsence <- swd %>%
    filter(pa == 0) #%>%
  # Create input file of k fold:
  if(kFold != 'noCrossValidate'){
    swdTrain <- swd %>%
      filter(k != kFold & pa == 1) %>%
      bind_rows(swdAbsence) %>%
      select(-k) %>%
      data.frame
  } else {
    swdTrain <- swd %>%
      filter(pa == 1) %>%
      bind_rows(swdAbsence) %>%
      select(-k) %>%
      data.frame
  }
  # Set model arguments
  modArguments <- c('nothreshold','nohinge', 'noquadratic','noproduct',
                    str_c('betamultiplier=', betaMultiplier),
                    'addallsamplestobackground',
                    'writebackgroundpredictions',
                    'noautofeature','nooutputgrids',
                    'maximumiterations=10000', 'verbose')
  # Run maxent model with training and background data:
  maxentModel <- maxent(swdTrain[,-1], swdTrain[,1], args = modArguments)
  return(maxentModel)
}

# Function to run model twice, first with all variables, and then with reduced
# variables:

maxentRunReduced <- function(swd, betaMultiplier,
                             kFold = 'noCrossValidate', excludeVariables = NULL){
  modelFull <- maxentRun(swd, betaMultiplier)
  variablesToExclude <- getVariableContribution(modelFull) %>%
    filter(contribution <= 1 & variable != 'tmin') %>%
    .$variable
  modelReduced <- maxentRun(swd, betaMultiplier, excludeVariables = variablesToExclude)
  return(modelReduced)
}

# Function to get the contribution of environmental variables for a given model:

getVariableContribution <- function(model){
  modelResults <- model@results %>% data.frame
  modelResults$variable <- row.names(modelResults)
  names(modelResults) <- c('contribution','variable')
  variableContribution <- modelResults %>%
    mutate(includesContribution = str_detect(variable, 'contribution')) %>%
    filter(includesContribution == TRUE) %>%
    select(-includesContribution) %>%
    mutate(variable = str_replace_all(variable, '.contribution', ''))
  return(variableContribution)  
}

# Function to get a vector of variables to remove from a given model:

variablesToRemove <- function(model){
  getVariableContribution(model) %>%
  filter(contribution < 1 & variable != 'tmin') %>%
  .$variable
}

# Function to read the lambda file associated with a given model:

readLambdaFile <- function(model){
  lambdaData <- model@lambdas
  tf <- tempfile()
  writeLines(lambdaData, tf)
  read.csv(tf, fill = TRUE, header = F) %>%
    select(variable = V1, lambda = V2)
}

# Function to calculate AIC for a given swd file and beta multiplier:

calcAIC <- function(swd, betaMultiplier) {
  # Point data, all:
  swdPA <- swd %>%
    select(-c(cellAddress, pa))
  # Point data, presence:
  swdP <- swd %>%
    filter(pa == 1) %>%
    select(-c(cellAddress, pa))
  # Extract a  model:
  model <- maxentRunReduced(swd, betaMultiplier)
  # Extract lambdas file and convert to a data frame:
  lambdas <- model@lambdas
  lambda.df <- data.frame(do.call('rbind', strsplit(lambdas,',',fixed=TRUE)))
  lambda.df <- lambda.df[-((nrow(lambda.df)-4):(nrow(lambda.df))),]
  # Determing the number of parameters that had a lambda of zero:
  nRemoveParameters<- length(lambda.df[as.numeric(as.character(lambda.df[,2])) == 0,1])
  # Calculate the number of model parameters used in the model:
  nPar <- nrow(lambda.df)- nRemoveParameters
  presences <- swd %>% filter(pa == 1) %>% select(dev_hi:tmin2)
  absences <- swd %>% filter(pa == 0) %>% select(dev_hi:tmin2)
  # Model predictions across points:
  predictions <- predict(model, as.matrix(swd), args=c('outputformat=raw'))
  # Standardize predictinos such that they sum to 1:
  predStand <- predictions/sum(predictions)
  # Model predictions at presence points
  predictionsP <- predict(model, as.matrix(swdP), args=c('outputformat=raw'))
  # Get log likelihood:
  loglik <- sum(log(predictionsP), na.rm = TRUE)
  # How many points?
  n <- length(predictionsP)
  # Calculate the AICc
  AICc <- -2*loglik + 2*nPar + ((2*nPar*(nPar+1))/(n-nPar-1))
  # Output
  aicFrame <- data.frame(nParm = nPar, beta = betaMultiplier,AICc = AICc)
  return(aicFrame)
}

# Function to assess AICc values across beta values and output a table of
# results:

betaFinder <- function(swd, betaValues){
  out <- data.frame(matrix(nrow = length(betaValues),ncol = 3))
  for (i in 1:length(betaValues)){
    out[i,] <- calcAIC(swd, betaValues[i])
  }   
  colnames(out) <- c('nparm', 'beta','AICc')
  out
}

#-------------------------------------------------------------------------------*
# ---- Calibrate models ----
#-------------------------------------------------------------------------------*

# For a set of beta values from 1-15 and each beta period, get the best beta
# multiplier:

betaValues <- seq(1, 15, 1)

betaPeriods <- vector('list', length = 5)

for(i in 1:length(betaPeriods)){
  betaFlock <- vector('list', length = 3)
  names(betaFlock) <- c('Small', 'Medium', 'Large')
  betaFlock$Small <- betaFinder(swdList[[i]]$Small, betaValues)
  betaFlock$Medium <- betaFinder(swdList[[i]]$Medium, betaValues)
  betaFlock$Large <- betaFinder(swdList[[i]]$Large, betaValues)
  betaPeriods[[i]] <- betaFlock
}

# OUTPUT (because the above takes a while to run):

betaPeriods[[1]]$Large %>%
  arrange(AICc) %>%
  mutate(dAIC = AICc - min(AICc))

betaByPeriod <- vector('list', length = 5)

betaByPeriod[[1]] <- c(6, 3, 2)

betaByPeriod[[2]] <- c(4, 8, 5)

betaByPeriod[[3]] <- c(12, 6, 6)

betaByPeriod[[4]] <- c(11, 9, 3)

betaByPeriod[[5]] <- c(11, 10, 3)

# To look at AICc's, you can use the below:

betaPeriods[[1]]$Large %>%
  arrange(AICc) %>%
  mutate(dAIC = AICc - min(AICc))

# Save the disk image below, due to long run time:

save.image(diskImage)

#===============================================================================*
# ---- MODEL EVALUATION ----
#===============================================================================*

# Get best models:

bestModelPeriod <- vector('list', length = 5)

for(i in 1:length(bestModelPeriod)){
  swdPeriod <- swdList[[i]]
  betas <- betaByPeriod[[i]]
  bestModelFlock <- vector('list', length = 3)
  for(j in 1:length(bestModelFlock)){
    bestModelFlock[[j]] <- maxentRunReduced(swdPeriod[[j]], betaMultiplier = betas[j])
  }
  bestModelPeriod[[i]] <- bestModelFlock
}

# Save the disk image below, due to long run time:

save.image(diskImage)

#-------------------------------------------------------------------------------*
# ---- Get AUC values associated with test points ----
#-------------------------------------------------------------------------------*

# Function to run maxent model and return AUC for a given cross-validation fold:

runMaxentAUC <- function(swd, bestModel, betaMultiplier, kFold){
  swd <- swd %>%
    select(-cellAddress)
  # Get environmental variables to include from the best model:
  variablesToInclude <- getVariableContribution(bestModel) %>%
    .$variable
  # Remove environmental variables not used in this model:
  swdReduced <- swd %>%
    select(pa, k) %>%
    bind_cols(
      swd %>% select(one_of(variablesToInclude))
      )
  # Make background, training, and test points:
  swdAbsence <- swdReduced %>%
    filter(pa == 0)
  swdTrain <- swdReduced %>%
    filter(k != kFold & pa == 1) %>%
    bind_rows(swdAbsence) %>%
    select(-k) %>%
    data.frame
  swdTest <- swdReduced %>%
    filter(k == kFold & pa == 1) %>%
    select(-k) %>%
    data.frame
  # Set model arguments
  modArguments <- c('nothreshold','nohinge', 'noquadratic','noproduct',
                    str_c('betamultiplier=', betaMultiplier),
                    'addallsamplestobackground','writebackgroundpredictions',
                    'noautofeature','nooutputgrids',
                    'maximumiterations=10000', 'verbose')
  # Run maxent model with training and background data:
  maxentModel <- maxent(swdTrain[,-1], swdTrain[,1], args = modArguments)
  # Predict model values at test points:
  predictionPresence <- dismo::predict(maxentModel, as.matrix(swdTest))
  predictionAbsence <- dismo::predict(maxentModel,
                                      swdAbsence %>%
                                        select(-c(pa, k)) %>%
                                        as.matrix)
  # Evaluate model:
  evaluationObject <- dismo::evaluate(p = predictionPresence,
                               a = predictionAbsence)
  # Return data frame with auc and cor values:
  data.frame(auc = evaluationObject@auc, 
             cor = evaluationObject@cor,
             row.names = NULL)
}

# Function to get summary stats from auc and cor:

getAUC <- function(swd, bestModel, betaMultiplier){
  # For loop runs maxent evaluation for each value of K
  outMat <- matrix(nrow = 5, ncol = 2)
  for(i in 1:5){
    outMat[i,] <- runMaxentAUC(swd, bestModel, betaMultiplier, i) %>%
      as.matrix
  }
  colnames(outMat) <- c('auc', 'cor')
  outFrame <- data.frame(outMat) %>%
    summarize(
       meanAUC = mean(auc),
       seAUC = se(auc),
       meanCor = mean(cor),
       seCor = se(cor)
    )
  outFrame
}

# Get AUC summary stats across observation periods and flock sizes:

aucPeriod <- vector('list', length = 5)

for(i in 1:length(aucPeriod)){
  swdPeriod <- swdList[[i]]
  betas <- betaByPeriod[[i]]
  bestModels <- bestModelPeriod[[i]]
  flockSizes <- c('Small', 'Medium', 'Large')
  aucFlock <- vector('list', length = 3)
  for(j in 1:length(bestModelFlock)){
    aucFlock[[j]] <- getAUC(
      swdPeriod[[j]], bestModels[[j]], betaMultiplier = betas[j]) %>%
      mutate(flockSize = flockSizes[j],
             period = i)
  }
  aucPeriod[[i]] <- bind_rows(aucFlock)
}

# Bind summaries into a single data frame:

summaryAUC <- bind_rows(aucPeriod)

# Save the disk image below, due to long run time:

save.image(diskImage)

#-------------------------------------------------------------------------------*
# ---- Predict models to raster ----
#-------------------------------------------------------------------------------*
# Example = 2014, see plotting suitability rasters in the plotting section
# Set periods:

periodValues <- 1:5

rList2014 <- vector('list', length = 5)

for(i in 1:length(periodValues)){
  rList2014[[i]] <- stack(
    rStack, 
    summarizeClimate(rasterDirTmin, rasterDirPPt, periodValues[i], 2014)
  )
}

# For a given best model and year make logistic prediction:

getLogisticPrediction <- function(bestModel, rasterStack){
  predict(bestModel, rasterStack,
          args='outputformat=logistic', 
          progress='text')
}

bestModelPredictions <- vector('list', length = 5)

for(i in 1:5){
  envR <- rList2014[[i]]
  bestModelPredictionsFlock <- vector('list', length = 3)
  for(j in 1:3){
    bestModelPredictionsFlock[[j]] <- getLogisticPrediction(
      bestModelPeriod[[i]][[j]],
      envR
    )
  }
  bestModelPredictions[[i]] <- bestModelPredictionsFlock
}

bestModelPredictions[[5]][[3]] <- NULL

# Save the disk image below, due to long run time:

save.image(diskImage)

# To plot the above, go to script plotSuitabilityMaps.R

#-------------------------------------------------------------------------------*
# ---- Area, prevalence, threshold ----
#-------------------------------------------------------------------------------*

# Function to get the predicted area and logistic threshold values for a given
# model:

getAreaPrevalenceFrame <- function(model){
  results <- model@results %>%
    data.frame
  results$stat <- row.names(results)
  row.names(results) <- NULL
  names(results) <- c('value', 'stat')
  results %>%
    filter(
      str_detect(stat,'Prevalence')|
        str_detect(stat, 'Equal.training.sensitivity.and.specificity.logistic.threshold')|
        str_detect(stat, 'Equal.training.sensitivity.and.specificity.area')
      ) %>%
    mutate(
      stat = ifelse(str_detect(stat, 'Prevalence'), 'prevalence', stat),
      stat = ifelse(str_detect(stat, 'logistic'), 'threshold', stat),
      stat = ifelse(str_detect(stat, 'area'), 'area', stat)
      ) %>%
    spread(stat, value)
}

# Add the flock names for the best models by period:

for(i in 1:length(bestModelPeriod)){
  names(bestModelPeriod[[i]]) <- c('Small', 'Medium', 'Large')
}

# Calculate predicted area for each period and flock size class:

areaPrevalenceByPeriodL <- vector('list', length = 5)
flockNames <- c('Small', 'Medium', 'Large')

for(i in 1:5){
  areaPrevalenceByFlock <- vector('list', length = 3)
  for(j in 1:3){
    model <- bestModelPeriod[[i]][[j]]
    areaPrevalenceByFlock[[j]] <- getAreaPrevalenceFrame(model) %>%
      mutate(flockSize = flockNames[j],
             samplingPeriod = i)
  }
  areaPrevalenceByPeriodL[[i]] <- bind_rows(areaPrevalenceByFlock)
}

# Bind rows into a single data frame:

areaPrevalenceSummary <- bind_rows(areaPrevalenceByPeriodL)

# Save the disk image below, due to long run time:

save.image(diskImage)

#-------------------------------------------------------------------------------*
# ---- Variable contribution ----
#-------------------------------------------------------------------------------*

# Get variable contribution frome for each flock size class and sampling period:

flockNames <- c('Small', 'Medium', 'Large')

envVariables <- c('tmin', 'flood', 'rowcrop', 'ppt', 
                  'flood', 'pasture', 'dev_hi','weth', 'wetw')

variableContributionByPeriod <- vector('list', length = 5)

for(i in 1:5){
  variableContributionByFlock <- vector('list', length = 3)
  for(j in 1:3){
    model <- bestModelPeriod[[i]][[j]]
    variableContributionByFlock[[j]] <- getVariableContribution(model) %>%
      mutate(flockSize = flockNames[j],
             samplingPeriod = i)
  }
  variableContributionByPeriod[[i]] <- bind_rows(variableContributionByFlock)
}

variableContributionSummary <- bind_rows(variableContributionByPeriod)

# Turn summary into a plotting-friendly data frame:

dataStackedBar <- bind_rows(variableContributionByPeriod) %>%
  tbl_df %>%
  mutate(variable = ifelse(variable == 'tmin2', 'tmin', variable)) %>%
  # Make the tmin contribution the sum of tmin and tmin2:
  group_by(flockSize, samplingPeriod, variable) %>%
  summarize(contribution = sum(contribution)) %>%
  ungroup %>%
  # Set some variables to "other":
  mutate(variable = as.character(variable),
         variable = ifelse(
           !(variable %in% c('tmin', 'flood', 'rowcrop', 'ppt', 
                             'flood', 'pasture', 'dev_hi','weth', 'wetw')),
           'other', variable)
  ) %>%
  # Sum the other percentages for each flock size class:
  group_by(samplingPeriod, flockSize, variable) %>%
  summarize(varContribution = sum(contribution)) %>%
  ungroup %>%
  # Reset the factor levels and labels for plotting:
  mutate(
    variable = variable %>%
      factor(
        levels = c('flood','pasture','rowcrop','ppt','dev_hi','tmin','weth','wetw','other'),
        labels = c('Floodplain (+)    ', 'Pasture (+/-)   ',
                   'Row crop (+)    ','Precipitation (+)    ', 'Highly developed (-)    ',
                   'Minimum temperature (+/-)    ', 'Emergent wetland  (+)     ',
                   'Woody wetland (+)     ','Other    '),
        ordered = TRUE),
    samplingPeriod = samplingPeriod %>%
      factor(
        levels = c(1:5),
        labels = c("March 1-14", "March 15-28","March 29-April 11",
                   "April 12-25","April 26-May 9")
      )
  )

#-------------------------------------------------------------------------------*
#---- Kolmogorov-Smirnov tests ---- 
#-------------------------------------------------------------------------------*
# A simple test to observe whether distribution of samples about an environemntal
#variable is distinguishable between two flock size classes

# A sampling frame for K-S analysis and plotting (below):

springSamples <- swdSpring %>%
  filter(count > 0) %>%
  mutate(flockSize = ifelse(
    count < 20,'Small',
    ifelse(count < 100, 'Medium',
           'Large')
  ),
  samplingPeriod = samplingPeriod %>%
    factor(
      levels = c(1:5),
      labels = c("March 1-14", "March 15-28","March 29-April 11",
                 "April 12-25","April 26-May 9")
    ))

# Function to run K-S tests for a given environmental variable and time period:

ksFun <- function(envVar, periodValue){
  # Set small-medium, small-large, and medium-large comparisons:
  compareFrame <- data.frame(
    sp1 = c('Small', 'Small', 'Medium'),
    sp2 = c('Medium', 'Large', 'Large')
  )
  # For each flock size class comparison, get stats:
  for(i in 1:nrow(compareFrame)){
    env1 <- springSamples %>%
      filter(flockSize == compareFrame[i,1],
             samplingPeriod == periodValue) %>%
      select_(envVar) %>%
      as.matrix %>% as.vector
    env2 <- springSamples %>%
      filter(flockSize == compareFrame[i,2],
             samplingPeriod == periodValue) %>%
      select_(envVar) %>%
      as.matrix %>% as.vector
    testOut <- ks.test(env1, env2)
    compareFrame$D[i] = testOut$statistic
    compareFrame$p[i] = testOut$p.value
    compareFrame$meanSp1 <- mean(env1)
    compareFrame$meanSp2 <- mean(env2)
  }
  return(compareFrame)
}

# Function to run K-S tests across time periods for a given environmental variable:

ksListPeriodFun <- function(envVar){
  ksListPeriod <- vector('list', length = 5)
  for(i in 1:5){
    ksListPeriod[[i]] <- ksFun(envVar, i) %>%
      mutate(samplingPeriod = i,
             p = round(p, 3),
             envVariable = envVar)
  }
  bind_rows(ksListPeriod)
}

# Run KS tests across environmental variables and periods:

ksList <- vector('list', length = length(envVariable))

for(j in 1:length(envVariables)){
  ksList[[j]] <- ksListPeriodFun(envVariables[j])
}

# Make list into a data frame:

ksFrame <- bind_rows(ksList)

# To look at output for a given variable and time period, you use the following
# example as a guide:

ksFrame %>%
  filter(samplingPeriod == 1,
         envVariable == 'tmin')

#===============================================================================*
#---- PLOTTING ---- 
#===============================================================================*
# Plots include:
#   1. Suitability maps, 2014 example (due to temperature and precip rasters)
#   2. AUC
#   3. Area/prevalence
#   4. Variable contribution plots
#   5. Temperature distribution plots

#-------------------------------------------------------------------------------*
# ---- Plot suitability maps ----
#-------------------------------------------------------------------------------*

# Function to crop polygons by bounding box:

gClip <- function(shp, bb){
  if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else b_poly <- as(extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}

# Get US states polygons:

states <-maps::map("state", col="transparent", plot=FALSE, fill = TRUE)

# Clip US states polygons to the raster output extent:

us <- map2SpatialPolygons(states, IDs = states$names,
                          proj4string=CRS("+proj=longlat +datum=WGS84")) %>%
  gSimplify(tol = 0.0001) %>%
  gBuffer(byid = TRUE, width = 0) %>%
  gClip(bbox(extent(rStack)))

# Aggregate then mask rasters to shapefiles and make into a data frame:

prepRasterToPlot = function(r, flockSize){
  r %>%
    # aggregate(2) %>%
    mask(us) %>%
    rasterToPoints %>%
    data.frame %>%
    mutate(flockClass = flockSize)
}

# Prepare rasters for plotting in ggplot for each period:

rListPeriod <- vector('list', length = 5)

for(i in 1:5){
  rListPeriod[[i]] <- bind_rows(
    prepRasterToPlot(bestModelPredictions[[i]][[1]], 'Small'),
    prepRasterToPlot(bestModelPredictions[[i]][[2]], 'Medium'),
    prepRasterToPlot(bestModelPredictions[[i]][[3]], 'Large')
  ) %>%
    mutate(samplingPeriod = i) 
}

# List to data frame and provide plot-friendly labels:

rastersForPlotting <- bind_rows(rListPeriod) %>%
  mutate(
    flockClass = factor(flockClass,
                        levels = c('Small', 'Medium', 'Large')),
    samplingPeriod = factor(samplingPeriod,
                            labels = c("March 1-14", "March 15-28",
                                       "March 29-April 11",
                                       "April 12 - 25",
                                       "April 26 - May 9"))
  )

# Plot suitability maps:

rastersForPlotting %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = layer)) +
  facet_grid(samplingPeriod ~ flockClass) +
  scale_fill_gradientn(name = 'Habitat\nsuitability',
                       colours=c('navyblue', 'darkblue', 'blue' , 'yellow', 'red', 'darkred'), 
                       na.value = 'gray70') +
  geom_map(map = fortify(us), data = us, aes(map_id = id, x = long, y = lat, group = group), 
           fill = NA, color = 'black', size = .15) +
  labs(x = 'Longitude', y = 'Latitude') +
  theme_bw() + 
  theme(axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(.6)),
        axis.ticks.length = unit(1, 'mm'),
        strip.text = element_text(size = rel(0.6)),
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6)))

# Save plot output:

ggsave(paste0(outPlotsDir, 'springSuitabilityMaps2014.png'), 
       width = 8.5, height = 11, units = 'in')
dev.off()


#-------------------------------------------------------------------------------*
# ---- Plot AUC by observation type ----
#-------------------------------------------------------------------------------*

# Color blind palette:

cbPalette <- c("#E69F00", "#009E73", "#D55E00")

# Plot AUC:

aucPlot <- ggplot(
  summaryAUC %>%
    filter(!(flockSize == 'Large' & period == 5)),
  aes(x = period, y = meanAUC, color = flockSize)) +
  geom_errorbar(aes(ymin=meanAUC - 1.96*seAUC,
                    ymax=meanAUC + 1.96*seAUC),
                position=position_dodge(width=0.1), width = 0,
                size = .75) +
  geom_point(size = 4, position=position_dodge(width=0.1)) +
  theme_bw() +
  ylab('AUC') +
  xlab('Observation period') +
  scale_color_manual(values = cbPalette, guide = guide_legend(title = 'Flock size')) +
  geom_line(size = 1, position=position_dodge(width=0.1)) +
  theme(axis.title = element_text(margin=margin(0,10,0,0))) +
  theme(legend.key = element_blank()) +
  theme(axis.title = element_text(size = 15))

# Save plot output:

png(filename = paste0(outPlotsDir, 'aucSpring.png'),
    width = 9, height = 4.5, units = 'in', res = 300)
aucPlot + 
  scale_x_continuous(
    breaks = 1:5,
    labels = c("March 1-14", "March 15-28","March 29-April 11",
               "April 12-25","April 26-May 9")
  )
dev.off()

#-------------------------------------------------------------------------------*
# ---- Plot area and prevalence ----
#-------------------------------------------------------------------------------*

# Color blind palette:

cbPalette <- c("#E69F00", "#009E73", "#D55E00")

# Plot area:

areaPlot <- ggplot(
  areaPrevalenceSummary %>%
    filter(!(flockSize == 'Large' & samplingPeriod == 5)),
  aes(x = samplingPeriod, y = area, color = flockSize)) +
  geom_point(size = 4, position=position_dodge(width=0.1)) +
  theme_bw() +
  ylab('Fractional Predicted Area') +
  xlab('Observation period') +
  ylim(.18, .41) +
  scale_color_manual(values = cbPalette, guide = guide_legend(title = 'Flock size')) +
  geom_line(size = 1, position=position_dodge(width=0.1)) +
  theme(axis.title = element_text(margin=margin(0,10,0,0))) +
  theme(legend.key = element_blank()) +
  theme(axis.title = element_text(size = 15)) + 
  scale_x_continuous(
    breaks = 1:5,
    labels = c("Mar. 1-14", "Mar. 15-28","Mar. 29-Apr. 11",
               "Apr. 12 - 25","Apr. 26 - May 9")
  )

# Plot prevalence:

prevalencePlot <- ggplot(
  areaPrevalenceSummary %>%
    filter(!(flockSize == 'Large' & samplingPeriod == 5)),
  aes(x = samplingPeriod, y = prevalence, color = flockSize)) +
  geom_point(size = 4, position=position_dodge(width=0.1)) +
  theme_bw() +
  ylab('Prevalence') +
  xlab('Observation period') +
  ylim(.18, .41) +
  scale_color_manual(values = cbPalette, guide = guide_legend(title = 'Flock size')) +
  geom_line(size = 1, position=position_dodge(width=0.1)) +
  theme(axis.title = element_text(margin=margin(0,10,0,0))) +
  theme(legend.key = element_blank()) +
  theme(axis.title = element_text(size = 15)) + 
  scale_x_continuous(
    breaks = 1:5,
    labels = c("Mar. 1-14", "Mar. 15-28","Mar. 29-Apr. 11",
               "Apr. 12 - 25","Apr. 26 - May 9")
  )

# Arrange as grid and save plot output:

png(filename = paste0(outPlotsDir, 'areaPrevalenceSpring.png'),
    width = 11, height = 4.5, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(
    prevalencePlot + guides(color = FALSE) +
      geom_text(aes(x = 4.5, y = .4), label = 'A)', size = 5, color = 1),
    areaPlot  +
      geom_text(aes(x = 4.5, y = .4), label = 'B)', size = 5, color = 1)
    , ncol=2, widths = c(.95,1.1))
)
dev.off()

#-------------------------------------------------------------------------------*
# ---- Plot variable contribution ----
#-------------------------------------------------------------------------------*

# Set plot theme:

plot_theme <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 7),
          axis.title.x = element_text(size = 12, vjust = -0.4),
          axis.title.y = element_text(size = 12, vjust = 1.1),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 7),
          panel.border = element_blank())
}

# Colors for bars:

colorBlindPalette <- c("#E69F00",  "#009E73", "#F0E442", 
                       "#0072B2", "#D55E00","#56B4E9",
                       '#9C661F','#68228B',
                       "gray60")

# Plots including temperature:

plotWithTemp <- ggplot(
  dataStackedBar %>%
    arrange(variable), 
  aes(x = flockSize, y  = varContribution, 
      fill = variable, order = as.numeric(variable))) + 
  geom_bar(stat = 'identity') + 
  geom_bar(stat = 'identity',color = 'black',
           show.legend = FALSE, size = .5) +
  plot_theme() +
  labs(x = 'Flock size class', y = 'Variable contribution (%)') +
  scale_fill_manual(name = 'Variable', values = colorBlindPalette) +
  theme(legend.title = element_blank(),
        legend.key.height = unit(1,"line"),
        legend.key = element_rect(size=2, color = 'white'),
        legend.position = 'top') +
  facet_wrap(~samplingPeriod, nrow = 1)

colorBlindPalette2 <- c("#E69F00",  "#009E73", "#F0E442", 
                        "#0072B2", "#D55E00",'#9C661F','#68228B',"gray60")

# Plots excluding temperature:

plotNoTemp <- ggplot(
  dataStackedBar %>% 
    arrange(variable) %>%
    filter(!str_detect(variable, 'temperature')) %>%
    group_by(flockSize) %>% 
    mutate(TotalContribution = sum(varContribution)) %>% 
    ungroup %>% 
    mutate(relativeContribution = varContribution/TotalContribution), 
  aes(x = flockSize, y  = varContribution, 
      fill = variable, order = as.numeric(variable))) + 
  geom_bar(stat = 'identity') + 
  geom_bar(stat = 'identity',color = 'black', show.legend = FALSE, size = .5) +
  plot_theme() +
  labs(x = 'Flock size class', y = 'Variable contribution (%)') +
  scale_fill_manual(name = 'Variable', values = colorBlindPalette2) +
  theme(legend.key.height = unit(2,"line")) +
  facet_wrap(~samplingPeriod, nrow = 1)

# Make plots into a list:

plots <- list(plotWithTemp, plotNoTemp)

# Set plot grid and legend:

g <- ggplotGrob(plots[[1]] + theme(legend.position="top"))$grobs

legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

# Save plot output:

png(filename = paste0(outPlotsDir, 'VariableContributionSpring.png'), 
    width = 8.5, height = 5.5, units = 'in', res = 300)
gridExtra::grid.arrange(legend, 
                        gridExtra::arrangeGrob(
                          plotWithTemp +
                            theme(
                              legend.position = 'none',
                              axis.title.x = element_blank()
                            ),
                          plotNoTemp + theme(legend.position = 'none'),
                          nrow = 2, heights = c(4.5,6)),
                        nrow = 2, heights = c(1, 6))
dev.off()


#-------------------------------------------------------------------------------*
# ---- Plot temperature distributions ----
#-------------------------------------------------------------------------------*

tempPlot <- springSamples %>%
  dplyr::select(samplingPeriod, flockSize, tmin) %>%
  ggplot(aes(x = tmin)) + 
  facet_grid(flockSize~samplingPeriod) +
  geom_density(aes(y = ..scaled..,fill = 'darkgray',  alpha = 0.8)) + 
  scale_fill_manual(values = cbPallete) + theme_bw() + 
  ylab('Scaled density') +
  xlab(expression(Average~minimum~temperature~(degree ~ C))) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15)
  ) + theme(legend.position="none")


png(filename = paste0(outPlotsDir, 'tminSpring.png'), 
    width = 8.5, height = 5.5, units = 'in', res = 300)
tempPlot
dev.off()

#### END ####