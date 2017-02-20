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
                      'gridExtra'))

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
    dplyr::select(-cellAddress)
  if(!is.null(excludeVariables)){
    swd <- swd %>%
      dplyr::select(-one_of(excludeVariables))
  }
  swdAbsence <- swd %>%
    dplyr::filter(pa == 0) #%>%
  # Create input file of k fold:
  if(kFold != 'noCrossValidate'){
    swdTrain <- swd %>%
      dplyr::filter(k != kFold & pa == 1) %>%
      bind_rows(swdAbsence) %>%
      dplyr::select(-k) %>%
      data.frame
  } else {
    swdTrain <- swd %>%
      dplyr::filter(pa == 1) %>%
      bind_rows(swdAbsence) %>%
      dplyr::select(-k) %>%
      data.frame
  }
  # Set model arguments
  modArguments <- c('nothreshold','nohinge', 'noquadratic','noproduct',
                    str_c('betamultiplier=', betaMultiplier),
                    'addallsamplestobackground','writebackgroundpredictions',
                    'noautofeature','nooutputgrids',
                    'maximumiterations=10000', 'verbose')
  # Run maxent model with training and background data:
  maxentModel <- maxent(swdTrain[,-1], swdTrain[,1], args = modArguments)
  return(maxentModel)
}

# Function to run model twice, first with all variables, and then with reducted variables:

maxentRunReduced <- function(swd, betaMultiplier,
                             kFold = 'noCrossValidate', excludeVariables = NULL){
  modelFull <- maxentRun(swd, betaMultiplier)
  variablesToExclude <- getVariableContribution(modelFull) %>%
    dplyr::filter(contribution <= 1 & variable != 'tmin') %>%
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
    dplyr::filter(includesContribution == TRUE) %>%
    dplyr::select(-includesContribution) %>%
    mutate(variable = str_replace_all(variable, '.contribution', ''))
  return(variableContribution)  
}

# Function to get a vector of variables to remove from a given model:

variablesToRemove <- function(model){
  getVariableContribution(model) %>%
  dplyr::filter(contribution < 1 & variable != 'tmin') %>%
  .$variable
}

# Function to read the lambda file associated with a given model

readLambdaFile <- function(model){
  lambdaData <- model@lambdas
  tf <- tempfile()
  writeLines(lambdaData, tf)
  read.csv(tf, fill = TRUE, header = F) %>%
    dplyr::select(variable = V1, lambda = V2)
}

# Function to calculate AIC for a given swd file and beta multiplier:

calcAIC <- function(swd, betaMultiplier) {
  # Point data, all:
  swdPA <- swd %>%
    dplyr::select(-c(cellAddress, pa))
  # Point data, presence:
  swdP <- swd %>%
    filter(pa == 1) %>%
    dplyr::select(-c(cellAddress, pa))
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
  presences <- swd %>% dplyr::filter(pa == 1) %>% dplyr::select(dev_hi:tmin2)
  absences <- swd %>% dplyr::filter(pa == 0) %>% dplyr::select(dev_hi:tmin2)
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

# Save the disk image below, due to long run time:

save.image(diskImage)

# OUTPUT (because the above takes a while to run):
# Also note, there are some funky values in the above (due to sampling a subset, looking for overall pattern in AIC)

betaPeriods[[1]]$Large %>%
  arrange(AICc) %>%
  mutate(dAIC = AICc - min(AICc))

betaByPeriod <- vector('list', length = 5)

betaByPeriod[[1]] <- c(6, 3, 2)

betaByPeriod[[2]] <- c(4, 8, 5)

betaByPeriod[[3]] <- c(12, 6, 6)

betaByPeriod[[4]] <- c(11, 9, 3)

betaByPeriod[[5]] <- c(11, 10, 3)

# Save the disk image below, due to long run time:

save.image(diskImage)

# To look at AICc's, you can use the below:

betaPeriods[[1]]$Large %>%
  arrange(AICc) %>%
  mutate(dAIC = AICc - min(AICc))

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

# I save the disk image below, due to long run time:

save.image(diskImage)

#-------------------------------------------------------------------------------*
# ---- Get AUC values associated with test points ----
#-------------------------------------------------------------------------------*

# Function to run maxent model and return AUC for a given cross-validation fold:

runMaxentAUC <- function(swd, bestModel, betaMultiplier, kFold){
  swd <- swd %>%
    dplyr::select(-cellAddress)
  # Get environmental variables to include from the best model:
  variablesToInclude <- getVariableContribution(bestModel) %>%
    .$variable
  # Remove environmental variables not used in this model:
  swdReduced <- swd %>%
    dplyr::select(pa, k) %>%
    bind_cols(
      swd %>% dplyr::select(one_of(variablesToInclude))
      )
  # Make background, training, and test points:
  swdAbsence <- swdReduced %>%
    dplyr::filter(pa == 0)
  swdTrain <- swdReduced %>%
    dplyr::filter(k != kFold & pa == 1) %>%
    bind_rows(swdAbsence) %>%
    dplyr::select(-k) %>%
    data.frame
  swdTest <- swdReduced %>%
    dplyr::filter(k == kFold & pa == 1) %>%
    dplyr::select(-k) %>%
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
                                        dplyr::select(-c(pa, k)) %>%
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

# I save the disk image below, due to long run time:

save.image(diskImage)

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

# Save plot to outplots directory:

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
# ---- Predict models to raster ----
#-------------------------------------------------------------------------------*
# Example = 2014
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

# Take a cursory look at output:

for(i in 1:5) plot(bestModelPredictions[[i]][[2]])

save.image(diskImage)

# To plot these, go to script plotSuitabilityMaps.R

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
    dplyr::filter(
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

# Arrange as grid and save:

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