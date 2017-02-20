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

#-------------------------------------------------------------------------------*
# ---- Random swd's for testing ----
#-------------------------------------------------------------------------------*

# Function to generate a random sample of flock size 1 by flock size 2
# Note: Flock size 1 is the flocks size of interest, comparison to flock size 2

randomSWDpair <- function(swd1, swd2){
  nSamples <- nrow(dplyr::filter(swd1, pa == 1))
  # Get background points:
  bgPoints <- dplyr::filter(swd1, pa == 0)
  # Get random sample of presence points :
  bind_rows(
    dplyr::filter(swd1, pa == 1),
    dplyr::filter(swd2, pa == 1)
    ) %>%
    dplyr::sample_n(nSamples,
      replace = FALSE) %>%
    dplyr:: bind_rows(bgPoints)
}

# Make random prediction models. Because the result will be large, save as RDS
# rather than store in memory:

# Empty matrices for output

for(i in 1:10000){
  # Random samples of points:
  swdRandomLS <- randomSWDpair(swdLarge, swdSmall)
  swdRandomLM <- randomSWDpair(swdLarge, swdMedium)
  swdRandomMS <- randomSWDpair(swdMedium, swdSmall)
  # Get models of permuted samples:
  modLS <- runMaxentAreaPrevalence(swdRandomLS, bestModelLarge, 1)
  modLM <- runMaxentAreaPrevalence(swdRandomLM, bestModelLarge, 1)
  modMS <- runMaxentAreaPrevalence(swdRandomMS, bestModelMedium, 1.4)
  # Save random models:
  saveRDS(modLS, paste0(saveRDSlocation, 'modLS', i, '.rds'))
  saveRDS(modLM, paste0(saveRDSlocation, 'modLM', i, '.rds'))
  saveRDS(modMS, paste0(saveRDSlocation, 'modMS', i, '.rds'))
}

# Read in rds files, and fill matrices:

modCombos <- c('modLS', 'modLM','modMS')

outList <- vector('list', length = 3)

for(i in 1:length(outList)){
  rdsFilepaths <- list.files(RDSlocation, pattern = modCombos[i], full.names = TRUE)
  subOutList <- vector('list', length = length(rdsFilepaths))
  for(j in 1:length(rdsFilepaths)){
    subOutList[[j]] <- readRDS(rdsFilepaths[j]) %>%
      getAreaPrevalenceFrame(j) %>%
      data.frame %>%
      dplyr::mutate(permuted = modCombos[i])
  }
  outList[[i]] <- bind_rows(subOutList)
}

permutationFrame <- bind_rows(outList)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPalette <- c("#E69F00", "#009E73",
               "#0072B2", "#D55E00", "#CC79A7")

# Predicted area plots

predictedAreaLFplot <- ggplot(permutationFrame %>%
         mutate(
           permuted = str_replace(permuted, 'modLM', 'Large-Medium'),
           permuted = str_replace(permuted, 'modLS', 'Large-Small'),
           permuted = str_replace(permuted, 'modMS', 'Medium-Small'),
           Permutation = permuted
           ) %>%
         dplyr::filter(permuted != 'Medium-Small'),
       aes(x=area)) + 
  geom_density(aes(group=Permutation, fill=Permutation), alpha=0.7) +
  scale_fill_manual(values = cbPalette) +
  geom_segment(data = getAreaPrevalenceFrame(bestModelLarge,1),
               aes(x = area, xend = area, y = 0, yend = 35),
               size = 2, linetype = 1, lineend = 'round') +
  theme_bw() +
  labs(x = 'Fractional predicted area',
       y = 'Density') +
  theme(
    axis.title = element_text(size = 15)
  ) +
  xlim(.15, .4) +
  ylim(0, 38)


predictedAreaMFplot <- ggplot(permutationFrame %>%
                                mutate(
                                  permuted = str_replace(permuted, 'modLM', 'Medium-Large'),
                                  permuted = str_replace(permuted, 'modLS', 'Large-Small'),
                                  permuted = str_replace(permuted, 'modMS', 'Medium-Small'),
                                  Permutation = permuted
                                ) %>%
                                dplyr::filter(permuted != 'Large-Small'),
                              aes(x=area)) + 
  geom_density(aes(group=Permutation, fill=Permutation), alpha=0.7) +
  scale_fill_manual(values = cbPalette[-c(1:2)]) +
  geom_segment(data = getAreaPrevalenceFrame(bestModelMedium,1),
               aes(x = area, xend = area, y = 0, yend = 35),
               size = 1.5, linetype = 1, lineend = 'round') +
  theme_bw() +
  labs(x = 'Fractional predicted area',
       y = 'Density') +
  theme(
    axis.title = element_text(size = 15)
  ) +
  xlim(.15, .4) +
  ylim(0, 38)

# Predicted prevalence plots

predictedPrevalenceLFplot <- ggplot(permutationFrame %>%
                                mutate(
                                  permuted = str_replace(permuted, 'modLM', 'Large-Medium'),
                                  permuted = str_replace(permuted, 'modLS', 'Large-Small'),
                                  permuted = str_replace(permuted, 'modMS', 'Medium-Small'),
                                  Permutation = permuted
                                ) %>%
                                dplyr::filter(permuted != 'Medium-Small'),
                              aes(x=prevalence)) + 
  geom_density(aes(group=Permutation, fill=Permutation), alpha=0.7) +
  scale_fill_manual(values = cbPalette) +
  geom_segment(data = getAreaPrevalenceFrame(bestModelLarge,1),
               aes(x = prevalence, xend = prevalence, y = 0, yend = 35),
               size = 1.5, linetype = 1, lineend = 'round') +
  theme_bw() +
  labs(x = 'Prevalence',
       y = 'Density') +
  theme(
    axis.title = element_text(size = 15)
  ) +
  xlim(.15, .4) +
  ylim(0, 38)


predictedPrevalenceMFplot <- ggplot(permutationFrame %>%
                                mutate(
                                  permuted = str_replace(permuted, 'modLM', 'Medium-Large'),
                                  permuted = str_replace(permuted, 'modLS', 'Large-Small'),
                                  permuted = str_replace(permuted, 'modMS', 'Medium-Small'),
                                  Permutation = permuted
                                ) %>%
                                dplyr::filter(permuted != 'Large-Small'),
                              aes(x=prevalence)) + 
  geom_density(aes(group=Permutation, fill=Permutation), alpha=0.7) +
  scale_fill_manual(values = cbPalette[-c(1:2)]) +
  geom_segment(data = getAreaPrevalenceFrame(bestModelMedium,1),
               aes(x = prevalence, xend = prevalence, y = 0, yend = 35),
               size = 1.5, linetype = 1, lineend = 'round') +
  theme_bw() +
  labs(x = 'Prevalence',
       y = 'Density') +
  theme(
    axis.title = element_text(size = 15)
  ) +
  xlim(.15, .4) +
  ylim(0, 38)


png(filename = "C:/Users/Brian/Desktop/gits/blitzAnalysis/outPlots/areaPrevalence.png",
    width = 12, height = 4.5, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(
    predictedPrevalenceLFplot + guides(fill = FALSE) +
      geom_text(aes(x = .38, y = 28), label = 'Large\nflocks', size = 5) +
      geom_text(aes(x = .15, y = 35), label = 'A)', size = 5),
    predictedAreaLFplot  +
      geom_text(aes(x = .38, y = 28), label = 'Large\nflocks', size = 5) +
      geom_text(aes(x = .15, y = 35), label = 'B)', size = 5),
    predictedPrevalenceMFplot + guides(fill = FALSE) +
      geom_text(aes(x = .38, y = 28), label = 'Medium\nflocks', size = 5) +
      geom_text(aes(x = .15, y = 35), label = 'A)', size = 5),
    predictedAreaMFplot +
      geom_text(aes(x = .38, y = 28), label = 'Medium\nflocks', size = 5) +
      geom_text(aes(x = .15, y = 35), label = 'B)', size = 5), ncol=2)
    )
dev.off()


###################################################################

# Next: k-fold for each model for getting error in predictions

runMaxentK <- function(swd, bestModel, betaMultiplier, kFold){
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
    dplyr::select(-k)
  swdTest <- swdReduced %>%
    dplyr::filter(k == kFold & pa == 1) %>%
    bind_rows(swdAbsence) %>%
    dplyr::select(-k)
  # Set model arguments
  modArguments <- c('nothreshold', 'nohinge', 'noproduct','noquadratic',
                    str_c('betamultiplier=', betaMultiplier),
                    'addallsamplestobackground',
                    'writebackgroundpredictions',
                    'noautofeature','nooutputgrids',
                    'maximumiterations=10000', 'verbose')
  # Run maxent model with training and background data:
  maxentModel <- maxent(swdTrain[,-1], swdTrain[,1], args = modArguments)
  return(maxentModel)
}


outListSmall <- vector('list', length = 5)
outListMedium <- vector('list', length = 5)
outListLarge <- vector('list', length = 5)

for(i in 1:5){
  outListSmall[[i]] <- runMaxentK(swdSmall, bestModelSmall, 2.7, i)
  outListMedium[[i]] <- runMaxentK(swdMedium, bestModelMedium, 1.4, i)
  outListLarge[[i]] <- runMaxentK(swdLarge, bestModelLarge, 1, i)
}

getSummaryStatsK <- function(outList, flockSize){
  results <- outList@results %>% 
    data.frame
  results$stat <- row.names(results)
  row.names(results) <- NULL
  names(results) <- c('value', 'stat')
  resultsFrame <- results %>%
    dplyr::filter(str_detect(stat,'Prevalence')|
           str_detect(stat, 'Maximum.training.sensitivity.plus.')) %>%
    dplyr::filter(!str_detect(stat, '.cumulative'),
           !str_detect(stat, 'omission'))
  resultsFrame$stat <- c('prevalence', 'logThresh', 'area')
  resultsFrame$flockSize <- flockSize
  return(resultsFrame)
}

summaryListSmall <- vector('list', length = 5)
summaryListMedium <- vector('list', length = 5)
summaryListLarge <- vector('list', length = 5)

for(i in 1:5){
  summaryListSmall[[i]] <- getSummaryStatsK(outListSmall[[i]], 'small')
  summaryListMedium[[i]] <- getSummaryStatsK(outListMedium[[i]],  'medium')
  summaryListLarge[[i]] <- getSummaryStatsK(outListLarge[[i]], 'large')
}

summaryByFlock <- bind_rows(summaryListSmall, summaryListMedium, summaryListLarge) %>%
  group_by(flockSize, stat) %>%
  summarize(meanValue = mean(value),
            seValue = se(value),
            minCi = meanValue - seValue * 1.96,
            maxCi = meanValue + seValue * 1.96)




# Generate a random sample of flock size 1 by flock size 2
# Note: Flock size 1 is the flocks size of interest, comparison to flock size 2

random.swd.pair <- function(swd1, swd2){
  bgPoints <- dplyr::filter(swd1, pa == 0)
  swd1pres <- dplyr::filter(swd1, pa == 1)
  swd2pres <- dplyr::filter(swd2, pa == 1)
  # Remove unevaluated flock
  # Determine the sample size as the proportion of flock size 1:
  probFS1 = nrow(swd1pres)/(nrow(swd1pres) + nrow(swd2pres))
  # Bind flock size 1 and flock size 2 frames:
  flockData <- bind_rows(swd1pres, swd2pres)
  # Generate random value of 1 or 0 with the probability of obtaining
  # a 1 related to the number of s1 observations:
  flockData %>%
    mutate(s.rand = rbinom(nrow(.),1,probFS1)) %>%
    dplyr::filter(s.rand == 1) %>%
    dplyr::mutate(pa = 1) %>%
    dplyr::select(-s.rand) %>%
    dplyr:: bind_rows(bgPoints)
}

maxentRunRandomPair <- function(swd1, swd2, bestMod, betaMultiplier){
  swd <- random.swd.pair(swd1, swd2)
  # Get environmental variables to include from the best model:
  variablesToInclude <- getVariableContribution(bestModel) %>%
    .$variable
  # Remove environmental variables not used in this model:
  swdReduced <- swd %>%
    dplyr::select(pa) %>%
    bind_cols(
      swd %>% dplyr::select(one_of(variablesToInclude))
    ) 
  # Set model arguments
  modArguments <- c('nothreshold', 'nohinge', 'noproduct','noquadratic',
                    str_c('betamultiplier=', betaMultiplier),
                    'addallsamplestobackground',
                    'writebackgroundpredictions',
                    'noautofeature','nooutputgrids',
                    'maximumiterations=10000', 'verbose')
  # Run maxent model with training and background data:
  maxentModel <- maxent(swdReduced[,-1], swdReduced[,1], args = modArguments)
}

getSampledAIC <- function(swd1, swd2, bestMod, betaMultiplier){
  model <- maxentRunRandomPair(swd1, swd2, bestMod, betaMultiplier)
  as.numeric(model@results[5,1])
}

# Get sampled AUC for 100 random samples for large to medium, large to small,
# medium to small

sampleAUClistLM <- vector('list', length = 100)
sampleAUClistLS <- vector('list', length = 100)
sampleAUClistMS <- vector('list', length = 100)

for(i in 1:100){
  sampleAUClistLM[[i]] <- getSampledAIC(swdLarge, swdMedium, bestModelLarge, 2.7)
  sampleAUClistLS[[i]] <- getSampledAIC(swdLarge, swdSmall, bestModelLarge, 2.7)
  sampleAUClistMS[[i]] <- getSampledAIC(swdMedium, swdSmall, bestModelMedium, 1.4)
}


sampleAUClistSeb <- vector('list', length = 100)
sampleAUClistSblitz <- vector('list', length = 100)
sampleAUClistMeb <- vector('list', length = 100)
sampleAUClistMblitz <- vector('list', length = 100)
sampleAUClistLeb <- vector('list', length = 100)
sampleAUClistLblitz <- vector('list', length = 100)

for(i in 1:100){
  sampleAUClistSeb <- vector('list', length = 100)
  sampleAUClistSblitz <- vector('list', length = 100)
  sampleAUClistMeb <- vector('list', length = 100)
  sampleAUClistMblitz <- vector('list', length = 100)
  sampleAUClistLeb <- vector('list', length = 100)
  sampleAUClistLblitz <- vector('list', length = 100)
}



hist(unlist(sampleAUClistLM))

data.frame(I = unlist(sampleAUClistMS)) %>%
  tbl_df %>%
  ggplot(aes(I)) +
  geom_density(fill = 'gray90') +
  scale_x_continuous(limits = c(.75, 1)) +
  ylab('Density')+
  xlab ('Modified-Hellinger similarity (I)') +
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(1.5), vjust = .9),
        axis.title.x = element_text(size = rel(1.5), vjust = -0.4)) +
  # geom_vline(xintercept = I.lf.sf$I.actual, linewidth = 2.5, linetype = 'longdash')
  geom_segment(data = data.frame(I = bestModelMedium@results[5,1]),
               aes(x = I, y = 0, xend = I, yend = Inf))#,
# linewidth = 2.5, linetype = 'longdash')


