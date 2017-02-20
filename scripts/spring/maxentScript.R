# maxent script for spring blitz
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
                      'lubridate','raster', 'dismo', 'ggplot2', 'rJava'))

# Set outPlots directory (example):

outPlotsDir <- 'C:/Users/Brian/Desktop/gits/blitzAnalysis/outPlots/'

# Because several steps in the process takes a while to run, I save the disk image at several pts:

diskImage <- 'spring_10-13.RData'

# IMPORTANT! Prior to running the below, run the file "makeSWD.R"
# Also, in order to run this script, you need to place the maxent jar file in the folder: win-library/[version]/dismo/java

#----------------------------------------------------------------------------*
# Basic functions:
#----------------------------------------------------------------------------*

se <- function(x) {sd(x)/sqrt(length(x))}

conf95 <- function(x) se(x)*1.96 

cropRasterByMinLat <- function(rLayer, minLat){
  extentR <- extent(rLayer)
  extentR@ymin <- minLat
  crop(rLayer, extentR)
}

#----------------------------------------------------------------------------*
# Get swd (across observation types):
#----------------------------------------------------------------------------*

swdList <- vector('list', length = 5)

for(i in 1:length(swdList)){
  flockList <- vector('list', length = 3)
  names(flockList) <- c('Small', 'Medium', 'Large')
  flockList$Small <- prepSWD(swdSpring, 1,19, i)
  flockList$Medium <- prepSWD(swdSpring, 20,99, i)
  flockList$Large <- prepSWD(swdSpring, 100,Inf, i)
  swdList[[i]] <- flockList
}

#===================================================================================================*
# ---- MODEL RUNNING AND CALIBRATION ----
#===================================================================================================*

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

# Get the contribution of environmental variables for a given model:

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

# Get a vector of variables to remove from a given model:

variablesToRemove <- function(model){
  getVariableContribution(model) %>%
  dplyr::filter(contribution < 1 & variable != 'tmin') %>%
  .$variable
}


# Read the lambda file associated with a given model

readLambdaFile <- function(model){
  lambdaData <- model@lambdas
  tf <- tempfile()
  writeLines(lambdaData, tf)
  read.csv(tf, fill = TRUE, header = F) %>%
    dplyr::select(variable = V1, lambda = V2)
}

#----------------------------------------------------------------------------*
# Calculate AIC:
#----------------------------------------------------------------------------*

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

#----------------------------------------------------------------------------*
# Assess AICc values across beta values and output a table of results:
#----------------------------------------------------------------------------*

betaFinder <- function(swd, betaValues){
  out <- data.frame(matrix(nrow = length(betaValues),ncol = 3))
  for (i in 1:length(betaValues)){
    out[i,] <- calcAIC(swd, betaValues[i])
  }   
  colnames(out) <- c('nparm', 'beta','AICc')
  out
}

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

# I save the disk image below, due to long run time:

save.image(diskImage)

# To look at AICc's, you can use the below:

betaPeriods[[1]]$Large %>%
  arrange(AICc) %>%
  mutate(dAIC = AICc - min(AICc))

#===================================================================================================*
# ---- MODEL EVALUATION ----
#===================================================================================================*

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

#----------------------------------------------------------------------------*
# ---- Get AUC values associated with test points ----
#----------------------------------------------------------------------------*

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

#----------------------------------------------------------------------------*
# ---- Plot AUC by observation type ----
#----------------------------------------------------------------------------*

cbPallete <- c("#E69F00", "#009E73", "#D55E00")

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
  scale_color_manual(values = cbPallete, guide = guide_legend(title = 'Flock size')) +
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

#----------------------------------------------------------------------------*
# ---- Make model predictions (raster maps) ----
#----------------------------------------------------------------------------*

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

#----------------------------------------------------------------------------*
# ---- Area, prevalence, threshold ----
#----------------------------------------------------------------------------*

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

#----------------------------------------------------------------------------*
# ---- Plot area and prevalence ----
#----------------------------------------------------------------------------*

cbPallete <- c("#E69F00", "#009E73", "#D55E00")

areaPlot <- ggplot(
  areaPrevalenceSummary %>%
    filter(!(flockSize == 'Large' & samplingPeriod == 5)),
  aes(x = samplingPeriod, y = area, color = flockSize)) +
  geom_point(size = 4, position=position_dodge(width=0.1)) +
  theme_bw() +
  ylab('Fractional Predicted Area') +
  xlab('Observation period') +
  ylim(.18, .41) +
  scale_color_manual(values = cbPallete, guide = guide_legend(title = 'Flock size')) +
  geom_line(size = 1, position=position_dodge(width=0.1)) +
  theme(axis.title = element_text(margin=margin(0,10,0,0))) +
  theme(legend.key = element_blank()) +
  theme(axis.title = element_text(size = 15)) + 
  scale_x_continuous(
    breaks = 1:5,
    labels = c("Mar. 1-14", "Mar. 15-28","Mar. 29-Apr. 11",
               "Apr. 12 - 25","Apr. 26 - May 9")
  )

prevalencePlot <- ggplot(
  areaPrevalenceSummary %>%
    filter(!(flockSize == 'Large' & samplingPeriod == 5)),
  aes(x = samplingPeriod, y = prevalence, color = flockSize)) +
  geom_point(size = 4, position=position_dodge(width=0.1)) +
  theme_bw() +
  ylab('Prevalence') +
  xlab('Observation period') +
  ylim(.18, .41) +
  scale_color_manual(values = cbPallete, guide = guide_legend(title = 'Flock size')) +
  geom_line(size = 1, position=position_dodge(width=0.1)) +
  theme(axis.title = element_text(margin=margin(0,10,0,0))) +
  theme(legend.key = element_blank()) +
  theme(axis.title = element_text(size = 15)) + 
  scale_x_continuous(
    breaks = 1:5,
    labels = c("Mar. 1-14", "Mar. 15-28","Mar. 29-Apr. 11",
               "Apr. 12 - 25","Apr. 26 - May 9")
  )

grid.arrange(
  arrangeGrob(
    areaPlot + guides(color = FALSE) +
      geom_text(aes(x = 4.5, y = .4), label = 'A)', size = 5, color = 1),
    prevalencePlot  +
      geom_text(aes(x = 4.5, y = .4), label = 'B)', size = 5, color = 1)
    , ncol=2, widths = c(.95,1.1))
)

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

# Make random prediction models. Because the result will be large, save as RDS rather than store in memory:

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
  saveLoc <- 'C:/Users/Brian/Dropbox/randomRUBLmodels/'
  saveRDS(modLS, paste0(saveLoc, 'modLS', i, '.rds'))
  saveRDS(modLM, paste0(saveLoc, 'modLM', i, '.rds'))
  saveRDS(modMS, paste0(saveLoc, 'modMS', i, '.rds'))
}

# Read in rds files, and fill matrices:

# RDSlocation <- 'C:/Users/Brian/Dropbox/randomRUBLmodels'


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

cbPallete <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPallete <- c("#E69F00", "#009E73",
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
  scale_fill_manual(values = cbPallete) +
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
  scale_fill_manual(values = cbPallete[-c(1:2)]) +
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
  scale_fill_manual(values = cbPallete) +
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
  scale_fill_manual(values = cbPallete[-c(1:2)]) +
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


  
# library(gridExtra)
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








#############################################################################################
  



# 



# Run maxent model with training and background data:
maxentModel <- maxent(swdTrain[,-1], swdTrain[,1], args = modArguments)


# evaluationObject <- dismo::evaluate(p = swdTest %>%
#                                dplyr::filter(pa > 0) %>%
#                                .$predictions,
#                              a = swdTest %>%
#                                dplyr::filter(pa == 0) %>%
#                                .$predictions)

modelOutList <- function(minFlockSize, maxFlockSize, years, 
                         betaMultiplier, protocolChoice = 'all'){
  swd <- prepSWD(minFlockSize,maxFlockSize, years, protocolChoice)
  outList <- vector('list', 5)
  for (i in 1:5){
    outList[[i]] <- maxentRun(swd, betaMultiplier, i)
  }
  outList
}




maxentRun <- function(swd, betaMultiplier,
                      kFold = 'noCrossValidate', excludeVariables = NULL)



# Get auc value to evaluate models of each flock size class:

evaluationListS <- vector('list', length =5)
evaluationListM <- vector('list', length =5)
evaluationListL <- vector('list', length =5)


for(i in 1:5){
  evaluationListS[[i]] <- maxentCrossValidation(1, 20, 2009:2011, 
                                                i, betaMultiplier, protocolChoice = 'all')@auc
  evaluationListM[[i]] <- maxentCrossValidation(21, 99, 2009:2011, 
                                                i, betaMultiplier, protocolChoice = 'all')@auc
  evaluationListL[[i]] <- maxentCrossValidation(100, Inf, 2009:2011, 
                                    i, betaMultiplier, protocolChoice = 'all')@auc
}

aucFrame <- data.frame(
  flockSize = c('S', 'M', 'L'),
  meanAUC = c(
    mean(unlist(evaluationListS)),
    mean(unlist(evaluationListM)),
    mean(unlist(evaluationListL))),
  confidence = c(
    sd(unlist(evaluationListS))/sqrt(5)*1.96,
    sd(unlist(evaluationListM))/sqrt(5)*1.96,
    sd(unlist(evaluationListL))/sqrt(5)*1.96)
  )


# Run model across k


maxentRun <- function(minFlockSize, maxFlockSize, years, betaMultiplier, protocolChoice = 'all'){
  # Create input file of k fold::
  swd <- prepSWD(minFlockSize,maxFlockSize, years, protocolChoice) %>%
    data.frame %>%
    dplyr::select(-k)
  # Set model arguments
  betaR <- str_c('betamultiplier=', betaMultiplier)
  modArguments <- c('nothreshold', 'nohinge', #'noproduct','noquadratic',
                    betaR, 'addallsamplestobackground',
                    'writebackgroundpredictions',#'replicates=5','writeplotdata',
                    'noautofeature','nooutputgrids',
                    'maximumiterations=10000', 'verbose')
  # Run maxent model with training and background data:
  maxentModel <- maxent(swd[,-1], swd[,1], args = modArguments)
  return(maxentModel)
}

smallFlockModel <- maxentRun(1, 19, 2009:2011, 1)

mediumFlockModel <- maxentRun(20, 99, 2009:2011, 0)

largeFlockModel <- maxentRun(99, Inf, 2009:2011, 0)

# 

minFlockSize = 90
maxFlockSize = Inf
years = 2009:2011
kFold = 3
betaR = 2
protocolChoice = 'all'

test14 <- maxentCrossValidation(minFlockSize, maxFlockSize, years, kFold, betaR, protocolChoice)

modelEvaluate <- function(minFlockSize, maxFlockSize, years, protocolChoice = 'all'){
  presenceLC <- prepSWD(minFlockSize,maxFlockSize, years, protocolChoice = 'all') %>%
    dplyr::filter(pa == 1) %>%
    dplyr::select(dev_hi:tmin) %>%
    as.matrix
  absenceLC <- prepSWD(minFlockSize,maxFlockSize, years, protocolChoice = 'all') %>%
    dplyr::filter(pa == 0) %>%
    dplyr::select(dev_hi:tmin) %>%
    as.matrix
  evaluate(presenceLC, absenceLC, test14[[2]])
}

slotNames(test14[[2]]$models)



modelEvaluate(minFlockSize, maxFlockSize, years, protocolChoice = 'all')

testEvaluate <- evaluate(
  p = ,
  a = prepSWD(minFlockSize,maxFlockSize, years, protocolChoice = 'all') %>%
    dplyr::filter(pa == 0) %>%
    dplyr::select(dev_hi:tmin) %>%
    as.matrix,
  model = test14[[2]]
)

testPredict <- predict(test14$maxentModel, 
        prepSWD(minFlockSize,maxFlockSize, years, protocolChoice = 'all') %>%
          dplyr::filter(k == kFold) %>%
          data.frame %>%
          dplyr::select(dev_hi:tmin)
)
          

test15<- maxentRun(swd$y2015, flockSizeClass, kFold, beta.multiplier)
test15

test16 <- maxentRun(swd$y2016, flockSizeClass, kFold, beta.multiplier)
test16
# 
# mTest <- maxentEvaluate(observationClass, flockSizeClass, kFold, beta.multiplier)

maxentEvaluate <- function(observationClass, flockSizeClass, kFold, beta.multiplier){
  # Run model
  maxentModelOut <- maxentRun(observationClass, flockSizeClass, kFold, beta.multiplier)
  maxentModel <- maxentModelOut$maxentModel
  # Create test point files for evaluation
  swdTest <-   swd[[observationClass]][[flockSizeClass]] %>%
    dplyr::filter(k == kFold)
  testPresence <- swdTest %>%
    dplyr::filter(sp == 1) %>%
    dplyr::select(-c(sp,lon, lat,k)) # dplyr::select(lon, lat)%>% 
  #SpatialPoints(proj4string = CRS(projection(rStack2009)))
  testAbsence <- swdTest %>%
    dplyr::filter(sp == 0) %>%
    dplyr::select(-c(sp,lon, lat,k)) # dplyr::select(lon, lat)%>%%>%
  #SpatialPoints(proj4string = CRS(projection(rStack2009)))
  modelEvaluation <- evaluate(testPresence, testAbsence, maxentModel)
  return(list(maxentModel = maxentModel, modelEvaluation = modelEvaluation,
              swd = maxentModelOut$swd, swdTest = swdTest))
}


# Function to run across folds:

maxentAcrossFolds <- function(observationClass, flockSizeClass, beta.multiplier){
  kFold <- 1:5
  outList <- list(length = 5)
  aucList <- list(length = 5)
  corList <- list(length = 5)
  presenceList <- list(length = 5)
  absenceList <- list(length = 5)
  
  for(i in 1:5){
    outList[[i]]<- maxentEvaluate(observationClass, flockSizeClass, kFold[i],
                                  beta.multiplier)
    aucList[[i]] <- outList[[i]]$modelEvaluation@auc
    corList[[i]] <- outList[[i]]$modelEvaluation@cor
    presenceList[[i]] <- data.frame(pa = 1, 
                                    predicted = outList[[i]]$modelEvaluation@presence)
    absenceList[[i]] <- data.frame(pa = 0, 
                                   predicted = outList[[i]]$modelEvaluation@absence)
  }
  aucValues <- as.numeric(unlist(aucList))
  corValues <- as.numeric(unlist(corList))
  paFrame <- rbind(do.call('rbind', presenceList),
                   do.call('rbind', absenceList))
  return(list(modelOut = outList, aucValues = aucValues, corValues = corValues, paFrame = paFrame))
}

# Across data (all), flock sizes:

observationClass = 'all'
beta.multiplier = 0

allIndOut <- maxentAcrossFolds(observationClass, 'ind', beta.multiplier)
allSfOut <- maxentAcrossFolds(observationClass, 'sf', beta.multiplier)
allLfOut <- maxentAcrossFolds(observationClass, 'lf', beta.multiplier)

# Functions to extract variable contribution and lambda data:

se <- function(x) sd(x)/sqrt(length(x))

readLambdaFile <- function(model, i){
  lambdaData <- model[[1]][[i]][[1]]@lambdas[1:15]
  tf <- tempfile()
  writeLines(lambdaData, tf)
  read.csv(tf, fill = TRUE, header = F) %>%
    dplyr::select(variable = V1, lambda = V2)
}

makeLambdaContributionFrame <- function(model){
  dfContributionMat <- matrix(nrow = 15, ncol = 5)
  lambdaMat <- matrix(nrow = 15, ncol = 5)
  
  for(i in 1:ncol(lambdaMat)){
    lambdaMat[,i] <- readLambdaFile(model, i)[,2]
    dfContributionMat[,i] <- model[[1]][[i]][[1]]@results[7:21]
  }
  
  lambdaSummary <- data.frame(
    variable = as.character(readLambdaFile(model,1)$variable), 
    meanLambda = apply(lambdaMat, 1, mean),
    seLambda = apply(lambdaMat, 1, se))
  
  variableContributionSummary <- data.frame(
    variable = as.character(readLambdaFile(model,1)$variable), 
    meanContribution = apply(dfContributionMat, 1, mean),
    seContribution = apply(dfContributionMat, 1, se))
  
  outFrame <- plyr::join(variableContributionSummary, lambdaSummary) %>%
    arrange(desc(meanContribution))
  
  return(outFrame)
}

lambdaContributionFrame_allIndOut <- makeLambdaContributionFrame(allIndOut)

lambdaContributionFrame_allSfOut <- makeLambdaContributionFrame(allSfOut)

lambdaContributionFrame_allLfOut <- makeLambdaContributionFrame(allLfOut)
# 
# 
# df2 <- data.frame(variable = df1[,1], 
#                   percentContribution = apply(df1[,-1], 1, mean),
#                   standardError = apply(df1[,-1], 1, se)) %>%
#                   arrange(desc(percentContribution))
# 


prob.r.stack = function(model, outformat){
  r = stack()
  for (i in 1:5){
    r1 = predict(model[[1]][[i]][[1]],rStack2009,
                 args=c(paste('outputformat=',outformat, sep = '')), 
                 progress='text')
    r = stack(r, r1)
  }
  r
}

allIndOutRlogistic <- prob.r.stack(allIndOut, 'logistic')
allSfOutRlogistic <- prob.r.stack(allSfOut, 'logistic')
allLfOutRlogistic <- prob.r.stack(allLfOut, 'logistic')

indLogisticMean <- mean(allIndOutRlogistic)
sfLogisticMean <- mean(allSfOutRlogistic)
lfLogisticMean <- mean(allLfOutRlogistic)

#-------------------------------------------------------------------------------*
# Niche equivalency:
#-------------------------------------------------------------------------------*

flockData <- swdRUBL[[1]] 

# Generate a random sample of species 1 by species 2
# Note: Species 1 is the species of interest, comparison to species 2

random.swd.pair = function(sp1, sp2){
  # Remove unevaluated flock
  flockData <- flockData %>%
    dplyr::filter(sp == sp1|sp == sp2)
  nSp1 <- nrow(dplyr::filter(flockData, sp == sp1))
  # Determine the sample size as the proportion of species 1:
  prob.s1 = nSp1/nrow(flockData)
  # Generate random value of 1 or 0 with the probability of obtaining
  # a 1 related to the # of s1 observations:
  flockData %>%
    mutate(s.rand = rbinom(nrow(.),1,prob.s1)) %>%
    dplyr::filter(s.rand == 1) %>%
    dplyr::mutate(sp = 1) %>%
    dplyr::select(-s.rand)
}

# Function to run maxent for a given dataset:

maxentRunRawPlot = function(inFlockData, beta.multiplier = 0){
  # Model input:
  swdFlock <- inFlockData
  max.in <- rbind(swdFlock, swdBG) %>%
    dplyr::select(-c(lat, lon))
  # Set model arguments
  beta.r <- str_c('betamultiplier=', beta.multiplier)
  mod.args <- c('nothreshold', 'nohinge', 'noproduct',#'noquadratic',
                beta.r, 'addallsamplestobackground',
                'writebackgroundpredictions','writeplotdata',
                'noautofeature','nooutputgrids',
                'maximumiterations=10000', 'verbose')
  # Run maxent model with training and background data:
  maxentModel <- maxent(max.in[,-1], max.in[,1], args = mod.args)
  return(predict(maxentModel, rStack2009, args = c('outputformat=raw')))
}

# Run models of empirical data:

lf <- maxentRunRawPlot(dplyr::filter(flockData, sp == 'lf') %>%
                         mutate(sp = 1))
sf <- maxentRunRawPlot(dplyr::filter(flockData, sp == 'sf') %>%
                         mutate(sp = 1))
ind <- maxentRunRawPlot(dplyr::filter(flockData, sp == 'ind') %>%
                          mutate(sp = 1))

# Run null models for flock size pairs
# Note: Each of the null surface lists are 382.5 Mb!
# 
# n.lf.sf <- list(length = 1000)
# n.sf.lf <- list(length = 1000)
# n.lf.ind <- list(length = 1000)
# n.ind.lf <- list(length = 1000)
# n.sf.ind <- list(length = 1000)
# n.ind.sf <- list(length = 1000)
# 
# for(i in 1:1000) n.sf.ind[[i]] <- maxentRunRawPlot(random.swd.pair('sf', 'ind'))
# for(i in 1:1000) n.ind.sf[[i]] <- maxentRunRawPlot(random.swd.pair('ind', 'sf'))
# for(i in 1:1000) n.lf.ind[[i]] <- maxentRunRawPlot(random.swd.pair('lf', 'ind'))
# for(i in 1:1000) n.ind.lf[[i]] <- maxentRunRawPlot(random.swd.pair('ind', 'lf'))
# for(i in 1:1000) n.lf.sf[[i]] <- maxentRunRawPlot(random.swd.pair('lf', 'sf'))
# for(i in 1:1000) n.sf.lf[[i]] <- maxentRunRawPlot(random.swd.pair('sf', 'lf'))


#-------------------------------------------------------------------------------
# Function to calculate modified Hellinger similarities for a given model run
#-------------------------------------------------------------------------------

I.dist = function(p.x, p.y){
  # Convert the rasters to SpatialGridDataFrame format:
  p.x = as(p.x, 'SpatialGridDataFrame')
  p.y = as(p.y, 'SpatialGridDataFrame')
  # Make a list of the probability surfaces:
  # p.list = list(p.x,p.y)
  # Calculate the modified-Hellinger similarity (Warren 2008, pg. 2870)
  # niche.overlap(p.list)[2,1]
  niche.overlap(list(p.y, p.x))
}

#-------------------------------------------------------------------------------
# Function to run niche equivalency analyses on two flock size classes
#-------------------------------------------------------------------------------
# This function generates a list of two elements, slot1 contains the modified-H
# distance (I) of the actual data and slot2 contains the modified-H of the
# null distribution.

run.nea <- function(sp1, sp2, iterations){ #, null.xy, null.yx){
  I.actual <- I.dist(
    maxentRunRawPlot(dplyr::filter(flockData, sp == sp1) %>%
                       mutate(sp = 1)),
    maxentRunRawPlot(dplyr::filter(flockData, sp == sp2) %>%
                       mutate(sp = 1))
  )[2]
  I.null <- rep(NA, iterations)
  for(i in 1:iterations){
    I.null[i] <- I.dist(
      maxentRunRawPlot(random.swd.pair(sp1, sp2)),
      maxentRunRawPlot(random.swd.pair(sp2, sp1))
    )[2]
  }
  nea.list <- list(I.actual, I.null)
  names(nea.list) <- c('I.actual','I.null')
  return(nea.list)
}

#-------------------------------------------------------------------------------
# Run niche equivalency analyses
#-------------------------------------------------------------------------------

# Large flock vs. small flock:

I.lf.sf <- run.nea('lf','sf',1000)

# Large flock vs. individual sightings (<20 individuals):

I.lf.ind <- run.nea('lf','ind',1000)

# Small flock vs. individual sightings (<20 individuals):

I.sf.ind <- run.nea('sf','ind',1000)

#-------------------------------------------------------------------------------
# Stats for niche equivalency analyses
#-------------------------------------------------------------------------------
# Statistic pairwise comparisons of null and actual modified-Hellinger 
# similarities for large vs. small flocks, large vs. individual sightings and 
# small flocks vs.individual sightings. 

out.stats = function(I.sp1.sp2){
  I.actual = I.sp1.sp2$I.actual
  I.null = I.sp1.sp2$I.null
  I.null.mean = mean(I.null)
  I.null.se = se(I.null)
  I.null.05 = quantile(I.null, probs = 0.05)
  I.null.10 = quantile(I.null, probs = 0.1)
  p = ecdf(I.null)(I.actual)
  l = list(I.actual, I.null,I.null.mean, I.null.se, I.null.05, 
           I.null.10, p)
  names(l) = c('I.actual','I.null', 'I.null.mean', 'I.null.se',
               'I.null.05','I.null.10','p') ; l
}

# LARGE FLOCK VS. SMALL FLOCK:

outStats_I.lf.sf <- out.stats(I.lf.sf)

# LARGE FLOCK VS. INDIVIDUALS:

outStats_I.lf.ind <- out.stats(I.lf.ind)

# SMALL FLOCK VS. INDIVIDUAL:

outStats_I.sf.ind <- out.stats(I.sf.ind)

#-------------------------------------------------------------------------------
# Plot the histograms
#-------------------------------------------------------------------------------
# Histograms compare the density distribution of the modified-Hellinger 
# similarities of the null models against the distance of the actual model for
# large vs. small flocks, large vs. individual sightings and small flocks vs. 
# individual sightings.

# Histogram function:

hist.mhd = function(I.sp1.sp2,main.label,out.name, leg){
  null.dist = out.stats(I.sp1.sp2)[[2]]
  emp.dist = out.stats(I.sp1.sp2)[[1]]
  plot.new()
  # setwd('C:/Users/Brian/Dropbox/rubl_12_15/scratch_out')
  # jpeg(out.name, 1200,1250, res = 300)
  hist(null.dist, breaks = seq(0,1, by = .005), freq = F,
       xlim = c(0.85,1), ylim = c(0,110),
       col = 'gray80', cex.lab = .9,
       main = '', ylab = '',xlab='',
       axes =F)
  axis(1, line = -0.4)
  axis(2, line = 0)
  title(line = 2.5, xlab = 'Modified-Hellinger similarity (I)', ylab = 'Density')
  title(main = main.label, line = .5, cex.main = 1)
  lines(c(emp.dist, emp.dist), c(0,100), lwd = 2, lty = 2)
  if (leg == T)
    legend(.855,90,'Null','gray80',1, bty = 'n',
           x.intersp = .95, cex = .9)
  if (leg == T)
    legend(.85,100,'Actual',lty = 2,lwd = 2, bty = 'n', 
           x.intersp = .4, cex = .9)
  # dev.off()
}

# Make plots:

library(ggplot2)

data.frame(I = I.lf.sf$I.null) %>%
  tbl_df %>%
  ggplot(aes(I)) +
  geom_density(fill = 'gray90') +
  scale_x_continuous(limits = c(.9, 1)) +
  ylab('Density')+
  xlab ('Modified-Hellinger similarity (I)') +
  theme_bw() +
  theme(axis.title.y = element_text(size = rel(1.5), vjust = .9),
        axis.title.x = element_text(size = rel(1.5), vjust = -0.4)) +
  # geom_vline(xintercept = I.lf.sf$I.actual, linewidth = 2.5, linetype = 'longdash')
  geom_segment(data = data.frame(I = I.lf.sf$I.actual),
               aes(x = I, y = 0, xend = I, yend = Inf))#,
# linewidth = 2.5, linetype = 'longdash')

data.frame(I = I.lf.ind$I.null) %>%
  tbl_df %>%
  ggplot(aes(I)) +
  geom_density(fill = 'gray90') +
  scale_x_continuous(limits = c(.9, 1)) +
  ylab('Density')+
  xlab ('Modified-Hellinger similarity (I)') +
  theme_bw() +
  # geom_vline(xintercept = I.lf.sf$I.actual, linewidth = 2.5, linetype = 'longdash')
  geom_segment(data = data.frame(I = I.lf.ind$I.actual),
               aes(x = I, y = 0, xend = I, yend = Inf))#,
#linewidth = 2.5, linetype = 'longdash')


data.frame(I = I.sf.ind$I.null) %>%
  tbl_df %>%
  ggplot(aes(I)) +
  geom_density(fill = 'gray90') +
  scale_x_continuous(limits = c(.9, 1)) +
  ylab('Density')+
  xlab ('Modified-Hellinger similarity (I)') +
  theme_bw() +
  # geom_vline(xintercept = I.lf.sf$I.actual, linewidth = 2.5, linetype = 'longdash')
  geom_segment(data = data.frame(I = I.sf.ind$I.actual),
               aes(x = I, y = 0, xend = I, yend = Inf))#,
#linewidth = 2.5, linetype = 'longdash')

hist.mhd(I.lf.sf, 'Large vs.medium flock sightings', 'mh_dist_lf_sf.jpg',T)

hist.mhd(I.lf.ind, 'Large flock vs. small flock sightings', 'mh_dist_lf_ind.jpg',F)

hist.mhd(I.sf.ind, 'Medium flock vs. small flock sightings', 'mh_dist_sf_ind.jpg',F)

#-------------------------------------------------------------------------------
# Predicted niche occupancy
#===============================================================================

#-------------------------------------------------------------------------------
# PNO functions:
#-------------------------------------------------------------------------------

# Calculate pno:

pno.df <- function(mod.x, mod.y, env.var){
  # Sum the raw probabilities about a given value of an environmental variable:
  pno.df <- data.frame(zonal(mod.x,rStack2009[[env.var]],'sum', digits = 2))
  pno.df[,3] <- zonal(mod.y,rStack2009[[env.var]],'sum', digits = 2)[,2]
  colnames(pno.df) <- c('env','pno.sp1','pno.sp2')
  pno.df
}


# Determine the modified-Hellinger similarity between two pnos:

pno.I <- function(mod.x, mod.y, env.var){
  df <- pno.df(mod.x, mod.y, env.var)
  # Calculate the modified-Hellinger similarity (I):
  niche.overlap(df)[2,1]
}


run.pno <- function(sp1, sp2, env.var, iterations){
  mod.x <- maxentRunRawPlot(dplyr::filter(flockData, sp == sp1) %>%
                              mutate(sp = 1))
  mod.y <- maxentRunRawPlot(dplyr::filter(flockData, sp == sp2) %>%
                              mutate(sp = 1))
  pno.actual <- pno.df(mod.x, mod.y, env.var)
  pno.I.actual <- pno.I(mod.x, mod.y, env.var)
  pno.I.null <- numeric()
  for (i in 1:iterations){
    null.xy <- maxentRunRawPlot(random.swd.pair(sp1, sp2))
    null.yx <- maxentRunRawPlot(random.swd.pair(sp2, sp1))
    pno.I.null[i] = pno.I(null.xy,null.yx, env.var)
  }
  pno.list <- list(pno.I.actual, pno.I.null, pno.actual) # pno.actual, 
  names(pno.list) = c('pno.I.actual','pno.I.null','pno.actual')
  return(pno.list)
}

# Run pno models, returns a list of the
# actual pno-I (one value) and a vector of 100 null pno-I's:

# run.pno = function(mod.x, mod.y, null.xy, null.yx, env.var){
#   pno.actual = pno.df(mod.x, mod.y, env.var)
#   pno.I.actual = pno.I(mod.x, mod.y, env.var)
#   pno.I.null = numeric()
#   for (i in 1:100){
#     pno.I.null[i] = pno.I(null.xy[[i]],null.yx[[i]], env.var)
#   }
#   pno.list = list(pno.I.actual, pno.I.null, pno.actual) # pno.actual, 
#   names(pno.list) = c('pno.I.actual','pno.I.null','pno.actual')
#   return(pno.list)
# }

#-------------------------------------------------------------------------------
# Run PNO
#-------------------------------------------------------------------------------
# THIS TAKES A VERY LONG TIME TO RUN!

pno.lf.sf <- vector('list', length = 15)
for(j in 1:15){
  pno.lf.sf[[j]] <- run.pno('lf', 'sf', j, 100)
}

pno.lf.ind <- vector('list', length = 15)
for(j in 1:15){
  pno.lf.ind[[j]] <- run.pno('lf', 'ind', j, 100)
}

pno.sf.ind <- vector('list', length = 15)
for(j in 1:15){
  pno.sf.ind[[j]] <- run.pno('sf', 'ind', j, 100)
}


# 
# pno.lf.sf = list()
# for(i in 1:15){
#   pno.lf.sf[[i]] = run.pno(lf, sf, n.lf.sf, n.sf.lf,i)
# }
# 
# pno.lf.ind = list()
# for(i in 1:15){
#   pno.lf.ind[[i]] = run.pno(lf, ind, n.lf.ind, n.ind.lf,i)
# }
# 
# pno.sf.ind = list()
# for(i in 1:15){
#   pno.sf.ind[[i]] = run.pno(sf, ind, n.sf.ind, n.ind.sf,i)
# }

names(pno.lf.sf) = names(rStack2009)
names(pno.lf.ind) = names(rStack2009)
names(pno.sf.ind) = names(rStack2009)


saveRDS(pno.lf.sf, 'pno.lf.sf.rds')
saveRDS(pno.lf.ind, 'pno.lf.ind.rds')
saveRDS(pno.sf.ind, 'pno.sf.ind.rds')

#-------------------------------------------------------------------------------
# Stat output of pno
#-------------------------------------------------------------------------------

pno.stats = function(pno.sp1.sp2.env){
  I.actual = pno.sp1.sp2.env$pno.I.actual
  I.null = pno.sp1.sp2.env$pno.I.null
  I.null.mean = mean(I.null)
  I.null.se = se(I.null)
  I.null.05 = quantile(I.null, probs = 0.05)
  I.null.10 = quantile(I.null, probs = 0.1)
  p = ecdf(I.null)(I.actual)
  l = list(I.actual, I.null,I.null.mean, I.null.se, I.null.05, 
           I.null.10, p)
  names(l) = c('I.actual','I.null', 'I.null.mean', 'I.null.se',
               'I.null.05','I.null.10','p') ; l
}

out.pno.stats = function(pno.sp1.sp2){
  pno.out = list()
  for (i in 1:15){
    pno.out[[i]] = pno.stats(pno.sp1.sp2[[i]])
  }
  names(pno.out) = names(rStack2009)
  pno.out
}

out.pno.lf.sf = out.pno.stats(pno.lf.sf)
out.pno.lf.ind = out.pno.stats(pno.lf.ind)
out.pno.sf.ind = out.pno.stats(pno.sf.ind)

#-------------------------------------------------------------------------------
# Plotting PNO
#-------------------------------------------------------------------------------

pno.lf.ind[['dev_hi']][[1]][,3]

names(rStack2009)[env]

plot.pno = function(env) {
  df = pno.lf.sf[[env]][[1]]
  ind = pno.lf.ind[[env]][[1]][,3]
  env = df[,1]
  lf = df[,2]
  sf = df[,3]
  plot(lf~env, type = 'l', xlim = c(0,1), lwd =2, bty ='l',
       main = names(rStack2009[[env]]), ylab = 'PNO')
  lines(env, sf, lty = 2, lwd = 2)
  lines(env, ind, lty = 3, lwd = 2)
}

plot.pno(3)


plot.pno = function(env) {
  df = pno.lf.sf[[env]][[3]]
  ind = pno.lf.ind[[env]][[3]][,3]
  lc = df[,1]
  lf = df[,2]
  sf = df[,3]
  plot(lf~lc, type = 'l', xlim = c(0,1), lwd =2, bty ='l',
       main = names(rStack2009)[env], ylab = 'PNO')
  lines(lc, sf, lty = 2, lwd = 2)
  lines(lc, ind, lty = 3, lwd = 2)
}

plot.pno(3)

#########################################################################################################
#########################################################################################################
#########################################################################################################
# Compare distributions
library(ggplot2)

ptsENV <- swdBG %>%
  dplyr::mutate(sp = 'bg', observationType = 'ebird') %>%
  dplyr::select(sp, observationType, lon:woodland) %>%
  dplyr::bind_rows(swdRUBL$eb %>%
                     dplyr::mutate(observationType = 'ebird') %>%
                     dplyr::select(sp, observationType, lon:woodland)) %>%
  dplyr::bind_rows(swdRUBL$bz %>%
                     dplyr::mutate(observationType = 'blitz') %>%
                     dplyr::select(sp, observationType, lon:woodland)) %>%
  dplyr::mutate(sp = factor(sp))

ptsENV$sp <- factor(sp, levels =  'bg','ind','sf', 'lf')


modeFun <- function(variable, rounding){
  d <- density(round_any(variable, rounding))
  dFrame <- data.frame(d$x, d$y)
  names(dFrame) <- c('x','y')
  dFrame %>% dplyr::filter(y == max(y)) %>%
    dplyr::select(x)
}

summarizeByFlockSize <- function(variable, rounding){
  outTable <- ptsENV %>%
    dplyr::select(sp, observationType, 
                  variable = matches(variable)) %>%
    dplyr::group_by(sp) %>%
    dplyr::summarize(min = min(variable),
                     max = max(variable),
                     mean = mean(variable),
                     variance = var(variable),
                     IQR = IQR(variable))
  if(variable == 'ppt' | variable == 'tmin'){
    outTable <- ptsENV %>%
      dplyr::select(sp, observationType, 
                    variable = matches(variable)) %>%
      dplyr::group_by(sp) %>%
      dplyr::summarize(min = min(variable),
                       max = max(variable),
                       mean = mean(variable),
                       variance = var(variable),
                       IQR = IQR(variable),
                       mode = as.numeric(modeFun(variable, 0.5)))
  }
  return(outTable)
}


ksTest <- function(variable, sp1, sp2){
  v1 <- ptsENV %>% 
    dplyr::filter(sp == sp1) %>%
    collect %>%
    .[[variable]]
  v2 <- ptsENV %>% 
    dplyr::filter(sp == sp2) %>%
    collect %>%
    .[[variable]]
  ks.test(v1, v2)
}

summarizeByFlockSize('tmin', 0.5)


ksTest('tmin','sf','lf')

ksTest('dev_hi','sf','ind')


ksTest('dev_li','sf','lf')

ksTest('flood','ind','bg')

ksTest('flood','lf','ind')

ksTest('flood','lf','sf')



ggplot(ptsENV, 
       aes(x = flood, fill = sp)) +
  geom_histogram(aes(y=0.5*..density..),
                 position = 'identity', binwidth=0.1) +
  facet_wrap(~sp,nrow=3) +
  theme(aspect.ratio = 1) +
  theme_bw()

ggplot(ptsENV, 
       aes(factor(sp,levels = c('bg','ind','sf','lf')),
           dev_hi, fill = sp)) +
  geom_violin(adjust = 0.5) +
  xlab('Point data class') +
  ylab('Environmental variable') +
  coord_flip() +
  theme_bw()


