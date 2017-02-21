# Maxent, running winter data
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

smartInstallAndLoad(c('dplyr', 'tidyr','stringi','stringr', 'sp',
                      'lubridate','raster', 'dismo', 'ggplot2'))

# Override some libraries for tidyverse functions:

filter <- dplyr::filter
select <- dplyr::select

# IMPORTANT! Prior to running the below, run the file "makeSWD.R"
# Also, in order to run this script, you need to place the maxent jar file in the folder: win-library/[version]/dismo/java

# Set outPlots directory (example, do not run):

outPlotsDir <- 'C:/Users/Brian/Desktop/gits/blitzAnalysis/outPlots/'

# Set path to climate data:

pathToClimateData <- '/Users/bsevans/Dropbox/rustyBlackbirdData/climateRasters/'

#-------------------------------------------------------------------------------*
# ---- Basic functions ---
#-------------------------------------------------------------------------------*

se <- function(x) {sd(x)/sqrt(length(x))}

#-------------------------------------------------------------------------------*
# Get swd (across observation types):
#-------------------------------------------------------------------------------*

swdSmall <- prepSWD(1,19, 2009:2011, protocolChoice = 'all')
swdMedium <- prepSWD(20,99, 2009:2011, protocolChoice = 'all')
swdLarge <- prepSWD(100,Inf, 2009:2011, protocolChoice = 'all')

#===============================================================================*
# ---- MODEL RUNNING AND CALIBRATION ----
#===============================================================================*

# Basic function to run a model:

maxentRun <- function(swd, betaMultiplier,
                      kFold = 'noCrossValidate', excludeVariables = NULL){
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
  # Set model arguments:
  modArguments <- c('nothreshold', 'nohinge', 'noproduct','noquadratic',
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

#-------------------------------------------------------------------------------*
# ---- Calculate AIC ----
#-------------------------------------------------------------------------------*

# Function to calculate AICc, given a specific SWD and beta multiplier:

calcAIC <- function(swd, betaMultiplier) {
  # Extract a  model:
  model <- maxentRunReduced(swd, betaMultiplier)
  # Extract lambdas file and convert to a data frame
  lambdas <- model@lambdas
  lambda.df <- data.frame(do.call('rbind', strsplit(lambdas,',',fixed=TRUE)))
  # Determing the number of parameters that had a lambda of zero 
  # Note: Using n.lz = Number of Lambdas that equals Zero
  n.lz <- length(lambda.df[as.numeric(as.character(lambda.df[,2])) == 0,1])
  # Calculate the number of model parameters
  kn <- length(lambdas)- n.lz - 4
  # Predict suitability:
  presencePredictions <- predict(model,
                                 swd %>% filter(pa == 1),
                                 args=c('outputformat=raw'))
  absencePredictions <- predict(model, 
                                swd %>% filter(pa == 0),
                                args=c('outputformat=raw'))
  probsum <- sum(presencePredictions, absencePredictions)
  # How many points?
  n <- length(presencePredictions)
  # Get log likelihood:
  loglik <- sum(log(presencePredictions / probsum))
  # Calculate the AICc
  AICc <- -2*loglik + 2*kn + ((2*kn*(kn+1))/(n-kn-1))
  # Output
  aicFrame <- data.frame(nParm = kn, beta = betaMultiplier,AICc = AICc)
  return(aicFrame)
}

#-------------------------------------------------------------------------------*
# ---- Assess aic values across beta values and output results table -----
#-------------------------------------------------------------------------------*

betaFinder <- function(swd, betaValues){
  out <- data.frame(matrix(nrow = length(betaValues),ncol = 3))
  for (i in 1:length(betaValues)){
    out[i,] <- calcAIC(swd, betaValues[i])
  }   
  colnames(out) <- c('nparm', 'beta','AICc')
  out
}

# Running betaFinder to determine the best model:

betaValues <- seq(0, 10, .1)

betaSmall <- betaFinder(swdSmall, betaValues)

betaMedium <- betaFinder(swdMedium, betaValues)
  
betaLarge <- betaFinder(swdLarge, betaValues)

# OUTPUT (because that takes a while to run):
# Small = 2.7, medium = 1.4, large = 1.0

#===============================================================================*
# ---- MODEL EVALUATION ----
#===============================================================================*

#-------------------------------------------------------------------------------*
# ---- Get AUC values associated with test points ----
#-------------------------------------------------------------------------------*

# Get best models for each flock size:

bestModelSmall <- maxentRunReduced(swdSmall, betaMultiplier = 2.7)
bestModelMedium  <- maxentRunReduced(swdMedium, betaMultiplier = 1.4)
bestModelLarge <- maxentRunReduced(swdLarge, betaMultiplier = 1)

# Function to run maxent model and return AUC for a given cross-validation fold:

runMaxentAUC <- function(swd, bestModel, betaMultiplier, kFold){
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
  modArguments <- c('nothreshold', 'nohinge', 'noproduct','noquadratic',
                    str_c('betamultiplier=', betaMultiplier),
                    'addallsamplestobackground',
                    'writebackgroundpredictions',
                    'noautofeature','nooutputgrids',
                    'maximumiterations=10000', 'verbose')
  # Run maxent model with training and background data:
  maxentModel <- maxent(swdTrain[,-1], swdTrain[,1], args = modArguments)
  # Predict model values at test points:
  predictionPresence <- dismo::predict(maxentModel, swdTest)
  predictionAbsence <- dismo::predict(maxentModel,swdAbsence %>%
                                 select(-c(pa, k)))
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

# Get AUC summary stats across observation types:

summarySmall <- getAUC(swdSmall, bestModelSmall, 2.7) %>%
  mutate(flockSize = 'Small')

summaryMedium <- getAUC(swdMedium, bestModelMedium, 1.4) %>%
  mutate(flockSize = 'Medium')

summaryLarge  <- getAUC(swdLarge, bestModelLarge, 1) %>%
  mutate(flockSize = 'Large')

summaryAUC <- bind_rows(summarySmall, summaryMedium, summaryLarge)

# Plot AUC (not by observation type):

ggplot(summaryAUC, aes(x = flockSize, y = meanAUC)) +
  geom_errorbar(aes(ymin=meanAUC - 1.96*seAUC,
                    ymax=meanAUC + 1.96*seAUC), width = 0) +#, position=pd)
  geom_point() +
  theme_bw() +
  ylab('AUC') +
  xlab('Flock size class') +
  geom_line(size = 1) + coord_flip()

# Calculate by auc summary stats by observation type:

summarySmallEb <- getAUC(
  prepSWD(1, 19, 2009:2011, 'eb'),
  bestModelSmall,
  2.7
) %>%
  mutate(flockSize = 'Small',
         encounterType = 'eBird')

summarySmallBlitz <- getAUC(
  prepSWD(1, 19, 2009:2011, 'blitz'),
  bestModelSmall,
  2.7
) %>%
  mutate(flockSize = 'Small',
         encounterType = 'Blitz') 

summaryMediumEb <- getAUC(
  prepSWD(1, 19, 2009:2011, 'eb'),
  bestModelMedium,
  1.4
) %>%
  mutate(flockSize = 'Medium',
         encounterType = 'eBird')

summaryMediumBlitz <- getAUC(
  prepSWD(1, 19, 2009:2011, 'blitz'),
  bestModelMedium,
  1.4
) %>%
  mutate(flockSize = 'Medium',
         encounterType = 'Blitz') 

summaryLargeEb <- getAUC(
  prepSWD(1, 19, 2009:2011, 'eb'),
  bestModelLarge,
  1
) %>%
  mutate(flockSize = 'Large',
         encounterType = 'eBird')

summaryLargeBlitz <- getAUC(
  prepSWD(1, 19, 2009:2011, 'blitz'),
  bestModelLarge,
  1
) %>%
  mutate(flockSize = 'Large',
         encounterType = 'Blitz') 

# Bind summaries into a single data frame:

summaryAUCbyObservationType <- summaryAUC %>%
  mutate(encounterType = 'All') %>%
  bind_rows(summarySmallEb, summarySmallBlitz,
            summaryMediumEb, summaryMediumBlitz,
            summaryLargeEb, summaryLargeBlitz)

# Example plot of AUC by observation:

dodge <- position_dodge(.4)
ann_text <- data.frame(mpg = 15,wt = 5,lab = "Text",
                       cyl = factor(8,levels = c("4","6","8")))

aucPlot <- ggplot(
  summaryAUCbyObservationType %>%
    mutate(
      encounterType = str_replace_all(encounterType,'All', 'Combined'),
      encounterType = str_replace_all(encounterType, 'Blitz',
                                      'Winter Blitz')),
  aes(x = encounterType, y = meanAUC)) +
  geom_errorbar(aes(ymin=meanAUC - 1.96*seAUC,
                    ymax=meanAUC + 1.96*seAUC), width = 0,
                size = .75) +
  geom_point(size = 3) +
  theme_bw() +
  ylab('AUC') +
  xlab('Observation class') +
  geom_line(size = 1) +
  theme(axis.title = element_text(margin=margin(0,10,0,0))) +
  theme(legend.key = element_blank()) + 
  facet_grid(flockSize ~ .)

#-------------------------------------------------------------------------------*
# ---- Get ROC and plot ----
#-------------------------------------------------------------------------------*

# Function to get a maxent evaluation frame:

maxentEvaluation <- function(swd, bestModel, betaMultiplier, kFold){
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
    filter(pa == 0) %>% 
    data.frame
  swdTrain <- swdReduced %>%
    filter(k != kFold & pa == 1) %>%
    bind_rows(swdAbsence) %>%
    select(-k) %>% 
    data.frame
  swdTest <- swdReduced %>%
    filter(k == kFold & pa == 1) %>%
    # bind_rows(swdAbsence) %>%
    select(-k) %>% 
    data.frame
  # Set model arguments
  modArguments <- c('nothreshold', 'nohinge', 'noproduct','noquadratic',
                    str_c('betamultiplier=', betaMultiplier),
                    'addallsamplestobackground',
                    'writebackgroundpredictions',
                    'noautofeature','nooutputgrids',
                    'maximumiterations=10000', 'verbose')
  # Run maxent model with training and background data:
  maxentModel <- maxent(swdTrain[,-1], swdTrain[,1], args = modArguments)
  # Predict model values at test points:
  predictionPresence <- dismo::predict(maxentModel, swdTest)
  predictionAbsence <- dismo::predict(maxentModel,swdAbsence %>%
                                        select(-c(pa, k)))
  # Evaluate model:
  evaluationObject <- dismo::evaluate(p = predictionPresence,
                                      a = predictionAbsence)
  # Return evaluation object:
  return(evaluationObject)
}

# Function to make an ROC frame (true and false positive rates):

getRocFrame <- function(eObject, k){
  logisticThresh <- seq(0, 1, .01) #%>% sort(decreasing = TRUE)
  rocFrame <- data.frame(logisticThresh)
  for(i in 1:length(logisticThresh)){
    fpr <- sum(eObject@absence >= logisticThresh[i])/length(eObject@absence)
    tpr <- sum(eObject@presence >= logisticThresh[i])/length(eObject@presence)
    rocFrame$fpr[i] <- fpr
    rocFrame$tpr[i] <- tpr
  }
  rocFrame$k <- k
  return(rocFrame)
}

# Get and plot roc for large flocks as an example:

rocList <- vector('list', length = 5)

for(i in 1:5){
  evaluationObject <- maxentEvaluation(swdLarge,bestModelLarge,1, i)
  rocList[[i]] <- getRocFrame(evaluationObject, i)
}

rocFrame <- bind_rows(rocList)

rocPlot <- ggplot(rocFrame %>%
         select(fpr, tpr, k) %>%
         mutate(fpr = round(fpr, 4),
                tpr = round(tpr, 4)) %>%
         distinct, 
       aes(y = tpr, x = fpr,group = k)) +
  geom_line(size = .6) +
  geom_line(data = data.frame(x = seq(0, 1, .5), k = 1), aes(x = x, y = x),
            linetype = 4, size = .6) +
  theme_bw() +
  ylab('Sensitivity') +
  xlab('1 - Specificity (Fractional predicted area)') +
  theme(axis.title = element_text(margin=margin(0,10,0,0)))

# Plot and save output:

png(filename = paste0(outPlotsDir, 'rocAUC.png'),
    width = 7.5, height = 4, units = 'in', res = 300)
grid.arrange(
  arrangeGrob(
    rocPlot +
      geom_text(
        data = data.frame(
          fpr = 0.05,
          tpr = .97,
          k = 1),
        label = "A)",
        size = 6), 
    aucPlot +
      geom_text(
        data = data.frame(
          flockSize = 'Large',
          encounterType = 3.1,
          meanAUC = .85),
        label = "B)",
        size = 6) +
      theme(axis.text.y = element_text(size = 6),
            axis.text.x = element_text(size = 8)),
    ncol=2, widths=c(1.6,.9))
)
dev.off()

#===============================================================================*
# ---- Area, prevalence, threshold ----
#===============================================================================*

# Function to run model for area and prevalence, given a threshold:

runMaxentAreaPrevalence <- function(swd, bestModel, betaMultiplier){
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
  swdTrain <- swdReduced %>%
    select(-k) %>%
    data.frame
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

# Function to get a fractional predicted area frame for a given model:

getAreaPrevalenceFrame <- function(model, i){
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
      iteration = i,
      stat = ifelse(str_detect(stat, 'Prevalence'), 'prevalence', stat),
      stat = ifelse(str_detect(stat, 'logistic'), 'threshold', stat),
      stat = ifelse(str_detect(stat, 'area'), 'area', stat)
      ) %>%
    spread(stat, value)
}

# Function to generate a random sample of flock size 1 by flock size 2 Note:
# Flock size 1 is the flocks size of interest, comparison to flock size 2

randomSWDpair <- function(swd1, swd2){
  nSamples <- nrow(filter(swd1, pa == 1))
  # Get background points:
  bgPoints <- filter(swd1, pa == 0)
  # Get random sample of presence points :
  bind_rows(
    filter(swd1, pa == 1),
    filter(swd2, pa == 1)
    ) %>%
    sample_n(nSamples,
      replace = FALSE) %>%
     bind_rows(bgPoints)
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
  saveLoc <- 'C:/Users/Brian/Dropbox/randomRUBLmodels/'
  saveRDS(modLS, paste0(saveLoc, 'modLS', i, '.rds'))
  saveRDS(modLM, paste0(saveLoc, 'modLM', i, '.rds'))
  saveRDS(modMS, paste0(saveLoc, 'modMS', i, '.rds'))
}

# Provide an example RDS location (example, do not run):

RDSlocation <- 'C:/Users/Brian/Dropbox/randomRUBLmodels'

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
      mutate(permuted = modCombos[i])
  }
  outList[[i]] <- bind_rows(subOutList)
}

permutationFrame <- bind_rows(outList)

# Example plot of predicted area by flock size (large):

predictedAreaLFplot <- ggplot(permutationFrame %>%
         mutate(
           permuted = str_replace(permuted, 'modLM', 'Large-Medium'),
           permuted = str_replace(permuted, 'modLS', 'Large-Small'),
           permuted = str_replace(permuted, 'modMS', 'Medium-Small'),
           Permutation = permuted
           ) %>%
         filter(permuted != 'Medium-Small'),
       aes(x=area)) + 
  geom_density(aes(group=Permutation, fill=Permutation), alpha=0.7) +
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

#===============================================================================*
# ---- Evaluation frames ----
#===============================================================================*

# Next: k-fold for each model for getting error in predictions

# Function to run maxent for a k-fold partition:

runMaxentK <- function(swd, bestModel, betaMultiplier, kFold){
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
    select(-k)
  swdTest <- swdReduced %>%
    filter(k == kFold & pa == 1) %>%
    bind_rows(swdAbsence) %>%
    select(-k)
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

# Run for small medium and large flocks:

outListSmall <- vector('list', length = 5)
outListMedium <- vector('list', length = 5)
outListLarge <- vector('list', length = 5)

for(i in 1:5){
  outListSmall[[i]] <- runMaxentK(swdSmall, bestModelSmall, 2.7, i)
  outListMedium[[i]] <- runMaxentK(swdMedium, bestModelMedium, 1.4, i)
  outListLarge[[i]] <- runMaxentK(swdLarge, bestModelLarge, 1, i)
}

# Function to get summary stats for the above:

getSummaryStatsK <- function(outList, flockSize){
  results <- outList@results %>% 
    data.frame
  results$stat <- row.names(results)
  row.names(results) <- NULL
  names(results) <- c('value', 'stat')
  resultsFrame <- results %>%
    filter(str_detect(stat,'Prevalence')|
           str_detect(stat, 'Maximum.training.sensitivity.plus.')) %>%
    filter(!str_detect(stat, '.cumulative'),
           !str_detect(stat, 'omission'))
  resultsFrame$stat <- c('prevalence', 'logThresh', 'area')
  resultsFrame$flockSize <- flockSize
  return(resultsFrame)
}

# Run to get summary stats for the above:

summaryListSmall <- vector('list', length = 5)
summaryListMedium <- vector('list', length = 5)
summaryListLarge <- vector('list', length = 5)

for(i in 1:5){
  summaryListSmall[[i]] <- getSummaryStatsK(outListSmall[[i]], 'small')
  summaryListMedium[[i]] <- getSummaryStatsK(outListMedium[[i]],  'medium')
  summaryListLarge[[i]] <- getSummaryStatsK(outListLarge[[i]], 'large')
}

# Get flock summaries:

summaryByFlock <- bind_rows(summaryListSmall, summaryListMedium, summaryListLarge) %>%
  group_by(flockSize, stat) %>%
  summarize(meanValue = mean(value),
            seValue = se(value),
            minCi = meanValue - seValue * 1.96,
            maxCi = meanValue + seValue * 1.96)

# Function to get an evaluation object:

modelEvaluate <- function(minFlockSize, maxFlockSize, years, protocolChoice = 'all'){
  presenceLC <- prepSWD(minFlockSize,maxFlockSize, years, protocolChoice = 'all') %>%
    filter(pa == 1) %>%
    select(dev_hi:tmin) %>%
    as.matrix
  absenceLC <- prepSWD(minFlockSize,maxFlockSize, years, protocolChoice = 'all') %>%
    filter(pa == 0) %>%
    select(dev_hi:tmin) %>%
    as.matrix
  evaluate(presenceLC, absenceLC, test14[[2]])
}

# Function to run model and get evaluation object:

maxentEvaluate <- function(observationClass, flockSizeClass, kFold, beta.multiplier){
  # Run model
  maxentModelOut <- maxentRun(observationClass, flockSizeClass, kFold, beta.multiplier)
  maxentModel <- maxentModelOut$maxentModel
  # Create test point files for evaluation
  swdTest <-   swd[[observationClass]][[flockSizeClass]] %>%
    filter(k == kFold)
  testPresence <- swdTest %>%
    filter(sp == 1) %>%
    select(-c(sp,lon, lat,k))
  testAbsence <- swdTest %>%
    filter(sp == 0) %>%
    select(-c(sp,lon, lat,k))
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

#===============================================================================*
# ---- Variable contribution ----
#===============================================================================*

# Functions to extract variable contribution for a given model:

readLambdaFileIteration <- function(model, i){
  lambdaData <- model[[1]][[i]][[1]]@lambdas[1:15]
  tf <- tempfile()
  writeLines(lambdaData, tf)
  read.csv(tf, fill = TRUE, header = F) %>%
    select(variable = V1, lambda = V2)
}

# Function to make lambda contribution frame:

makeLambdaContributionFrame <- function(model){
  dfContributionMat <- matrix(nrow = 15, ncol = 5)
  lambdaMat <- matrix(nrow = 15, ncol = 5)
  for(i in 1:ncol(lambdaMat)){
    lambdaMat[,i] <- readLambdaFileIteration(model, i)[,2]
    dfContributionMat[,i] <- model[[1]][[i]][[1]]@results[7:21]
  }
  # Lambda summary data:
  lambdaSummary <- data.frame(
    variable = as.character(readLambdaFile(model,1)$variable), 
    meanLambda = apply(lambdaMat, 1, mean),
    seLambda = apply(lambdaMat, 1, se))
  # Summary contribution frame
  variableContributionSummary <- data.frame(
    variable = as.character(readLambdaFile(model,1)$variable), 
    meanContribution = apply(dfContributionMat, 1, mean),
    seContribution = apply(dfContributionMat, 1, se))
  # Output:
  outFrame <- plyr::join(variableContributionSummary, lambdaSummary) %>%
    arrange(desc(meanContribution))
  return(outFrame)
}

lambdaContributionFrame_allIndOut <- makeLambdaContributionFrame(allIndOut)

lambdaContributionFrame_allSfOut <- makeLambdaContributionFrame(allSfOut)

lambdaContributionFrame_allLfOut <- makeLambdaContributionFrame(allLfOut)

#===============================================================================*
# ---- Niche equivalency analysis ----
#===============================================================================*

flockData <- swdRUBL[[1]] 

# Generate a random sample of species 1 by species 2
# Note: Species 1 is the species of interest, comparison to species 2

random.swd.pair <- function(sp1, sp2){
  # Remove unevaluated flock
  flockData <- flockData %>%
    filter(sp == sp1|sp == sp2)
  nSp1 <- nrow(filter(flockData, sp == sp1))
  # Determine the sample size as the proportion of species 1:
  prob.s1 = nSp1/nrow(flockData)
  # Generate random value of 1 or 0 with the probability of obtaining
  # a 1 related to the # of s1 observations:
  flockData %>%
    mutate(s.rand = rbinom(nrow(.),1,prob.s1)) %>%
    filter(s.rand == 1) %>%
    mutate(sp = 1) %>%
    select(-s.rand)
}

# Function to run maxent for a given dataset:

maxentRunRawPlot <- function(inFlockData, beta.multiplier = 0){
  # Model input:
  swdFlock <- inFlockData
  max.in <- rbind(swdFlock, swdBG) %>%
    select(-c(lat, lon))
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

lf <- maxentRunRawPlot(filter(flockData, sp == 'lf') %>%
                         mutate(sp = 1))
sf <- maxentRunRawPlot(filter(flockData, sp == 'sf') %>%
                         mutate(sp = 1))
ind <- maxentRunRawPlot(filter(flockData, sp == 'ind') %>%
                          mutate(sp = 1))

# Function to calculate modified Hellinger similarities for a given model run:

I.dist <- function(p.x, p.y){
  # Convert the rasters to SpatialGridDataFrame format:
  p.x <- as(p.x, 'SpatialGridDataFrame')
  p.y <- as(p.y, 'SpatialGridDataFrame')
  # Calculate the modified-Hellinger similarity (Warren 2008, pg. 2870)
  niche.overlap(list(p.y, p.x))
}

# Function to run niche equivalency analyses on two flock size classes. 
# Note: This function generates a list of two elements, slot1 contains the 
# modified-H  distance (I) of the actual data and slot2 contains the modified-H
# of the null distribution.

run.nea <- function(sp1, sp2, iterations){
  I.actual <- I.dist(
    maxentRunRawPlot(filter(flockData, sp == sp1) %>%
                       mutate(sp = 1)),
    maxentRunRawPlot(filter(flockData, sp == sp2) %>%
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

# Run niche equivalency analyses:

# Large flock vs. small flock:

I.lf.sf <- run.nea('lf','sf',1000)

# Large flock vs. individual sightings (<20 individuals):

I.lf.ind <- run.nea('lf','ind',1000)

# Small flock vs. individual sightings (<20 individuals):

I.sf.ind <- run.nea('sf','ind',1000)

# Stats for niche equivalency analyses.
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

#-------------------------------------------------------------------------------*
# ---- Plot the niche equivalency histograms ----
#-------------------------------------------------------------------------------*
# Histograms compare the density distribution of the modified-Hellinger 
# similarities of the null models against the distance of the actual model for
# large vs. small flocks, large vs. individual sightings and small flocks vs. 
# individual sightings.

# Plot Niche equivalency analyses:

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
  geom_segment(data = data.frame(I = I.lf.sf$I.actual),
               aes(x = I, y = 0, xend = I, yend = Inf))

data.frame(I = I.lf.ind$I.null) %>%
  tbl_df %>%
  ggplot(aes(I)) +
  geom_density(fill = 'gray90') +
  scale_x_continuous(limits = c(.9, 1)) +
  ylab('Density')+
  xlab ('Modified-Hellinger similarity (I)') +
  theme_bw() +
  geom_segment(data = data.frame(I = I.lf.ind$I.actual),
               aes(x = I, y = 0, xend = I, yend = Inf))

data.frame(I = I.sf.ind$I.null) %>%
  tbl_df %>%
  ggplot(aes(I)) +
  geom_density(fill = 'gray90') +
  scale_x_continuous(limits = c(.9, 1)) +
  ylab('Density')+
  xlab ('Modified-Hellinger similarity (I)') +
  theme_bw() +
  geom_segment(data = data.frame(I = I.sf.ind$I.actual),
               aes(x = I, y = 0, xend = I, yend = Inf))

#===============================================================================*
# ---- Make model predictions (for raster maps) ----
#===============================================================================*

# Function to get environmental rasters for a given year:

getRstack <- function(year){
  # Add tmin and ppt for a given year:
  rStack[['tmin']] <- raster(paste0(pathToClimateData, 'tmin',year))
  rStack[['tmin2']] <- rStack[['tmin']]^2
  rStack[['ppt']] <-  raster(paste0(pathToClimateData, 'ppt',year))
  return(rStack)
}

# Function to get prediction of a given best model and year:

getLogisticPrediction <- function(bestModel, rasterStack){
  predict(bestModel, rasterStack,
          args='outputformat=logistic', 
          progress='text')
}

# Get logistic predictions (2009 example):

rStack2009 <- getRstack(2009)

small2009 <- getLogisticPrediction(bestModelSmall, rStack2009)
medium2009 <- getLogisticPrediction(bestModelMedium, rStack2009)
large2009 <- getLogisticPrediction(bestModelLarge, rStack2009)

# To plot these, go to script plotSuitabilityMaps.R

#### END ####