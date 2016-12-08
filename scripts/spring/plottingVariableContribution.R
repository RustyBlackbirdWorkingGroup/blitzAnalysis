#-----------------------------------------------------------------------------------*
# ---- SET-UP ----
#===================================================================================*

library(dplyr) ; library(tidyr) ; library(ggplot2) 
library(grid) ; library(gridExtra) ; library(stringr)

# STACKED BAR:

# Theme:

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

#-----------------------------------------------------------------------------------*
# ---- DATA ----
#===================================================================================*
# Data are derived in script: maxentScript.R
# My working directory: C:/Users/Brian/Desktop/gits/RUBL

flockNames <- c('Small', 'Medium', 'Large')

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

dataStackedBar <- bind_rows(variableContributionByPeriod) %>%
  tbl_df %>%
  mutate(variable = ifelse(variable == 'tmin2', 'tmin', variable)) %>%
  # Make the tmin contribution the sum of tmin and tmin2:
  group_by(flockSize, samplingPeriod, variable) %>%
  summarize(contribution = sum(contribution)) %>%
  ungroup %>%
  #---------------------------------------------------------------------------------*
  # Set some variables to "other":
  mutate(variable = as.character(variable),
         variable = ifelse(
    !(variable %in% c('tmin', 'flood', 'rowcrop', 'ppt', 
                      'flood', 'pasture', 'dev_hi','weth', 'wetw')),
    'other', variable)
  ) %>%
  #---------------------------------------------------------------------------------*
  # Sum the other percentages for each flock size class:
  group_by(samplingPeriod, flockSize, variable) %>%
  summarize(varContribution = sum(contribution)) %>%
  ungroup %>%
  #---------------------------------------------------------------------------------*
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

#-----------------------------------------------------------------------------------*
# ---- PLOTTING ----
#===================================================================================*

# names(colorBlindPalette) <- unique(dataStackedBar$variable)

plotWithTemp <- ggplot(dataStackedBar %>%
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
        # legend.margin = unit(0, 'line'),
        # legend.key.size = unit(2, 'lines'),
        legend.position = 'top') +
  facet_wrap(~samplingPeriod, nrow = 1)

colorBlindPalette2 <- c("#E69F00",  "#009E73", "#F0E442", 
                         "#0072B2", "#D55E00",'#9C661F','#68228B',"gray60")

plotNoTemp <- ggplot(dataStackedBar %>% 
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

plots <- list(plotWithTemp, plotNoTemp)

g <- ggplotGrob(plots[[1]] + theme(legend.position="top"))$grobs

legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]


# Plot output:

png(filename = "C:/Users/Brian/Desktop/gits/blitzAnalysis/outPlots/VariableContributionSpring.png", 
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


