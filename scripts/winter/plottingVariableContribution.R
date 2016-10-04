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
                       "#0072B2", "#D55E00","#56B4E9", "gray60")

#-----------------------------------------------------------------------------------*
# ---- DATA ----
#===================================================================================*
# Data are derived in script: maxentScript.R
# My working directory: C:/Users/Brian/Desktop/gits/RUBL

dataStackedBar <- bind_rows(
  # Bind the three contribution frames for each flock size class:
  getVariableContribution(bestModelSmall) %>%
    mutate(flockSize = 'Small'),
  getVariableContribution(bestModelMedium) %>%
    mutate(flockSize = 'Medium'),
  getVariableContribution(bestModelLarge) %>%
    mutate(flockSize = 'Large')
  ) %>%
  mutate(variable = ifelse(variable == 'tmin2', 'tmin', variable)) %>%
  # Make the tmin contribution the sum of tmin and tmin2:
  group_by(flockSize, variable) %>%
  summarize(contribution = sum(contribution)) %>%
  ungroup %>%
  #---------------------------------------------------------------------------------*
  # Set some variables to "other":
  mutate(variable = as.character(variable),
         variable = ifelse(
    !(variable %in% c('tmin', 'flood', 'rowcrop', 'ppt', 
                      'flood', 'pasture', 'dev_hi')),
    'other', variable)
  ) %>%
  #---------------------------------------------------------------------------------*
  # Sum the other percentages for each flock size class:
  group_by(flockSize, variable) %>%
  summarize(varContribution = sum(contribution)) %>%
  ungroup %>%
  #---------------------------------------------------------------------------------*
  # Reset the factor levels and labels for plotting:
  mutate(
    variable = variable %>%
      factor(
        levels = c('flood','pasture','rowcrop','ppt','dev_hi','tmin','other'),
        labels = c('Floodplain    ', 'Pasture    ',
                   'Row crop    ','Precipitation    ', 'Highly developed    ',
                   'Minimum temperature    ', 'Other    '),
             ordered = TRUE)
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
        legend.position = 'top')

colorBlindPalette2 <- c("#E69F00",  "#009E73", "#F0E442", 
                         "#0072B2", "#D55E00","gray60")

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
  theme(legend.key.height = unit(2,"line"))

plots <- list(plotWithTemp, plotNoTemp)

g <- ggplotGrob(plots[[1]] + theme(legend.position="top"))$grobs

legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]


# Plot output:

png(filename = "C:/Users/Brian/Desktop/gits/blitzAnalysis/outPlots/VariableContributionMultiPlot.png", 
    width = 7.5, height = 4, units = 'in', res = 300)
gridExtra::grid.arrange(legend, gridExtra::arrangeGrob(plotWithTemp + theme(legend.position = 'none'),
                                 plotNoTemp + theme(legend.position = 'none'),
                                 nrow = 1),
             nrow = 2, heights = c(1, 6))
dev.off()


