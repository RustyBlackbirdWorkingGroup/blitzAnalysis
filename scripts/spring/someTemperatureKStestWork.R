pnoSpring <- swdSpring %>%
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

tempPlot <- pnoSpring %>%
  dplyr::select(samplingPeriod, flockSize, tmin) %>%
  ggplot(
    aes(x = tmin)
  ) + 
  facet_grid(flockSize~samplingPeriod) +
  geom_density(aes(y = ..scaled..,fill = 'darkgray',  alpha = 0.8)) + 
  scale_fill_manual(values = cbPallete) + theme_bw() + 
  ylab('Scaled density') +
  xlab(expression(Average~minimum~temperature~(degree ~ C))) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15)
  ) + theme(legend.position="none")

png(filename = "C:/Users/Brian/Desktop/gits/blitzAnalysis/outPlots/tminSpring.png", 
    width = 8.5, height = 5.5, units = 'in', res = 300)
tempPlot
dev.off()

pnoSpring1 <- swdSpring %>%
  filter(count > 0) %>%
  mutate(flockSize = ifelse(
    count < 20,'Small',
    ifelse(count < 100, 'Medium',
           'Large')
  ))

pnoSpring1 %>%
  group_by(flockSize, samplingPeriod) %>%
  summarize(tmin = mean(tmin),
            se = se(tmin))

ks.test(pnoSpring1 %>%
          filter(samplingPeriod == 5) %>%
          filter(flockSize == 'Small') %>%
          .$tmin,
        pnoSpring1 %>%
          filter(samplingPeriod == 5) %>%
          filter(flockSize == 'Medium') %>%
          .$tmin
        )


ks.test(pnoSpring1 %>%
          filter(samplingPeriod == 5) %>%
          filter(flockSize == 'Small') %>%
          .$tmin,
        pnoSpring1 %>%
          filter(samplingPeriod == 5) %>%
          filter(flockSize == 'Large') %>%
          .$tmin
)

ks.test(pnoSpring1 %>%
          filter(samplingPeriod == 5) %>%
          filter(flockSize == 'Medium') %>%
          .$tmin,
        pnoSpring1 %>%
          filter(samplingPeriod == 5) %>%
          filter(flockSize == 'Large') %>%
          .$tmin
)


ks.test(pnoSpring1 %>%
          filter(samplingPeriod == 2) %>%
          filter(flockSize == 'Small') %>%
          .$tmin,
        pnoSpring1 %>%
          filter(samplingPeriod == 4) %>%
          filter(flockSize == 'Small') %>%
          .$tmin
)

ks.test(pnoSpring1 %>%
          #           filter(samplingPeriod == 2) %>%
          filter(flockSize == 'Small') %>%
          .$rowcrop,
        pnoSpring1 %>%
          # filter(samplingPeriod == 4) %>%
          filter(flockSize == 'Medium') %>%
          .$rowcrop
)



ks.test(pnoSpring1 %>%
#           filter(samplingPeriod == 2) %>%
          filter(flockSize == 'Small') %>%
          .$rowcrop,
        pnoSpring1 %>%
          # filter(samplingPeriod == 4) %>%
          filter(flockSize == 'Large') %>%
          .$rowcrop
)

ks.test(pnoSpring1 %>%
          #           filter(samplingPeriod == 2) %>%
          filter(flockSize == 'Medium') %>%
          .$rowcrop,
        pnoSpring1 %>%
          # filter(samplingPeriod == 4) %>%
          filter(flockSize == 'Large') %>%
          .$rowcrop
)

ksfun <- function(envVar, periodValue){
  compareFrame <- data.frame(
    sp1 = c('Small', 'Small', 'Medium'),
    sp2 = c('Medium', 'Large', 'Large')
  )
  for(i in 1:nrow(compareFrame)){
    env1 <- pnoSpring1 %>%
      filter(flockSize == compareFrame[i,1],
             samplingPeriod == periodValue) %>%
      select_(envVar) %>%
      as.matrix %>% as.vector
    env2 <- pnoSpring1 %>%
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

ksListPeriod <- vector('list', length = 5)

for(i in 1:5){
  ksListPeriod[[i]] <- ksfun('pasture', i) %>%
    mutate(samplingPeriod = i,
           p = round(p, 3))
}

bind_rows(ksListPeriod) %>%
  filter(p < .06)
