library(ggplot2)
library(maps)
library(maptools)
library(rgeos)

gClip <- function(shp, bb){
  if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else b_poly <- as(extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}

states<-maps::map("state", col="transparent", plot=FALSE, fill = TRUE)



# Get us shape file and clip it to the raster output extent:

us <- map2SpatialPolygons(states, IDs = states$names,
                             proj4string=CRS("+proj=longlat +datum=WGS84")) %>%
  gSimplify(tol = 0.0001) %>%
  gBuffer(byid = TRUE, width = 0) %>%
  gClip(bbox(extent(small2009)))

# Aggregate then mask rasters to shapefiles and put in format that ggplot can read:

prepRasterToPlot = function(r, flockSize){
  r %>%
    # aggregate(2) %>%
    mask(us) %>%
    rasterToPoints %>%
    data.frame %>%
    mutate(flockClass = flockSize)
}

# Prepare rasters for plotting in ggplot:

rastersForPlotting <- bind_rows(
  prepRasterToPlot(small2009, 'Small'),
  prepRasterToPlot(medium2009, 'Medium'),
  prepRasterToPlot(large2009, 'Large')) %>%
  mutate(flockClass = factor(flockClass, levels = c('Small', 'Medium', 'Large')),
         layer = ifelse(layer < 0.00, NA, layer))

# Plot

rastersForPlotting %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = layer)) +
  # geom_polygon(aes(data = usGGshape, x = long, y = lat, group = group), fill = NA, color = 'black')
  facet_grid(. ~ flockClass) +
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

outDirectory <- 'C:/Users/Brian/Desktop/gits/blitzAnalysis/outPlots/'

ggsave(paste0(outDirectory, 'winterSuitabilityMaps2009.png'), 
       width = 6.5, height = 2.3, units = 'in')

# 
# bind_rows((sfLogisticMean/indLogisticMean) %>% aggregate(2) %>% rasterToPoints %>% data.frame %>%
#                  mutate(d = 'IndSF'),
#                (lfLogisticMean/indLogisticMean) %>% aggregate(2) %>% rasterToPoints %>% data.frame %>%
#                  mutate(d = 'IndLF')) %>%
#   mutate(d = factor(d, levels = c('IndSF', 'IndLF'), 
#                     labels = c('Small flocks - individuals)', 'Large flocks - individuals')),
#          layer = ifelse(layer < 0.0001 & layer > - 0.0001, NA, layer)) %>%
#   ggplot(aes(x, y)) + 
#   geom_raster(aes(fill = layer)) +
#   facet_grid(. ~ d) +
#   theme_bw() + 
#   labs(x = 'Latitude', y = 'Longitude') +
#   theme(axis.title.y = element_text(vjust = 1.5, size = rel(1.75)),
#         axis.title.x = element_text(vjust = -0.25, size = rel(1.75)),
#         axis.text = element_text(size = rel(1.25)),
#         strip.text = element_text(size = rel(2), face = 'bold')) + 
#   scale_fill_gradientn(colours=c('navyblue', 'darkblue', 'blue' ,'blue', 'gray50', 'yellow', 'red', 'darkred'), 
#                          na.value = 'gray20')

# e <- extent(c(-93,-85,28.49803,37))
# 
# bind_rows((crop(sfLogisticMean,e)) %>%
#             rasterToPoints %>% data.frame %>%
#             mutate(d = 'sf'),
#           (crop(lfLogisticMean,e)) %>%
#             rasterToPoints %>% data.frame %>%
#             mutate(d = 'lf')) %>%
#   filter(d!= 'LFSF') %>%
#   mutate(d = factor(d, levels = c('sf', 'lf'), labels = c('Small flock', 'Large flock')),
#          layer = ifelse(layer < 0.0001 & layer > - 0.0001, NA, layer)) %>%
#   ggplot(aes(x, y)) + 
#   geom_raster(aes(fill = layer)) +
#   labs(x = 'Latitude', y = 'Longitude') +
#   facet_grid(. ~ d) +
#   theme_bw() + 
#   theme(axis.title.y = element_text(vjust = 1.5, size = rel(1.75)),
#         axis.title.x = element_text(vjust = -0.25, size = rel(1.75)),
#         axis.text = element_text(size = rel(1.25)),
#         strip.text = element_text(size = rel(2), face = 'bold')) + 
#   scale_fill_gradientn(colours=c('navyblue', 'darkblue', 'blue' ,'blue', 'gray50', 'yellow', 'red', 'darkred'), 
#                                     na.value = 'gray20')







# setwd('C:/Users/Brian/Desktop/gits/RUBL')


