#===================================================================================================*
# ---- SET-UP ----
#===================================================================================================*

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
                      'lubridate','raster', 'dismo', 'ggplot2',
                      'maps', 'maptools', 'rgeos'))

# Set out directory (example, do not run):

outDirectory <- 'C:/Users/Brian/Desktop/gits/blitzAnalysis/outPlots/'


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

#===================================================================================================*
# ---- PLOT DATA ----
#===================================================================================================*

rastersForPlotting %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = layer)) +
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

ggsave(paste0(outDirectory, 'winterSuitabilityMaps2009.png'), 
       width = 6.5, height = 2.3, units = 'in')