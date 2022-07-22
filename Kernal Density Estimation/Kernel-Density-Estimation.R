#################################
# Non-breeding movements of Tufted Puffins in the eastern North Pacific Ocean: 
# temporal and spatial patterns
#
# Seasonal kernel density estimation and mapping
# 
# This code uses geolocator tag "BT132" as example. 
# 
# used: R version 4.1.3
####################################

#-----


## Load packages

# Kernel density estimation
require(adehabitatHR)
require(CircStats)
require(sp)
require(rgdal)
require(rgeos)
library(rgdal)
library(sf)


# Mapping
library(ggplot2)
library(maps)
library(mapdata)
library(mapproj)

#-----

# Set working directory 
  # Files needed for code to run:
    # .csv file: output from probGLS analysis with twice-daily "most probable locations"
    # with additional column "BroadSiteType" (FALL, WINTER, SPRING) for residency periods 
    # derived from change point analysis

#set working directory
setwd("C:/GeolocatorData")
getwd()
list.files()


#-----

# Load and format data
  # Example tag: BT132
BT132<- read.csv("BT132_M_19-20_withSitesRefined.csv",header=T)

#Plotting will only recognize long in the 360 degree format. This 
  #creates a new column in each file with the 360 version
BT132$lon2 <- ifelse(BT132$lon < 0, BT132$lon + 360, BT132$lon)


# Subset by season, if desired. Can also combine (rbind) multiple tracks together 
#(e.g., males, females, by year, etc), however you want to calculate the utilization distributions. 
# In this example, we'll subset to the "WINTER" season, as defined by change-point analysis. 

BT132.winter <- BT132[BT132$BroadSiteType == "WINTER1",]


# Format desired KD data into a spatial data frame 
BT132.winter.sp <- SpatialPointsDataFrame(coords=as.data.frame(cbind(BT132.winter$lon,BT132.winter$lat)),
                                       data=BT132.winter, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# Transform spatial data frame
BT132.winter.sp.p <- spTransform(BT132.winter.sp, CRS("+init=epsg:3571"))


#-----

# Kernel density/utilization distribution estimation

# Calculate the kernel density distribution using adehabitatR package. 
kud.BT132winter <- kernelUD(BT132.winter.sp.p, h = "href", grid = 1000)


#Calculate the % utilization distributions of interest (e.g., 25%, 50%, 75%, 95%)

ud25.BT132winter <- getverticeshr(kud.BT132winter, percent=25, standardize=TRUE)
ud50.BT132winter <- getverticeshr(kud.BT132winter, percent=50, standardize=TRUE)
ud75.BT132winter <- getverticeshr(kud.BT132winter, percent=75, standardize=TRUE)
ud95.BT132winter <- getverticeshr(kud.BT132winter, percent=95, standardize=TRUE)


#-----

# Map seasonal utilization distributions

# Transform polygons
ud25.BT132winter <- spTransform(ud25.BT132winter, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))  
ud50.BT132winter <- spTransform(ud50.BT132winter, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))  
ud75.BT132winter <- spTransform(ud75.BT132winter, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))  
ud95.BT132winter <- spTransform(ud95.BT132winter, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))  


# Create Map

# Define map area (this is far beyond the needed extent)
map1<-map_data(map("world2", fill=T, xlim=c(110,240), ylim=c(25,75), col="grey46"))

# Middleton Island Lat/lon - location of tagging
mdo <- c(-146.3414,59.40586014)


kdmap.BT132winter <- ggplot()+
  geom_polygon(data = ud95.BT132winter, aes(x=long, y = lat, group = group), fill = "lightblue1", show.legend = T)+
  geom_polygon(data = ud75.BT132winter, aes(x=long, y = lat, group = group), fill = "lightskyblue", show.legend = T)+
  geom_polygon(data = ud50.BT132winter, aes(x=long, y = lat, group = group), fill = "steelblue2", show.legend = T)+
  geom_polygon(data = ud25.BT132winter, aes(x=long, y = lat, group = group), fill = "blue", show.legend = T)+
  geom_polygon(data=map1, aes(x=long, y=lat, group=group), color = "grey84", fill = "snow2")+
  coord_map("lambert",40,75,xlim = c(200,230),ylim=c(45,65))+
  scale_x_continuous(breaks=c(200, 210, 220, 230),
                     labels=c("-160", "-150", "-140", "-130"))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()
kdmap.BT132winter

