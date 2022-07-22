################################
# Non-breeding movements of Tufted Puffins in the eastern North Pacific Ocean: 
# temporal and spatial patterns
# 
# Processing of raw geolocator data 
#
# used: R version 4.1.3 
################################



# Use this code to process raw geolocator light data, estimate most probable 
# twice day locations, and display most probable locations on map.  

# This code uses geolocator tag "BT132" as example. 


#-----

## Twilight determination and light data processing

# Install and load packages
install.packages("devtools")
library(devtools)

devtools::install_github("SWotherspoon/BAStag") 
library(BAStag)

# Format light data for probGLS analysis
devtools::install_github("SLisovski/TwGeos")
library(TwGeos)

# probGLS analysis
library(plyr)
library(rgdal)
library(stringr)
library(sp)

install_github("benjamin-merkel/probGLS")
library(probGLS)
library(ncdf4)
library(tidyverse)
library(lubridate)

# Display most probable locations on map
library(rgeos)


#------



# Set working directory to folder containing raw geolocator files
  # Files needed for code to run:
      # .lux (light data), .deg (activity data), .sst (sea surface temperature data)
setwd("C:/Users/GeolocatorData")
getwd()
list.files()


###-----
### Twilight determination for Migrate Technology Tags using BAStag and TwGeos

# Define parameters
zenith <- 96
threshold <-1.5 # Recommended threshold for Integeo Migrate Technology tags (Rakhimberdiev et al. 2016)
  # below is night, above is day 
dt <- 600 
offset <- 4

# Define calibration location 
cal.loc <- c(-145.9748, 60.39198)
locations <- matrix(c(cal.loc),1,2,byrow=T)

# Read in drift-adjusted raw light file (.lux)
  # Example tag: BT132
d.lig<-readMTlux("BT132_13Jul20_214708driftadj.lux", skip=20)
lightImage(d.lig,offset=offset)

# Extract calibration segment
  #BT132 calibration period: 6/19/2019 - 6/27/2019
d.calib <- subset(d.lig,d.lig$Date > min(d.lig$Date)-0*24*60*60 & d.lig$Date < min(d.lig$Date)+8*24*60*60)
lightImage(d.calib,offset=offset)
thresholdCalibrate(d.calib,locations[1,1],locations[1,2],
                   xlim=c(92,98),ylim=c(0,15),pch=16)
abline(h=threshold,v=zenith) 
#-----


## Interactive processing. 
## The process has 4 stages:
##
## 1. The subset of data is selected (left click=start, right click = end)
## 2. Twilight search (left click = search point, right click = prevent search)
## 3. Insert missing twilights (left click = sunset, right click = sunrise)
## 4. Edit individual twilight (left click = select time, right click = delete)
##

## Process light data
twl <- preprocessLight(d.lig,threshold,offset=offset,zenith=zenith,
                       fixed=locations)
head(twl)

## Investigate twilights less than 1 hrs apart
which(diff(as.numeric(twl$Twilight)) < 60*60)

# Investigate rise/set out of sequence
which(diff(twl$Rise)==0) 
twl <- twl[!twl$Deleted,]

## Store twilights
save(twl,file="BT132_twlraw.RData")


## Use TwGeos to format for ProbGLS analysis
twl<-export2GeoLight(twl)



# save raw twilights as CSV file. 
load("BT132_twlraw.RData")
twl
write.csv(twl, file = "BT132_twlraw.csv")
twl<-export2GeoLight(twl)


#-----


## ProbGLS analysis



# Model definition

# Set wet dry data resolution (in sec)
wdr <- 6 

# Load twilight events derived above
trn<-twl

# Load drift-adjusted activity file (wet/dry data) and format date/time
act <- read.table('BT132_13Jul20_214708driftadj.deg',sep="\t",skip=19,header=T) 
act$dtime <- dmy_hms(act$DD.MM.YYYY.HH.MM.SS)
act$dtime <- as.POSIXct(act$dtime, tz = "UTC")

# Rename activity column - wet dry data column must be called 'wetdry'
act$wetdry<- act$wets0.50 
act <- act[!is.na(act$dtime),]

# Load tag SST file and format date/time
td <- read.table('BT132_13Jul20_214708.sst',sep="\t",skip=19,header=T) 
td$dtime <- as.POSIXct(strptime(td[,1], format = "%d/%m/%Y %H:%M:%S"),tz='UTC')

sen <- sst_deduction(datetime = td$dtime, temp = td[,4], temp.range = c(-2,60))
sen <- sen[sen$SST.remove==F,]

# run algorithm
tw <- twilight_error_estimation()

# Calc stdev(td)
td
head(td)
sd <- sd(td$wet.max.C..wet.mean.C.)

# Define algorithm parameters: 
pr <- prob_algorithm(particle.number  = 5000,
                     iteration.number            = 500,
                     trn                         = trn, #[strftime(trn$tFirst,"%m") %in% c(8:12,1:4),], 
                     sensor                      = sen[sen$SST.remove==F,],
                     act                         = act, #act[act$X=="wet",], 
                     loess.quartile              = 10,
                     tagging.location            = c(-146.3420706, 59.40751557), 
                     tagging.date                = as.Date("2019-07-19"), #GMT
                     retrieval.date              = as.Date("2020-06-22"), #GMT
                     wetdry.resolution           = wdr,
                     sunrise.sd                  = tw,
                     sunset.sd                   = tw,
                     range.solar                 = c(-7, -1),
                     speed.wet                   = c(1.2 , 0.34, 1.8), #Lapsansky et al. 2020, 
                     speed.dry                   = c(16.5, 2.5, 22), #Spear & Ainley 1997
                     sst.sd                      = sd,       
                     max.sst.diff                = 3,  #Difference allowed between tag & satellite temp readings        
                     days.around.spring.equinox  = c(21,14),   
                     days.around.fall.equinox    = c(14,21), 
                     ice.conc.cutoff             = 0,
                     boundary.box                = c(121,-116,30,90),
                     land.mask                   = T,
                     med.sea                     = T,        
                     black.sea                   = T,        
                     baltic.sea                  = T,      
                     caspian.sea                 = T,    
                     east.west.comp              = F,
                     NOAA.OI.location            = "C:/Users/Environmental Variables for probGLS") 

alarm()


# Plot of latitude, longitude, and SST over time
plot_timeline(pr,degElevation = NULL)


# Display most probable locations on map
cs<-readOGR("C:/Users/Environmental Variables for probGLS", "ne_50m_land")
proj.stereo <- "+proj=stere +lat_0=90 +lat_ts=75 +lon_0=-170 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

pr1 <-pr
pr1[[1]]<- spTransform(pr[[1]],CRS(proj.stereo))
pr1[[2]]<- spTransform(pr[[2]],CRS(proj.stereo))
cs2<- spTransform(crop(cs,extent(-180,180,0,90)),CRS(proj.stereo))

opar <- par(mfrow=c(1,1),mar=c(4,4,0,0))
plot(pr1[[1]],col="white")
for(s in 1:length(unique(pr[[2]]$step))){
  plot(pr1[[1]][pr1[[1]]$step==unique(pr1[[1]]$step)[s],],col=colorRampPalette(c('grey90','grey50'))(nrow(pr[[2]]))[s],
       add=T,pch=19,cex=0.3)
}
plot(cs2,add=T,border=8)
mm2 <-pr1[[2]]
lines(mm2$lon,mm2$lat,col=1)
points(mm2$lon,mm2$lat,cex=1,pch=21,bg=colorRampPalette(c('yellow','darkred'))(nrow(pr[[2]])))
mm3 <- mm2[is.na(mm2$median.sun.elev),]
points(mm3$lon,mm3$lat,cex=0.7,pch=3)

par(opar) 


# Write data to file
write.table(pr$`most probable track`, file="BT132_mostprob.csv", col.names=NA, sep=",")







