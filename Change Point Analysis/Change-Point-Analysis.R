#################################
# Non-breeding movements of Tufted Puffins in the eastern North Pacific Ocean: 
# temporal and spatial patterns
#
# Residency analysis using a changepoint model
#
# This code uses geolocator tag "BT132" as example. 
# 
# used: R version 4.1.3
####################################



#-----

## Install and load packages
#install.packages("GeoLight", "devtools")
#devtools::install_github("SLisovski/GeoLight")

library("GeoLight")
library("devtools")

#-----


# Set working directory to folder containing raw .lux file
  # Files needed for code to run:
      # .lux (light data)
setwd("C:/GeolocatorData")
getwd()
list.files()


#-----

# Format data for use in changelight analysis

# EXAMPLE TAG: BT132
# Load and transform lux file to format suitable for GeoLight
BT132 <- luxTrans("BT132_13Jul20_214708driftadj.lux")

# Subset to deployment & retrieval dates GMT
BT132_c1 <- BT132[8429:106061,]  #BT132: 7/19 - 6/22/2020

# Calculate twilight events (sunrise/sunset) from light intensity measurements over time 
BT132_c2 <- twilightCalc(BT132_c1$datetime, light = BT132_c1$light,LightThreshold = TRUE,
                         preSelection = TRUE, maxLight = NULL, ask = FALSE, 
                         allTwilights = TRUE)

# Format the "consecTwilights" column into a dataframe
BT132_c3 <- as.data.frame(BT132_c2$consecTwilights)

# Specify time zone for each time column
BT132_c3$tFirst <- as.POSIXct(BT132_c3$tFirst, tz = "GMT")
BT132_c3$tSecond <- as.POSIXct(BT132_c3$tSecond, tz = "GMT")


#-----

# Calculate residencies for each tag using change point analysis 

BT132_residency <- changeLight(BT132_c3$tFirst, BT132_c3$tSecond, 
                               BT132_c3$type, rise.prob = 0.15, 
                               set.prob = 0.15, plot = TRUE, 
                               summary = TRUE)
capture.output(BT132_residency$migTable, file = "BT132_residency.txt")





