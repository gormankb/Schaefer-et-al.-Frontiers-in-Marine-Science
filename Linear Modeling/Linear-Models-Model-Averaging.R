################################
# Non-breeding movements of Tufted Puffins in the eastern North Pacific Ocean: 
# temporal and spatial patterns
# 
# Linear models and model averaging to test for differences in migratory patterns 
# of geo-tagged Tufted Puffins
#
# used: R version 4.1.3 
################################



# Testing for differences in multiple response variables
# using LMs, with sex, year, body size, and body condition as explanatory variables.  


#-----


# Load packages
library(AICcmodavg)
library(MuMIn)
library(sjPlot)


#-----


# Set working directory 
  # File needed for code to run: .csv file containing full migration summaries
  # for each tag


setwd("C:/MigrationSummaries")
getwd()
list.files()


# Load and format data

df <- read.csv("MASTER_TUPU_MigrationSummaries_2018-2021_forR.csv",header=T)
dim(df)
head(df)

# Remove 2020-2021 tracks (n = 2) (Year == 2020/2021)
nrow(df[df$Year == "2020/2021",]) #2
data<- df[!(df$Year == "2020/2021"),]
dim(data)


# Define variables to be used in models
  # Response variables
mdo.departure <- data$MDO_Departure_DOY
fall.dist <- data$CumulativeDistToWinter_m
weekly.fall.dist <- data$Weekly_Fall_Dist_km
mdo2winter.dist <- data$StraightDistMDOtoWinter
days2winter <- data$Days_to_Winter_Site
num.fall.sites <- data$Num_Fall_Sites
winter.lon <- data$Mean_Winter_Lon
winter.lat <- data$Mean_Winter_Lat
winter.arrival <- data$Winter_Arrival_DOY
winter.duration <- data$Winter_Num_Days
winter.departure <- data$Winter_Departure_DOY
weekly.spring.dist <- data$Weekly_Spring_Dist_km
spring.dist <- data$CumulativeDistToMDO_m
num.spring.sites <- data$Num_Spring_Sites
days2mdo <- data$Days_to_Breeding_Site
mdo.arrival <- data$MDO_Arrival_DOY

  # Explanatory variables
sex <- as.factor(data$Sex)
year <- as.factor(data$Year)
condition <- data$BodyCondition

# create dataframe with terms of interest
df2 <- cbind(mdo.departure, fall.dist, weekly.fall.dist, mdo2winter.dist, days2winter, num.fall.sites, 
             winter.lon, winter.lat, winter.arrival, winter.duration, winter.departure, 
             spring.dist, weekly.spring.dist, num.spring.sites, days2mdo, mdo.arrival, sex, year, 
             condition)
df3 <- as.data.frame(df2)


#-----

# Examine variable distributions
hist(mdo.departure)
hist(fall.dist)
hist(weekly.fall.dist)
hist(mdo2winter.dist)
hist(days2winter)
hist(num.fall.sites)
hist(winter.lon)
hist(winter.lat)
hist(winter.arrival)
hist(winter.duration)
hist(winter.departure)
hist(spring.dist)
hist(weekly.spring.dist)
hist(num.spring.sites)
hist(days2mdo)
hist(mdo.arrival)

hist(condition)


#------

# Linear Models: 

# Candidate Model Set: 
  # 1. null
  # 2. ~ sex
  # 3. ~ year
  # 4. ~ condition
  # 5. ~ sex + year
  # 6. ~ sex + condition
  # 7. ~ year + condition
  # 8. ~ sex + year + condition
  # 9. ~ sex + condition + sex * condition
  # 10. ~ year + condition + year * condition
  # 11. ~sex + year + sex * year


# Example for one response variable: Departure date from MDO (breeding grounds)

### Departure date from MDO
mdo.departure1<- lm(mdo.departure ~ 1)

mdo.departure2<- lm(mdo.departure ~ sex)

mdo.departure3<- lm(mdo.departure ~ year)

mdo.departure4<- lm(mdo.departure ~ condition)

mdo.departure5<- lm(mdo.departure ~ sex + year)

mdo.departure6<- lm(mdo.departure ~ sex + condition)

mdo.departure7<- lm(mdo.departure ~ year + condition)

mdo.departure8<- lm(mdo.departure ~ sex + year + condition)

mdo.departure9<- lm(mdo.departure ~ sex*condition)

mdo.departure10<- lm(mdo.departure ~ year*condition)

mdo.departure11<- lm(mdo.departure ~ sex*year)

mdo.departure.models <- list(mdo.departure1, mdo.departure2, mdo.departure3, 
                             mdo.departure4, mdo.departure5, mdo.departure6, 
                             mdo.departure7, mdo.departure8, mdo.departure9, 
                             mdo.departure10, mdo.departure11)

# AICc table
mdo.departure.aictab <- aictab(mdo.departure.models, second.ord = TRUE, sort = TRUE)
mdo.departure.aictab

# Model selection table
mdo.departure.output <- model.sel(mdo.departure.models)
mdo.departure.output

# Model averaging using all candidate models
mdo.departure.modavg <- model.avg(mdo.departure.output, revised.var = TRUE)
summary(mdo.departure.modavg)

# 95% confidence intervals for explanatory variables
confint(mdo.departure.modavg)

# sum of weights for explanatory variables
sw(mdo.departure.modavg)


