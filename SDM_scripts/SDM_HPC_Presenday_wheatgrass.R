# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

#SDM WORLDCLIM 1970-2000 30 seconds
library(dismo)
library(raster)
library(maptools)
library(rasterVis)
data(wrld_simpl)
library(readr)

################
# INPUT: Lat & Long Data
################

#LOCAL
#datapoints <- read.csv("~/R/Intermediate_Wheatgrass_collab/IWG_lat_long.csv", header=TRUE, sep=",")

#HPC file path
datapoints <- read.csv("/home/ahmccorm/kantar_koastore/anna/Intermediate_wheatgrass_collab/Present_day/IWG_lat_long.csv", header=TRUE, sep=",")
head(datapoints)

# Remove rows with NA values in Longitude and Latitude
datapoints <- datapoints[complete.cases(datapoints$Longitude), ]
datapoints <- datapoints[complete.cases(datapoints$Latitude), ]

# Summary
summary(datapoints)

# Identify duplicates
dups <- duplicated(datapoints[, c("Longitude", "Latitude")])
sum(dups)

# Keep only non-duplicated records
datapoints <- datapoints[!dups, ]

# Keep only the Longitude and Latitude columns
datapoints <- datapoints[, c("Longitude", "Latitude")]

#have one rogue Canadian value to remove
datapoints <- datapoints[datapoints$Longitude > -25, ]

# View first few rows
head(datapoints)

# Plot world map with datapoints
library(maptools)
data(wrld_simpl)
plot(wrld_simpl, axes=TRUE, col="lightgreen")
box()
points(datapoints$Longitude, datapoints$Latitude, col='blue', pch=20, cex=0.75)


######################################################
# Load climate data (all variables) wc2.0_30s_bio
######################################################
#HPC filepath
files <- list.files(path="~/kantar_koastore/anna/Barley_collab/barley_parental/SDM/wc2.0_30s_bio/",
                    pattern='tif', full.names=TRUE)

#LOCAL
#files <- list.files(path="/Users/annamccormick/R/Data/climate_data/wc2.0_30s_bio/", pattern='tif', full.names=TRUE)


predictors <- stack(files)

# Plot the first layer of the RasterStack
plot(predictors, 1)
plot(wrld_simpl, add=TRUE)
points(datapoints$Longitude, datapoints$Latitude, col='blue')

# Extract climate values for the datapoints
presvals <- extract(predictors, datapoints)

# Generate random background points
set.seed(0)
backgr <- randomPoints(predictors, 500)

# Extract background climate values
absvals <- extract(predictors, backgr)

# Presence/absence labels
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))

# Create dataset for SDM
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
head(sdmdata)

# Drop 'biome' layer from predictors
pred_nf <- dropLayer(predictors, 'biome')

# Perform k-fold cross-validation (split into 3 groups)
group <- kfold(datapoints, 3)

# Training and testing data
pres_train <- datapoints[group != 1, ]
pres_test <- datapoints[group == 1, ]

# Generate background points
backg <- randomPoints(pred_nf , n=1000,  extf = 1.25)
colnames(backg) = c('Longitude', 'Latitude')

group <- kfold(backg, 3)
backg_train <- backg[group != 1, ]
backg_test <- backg[group == 1, ]

# Final visualization to check correct
r = raster(pred_nf, 1)
plot(!is.na(r), col=c('white', 'light grey'), legend=FALSE)
plot(wrld_simpl, add=TRUE, col='red', lwd=2)
points(backg_train, pch='-', cex=0.5, col='yellow')
points(backg_test, pch='-', cex=0.5, col='black')
points(pres_train, pch= '+', col='green')
points(pres_test, pch='+', col='blue')

################
#Sanity checks for splits
################
# Count number of training and testing points
n_train <- nrow(pres_train)
n_test <- nrow(pres_test)

# Print counts
print(paste("Number of training points:", n_train))
print(paste("Number of testing points:", n_test))

# Count number of background training and testing points
n_backg_train <- nrow(backg_train)
n_backg_test <- nrow(backg_test)

# Print counts
print(paste("Number of background training points:", n_backg_train))
print(paste("Number of background testing points:", n_backg_test))

head(pres_train)
head(backg_train)


################
#Maxent
################
library(rJava)
install.packages("snow")

jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

#Variable Contribution
if (file.exists(jar)) {
  xm <- maxent(predictors, pres_train)
  plot(xm)
} else {
  cat('cannot run this example because maxent is not available')
  plot(1)
}

#World Suitability Map
if (file.exists(jar)) {
  response(xm)
} else {
  cat('cannot run this example because maxent is not available')
  plot(1)
}




if (file.exists(jar)) {
  xme <- evaluate(pres_test, backg_test, xm, predictors)
  px <- predict(predictors, xm, progress='')
  trxm <- threshold(xme, 'spec_sens')

  par(mfrow=c(1,2))
  plot(px, main='Maxent, raw values')#xlim=c(-160,-55), ylim=c(0,80), main='Maxent, raw values')
  plot(wrld_simpl, add=TRUE, border='dark grey')

  plot(px > trxm, main='presence/absence')#xlim=c(-160,-55), ylim=c(0,80), main='presence/absence')
  plot(wrld_simpl, add=TRUE, border='dark grey')
  points(pres_train, pch='+')
} else {
  plot(1)
}


writeRaster(px,'/home/ahmccorm/kantar_koastore/anna/Intermediate_wheatgrass_collab/Present_day/IWG_SDM_presentday.tif',options=c('TFW=YES'))

# Plot ROC curve (AUC)
par(mfrow=c(1, 1)) # Resetting to default to ensure the plot is not side by side
plot(xme, 'ROC', cex.lab = 1, cex.axis =1, cex.main =2)

# Now the next plot will be separate
plot(xm, cex.lab = 2, cex.axis = 4, cex.main = 2, cex.sub = 3) # Variable contribution

# Save the workspace to a .RData file
save.image(file="/home/ahmccorm/kantar_koastore/anna/Intermediate_wheatgrass_collab/Present_day/my_workspace_Present.RData")
