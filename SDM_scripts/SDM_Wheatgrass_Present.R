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
datapoints <- read.csv("~/R/Intermediate_Wheatgrass_collab/IWG_lat_long.csv", header=TRUE, sep=",")

#HPC file path
#datapoints <- read.csv("/home/ahmccorm/kantar_koastore/anna/Intermediate_wheatgrass_collab/Present_day/IWG_lat_long.csv", header=TRUE, sep=",")
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
#files <- list.files(path="~/kantar_koastore/anna/Barley_collab/barley_parental/SDM/wc2.0_30s_bio/", pattern='tif', full.names=TRUE)

#LOCAL
files <- list.files(path="/Users/annamccormick/R/Data/climate_data/wc2.0_30s_bio/", pattern='tif', full.names=TRUE)

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


####################################################################################################################################


#################################
#Plots - whole world
#################################
library(ggplot2)
library(viridis)
library(raster)
species_distribution <- raster("~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_presentday.tif")

species_distribution_reduced <- aggregate(species_distribution, fact=5, fun=mean) # 'fact' is the aggregation factor

# Convert to data frame
df_reduced <- as.data.frame(rasterToPoints(species_distribution_reduced))
names(df_reduced) <- c("Longitude", "Latitude", "Value")

p5<- ggplot(data = df_reduced, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Species Distribution") +
  labs(title = "Present 1970-2000") +
  theme_minimal()

ggsave("~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_present.pdf", plot = p5, width = 12, height = 6.5, device = "pdf")

#################################
# Plot - Minnesota Only
#################################
library(ggplot2)
library(viridis)
library(raster)
library(sf)

# Load the raster (species distribution)
species_distribution <- raster("~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_presentday.tif")

# Load the U.S. shapefile and extract Minnesota
us_shape <- st_read("~/R/states_21basic/states.shp")
minnesota_shape <- us_shape[us_shape$STATE_NAME == "Minnesota", ]

# Convert to Spatial object for raster operations
minnesota_shape_sp <- as(minnesota_shape, "Spatial")

# Crop and mask the raster to Minnesota's extent
species_mn_crop <- crop(species_distribution, extent(minnesota_shape_sp))
species_mn_masked <- mask(species_mn_crop, minnesota_shape_sp)

# Save the cropped raster as a TIFF file
writeRaster(species_mn_masked, 
            "~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_present_MN.tif", 
            format = "GTiff")

# Reduce resolution for plotting
species_mn_reduced <- aggregate(species_mn_masked, fact=5, fun=mean)

# Convert to data frame for ggplot
df_mn <- as.data.frame(rasterToPoints(species_mn_reduced))
names(df_mn) <- c("Longitude", "Latitude", "Value")

# Plot for Minnesota
p_mn <- ggplot(data = df_mn, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Species Distribution") +
  labs(title = "Intermediate Wheatgrass - Minnesota (1970-2000)") +
  theme_minimal()

# Save the plot as a PDF
ggsave("~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_present_MN.pdf", 
       plot = p_mn, 
       width = 8, height = 6.5, device = "pdf")

# Display the plot
p_mn


#################################
# Plot - USA Only
#################################
library(ggplot2)
library(viridis)
library(raster)
library(sf)

# Load the raster (species distribution)
species_distribution <- raster("~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_presentday.tif")

# Load the U.S. shapefile
us_shape <- st_read("~/R/states_21basic/states.shp")

# Convert to Spatial object for raster operations
us_shape_sp <- as(us_shape, "Spatial")

# Crop and mask the raster to the USA's extent
species_us_crop <- crop(species_distribution, extent(us_shape_sp))
species_us_masked <- mask(species_us_crop, us_shape_sp)

# Save the cropped raster as a TIFF file
writeRaster(species_us_masked, 
            "~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_present_USA.tif", 
            format = "GTiff")

# Reduce resolution for plotting
species_us_reduced <- aggregate(species_us_masked, fact=5, fun=mean)

# Convert to data frame for ggplot
df_us <- as.data.frame(rasterToPoints(species_us_reduced))
names(df_us) <- c("Longitude", "Latitude", "Value")

# Plot for the USA
p_us <- ggplot(data = df_us, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Species Distribution") +
  labs(title = "Intermediate Wheatgrass - USA (1970-2000)") +
  theme_minimal()

# Save the plot as a PDF
ggsave("~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_present_USA.pdf", 
       plot = p_us, 
       width = 12, height = 6.5, device = "pdf")

# Display the plot
p_us



################################################################################################################################################################
#FUTURE
#####################################################################################################################
# Bring 13 CMIP6 models together
######################################
library(raster)

# List of raster inputs
raster_files <- list.files(path = "~/R/Intermediate_Wheatgrass_collab/SDM_Future_2050", 
                           pattern = "^Future_2050_.*\\.tif$", 
                           full.names = TRUE)

# Check the loaded files
print(raster_files)

# Load rasters as a stack
raster_stack <- stack(raster_files)
print(raster_stack)

# Compute the mean raster
mean_raster <- calc(raster_stack, mean, na.rm = TRUE)
# Plot the mean raster
#plot(mean_raster, main = "Mean SDM across Models")
writeRaster(mean_raster, "~/R/Intermediate_Wheatgrass_collab/IWG_13_model_SDM_future_2050.tif", format = "GTiff", overwrite = TRUE)


# Compute the standard deviation raster
#sd_raster <- calc(raster_stack, sd, na.rm = TRUE)
# Plot the standard deviation raster
#plot(sd_raster, main = "Standard Deviation of SDMs")
#writeRaster(sd_raster, "~/R/Intermediate_Wheatgrass_collab/IWG_SD_13_model_SDM_future_2050.tif", format = "GTiff", overwrite = TRUE)


#Plots - whole world
#################################
library(ggplot2)
library(viridis)
library(raster)
species_distribution <- raster("~/R/Intermediate_Wheatgrass_collab/IWG_13_model_SDM_future_2050.tif")
#species_distribution <- raster("~/R/barley_collab/Collab_Proposal_SDMs/Barley_SDM_ssp585_2070_merged.tif")

species_distribution_reduced <- aggregate(species_distribution, fact=5, fun=mean) # 'fact' is the aggregation factor

# Convert to data frame
df_reduced <- as.data.frame(rasterToPoints(species_distribution_reduced))
names(df_reduced) <- c("Longitude", "Latitude", "Value")

p5<- ggplot(data = df_reduced, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Species Distribution") +
  labs(title = "Future 2041-2060") +
  theme_minimal()

ggsave("~/R/Intermediate_Wheatgrass_collab/IWG_13_model_SDM_future_2050_world.pdf", plot = p5, width = 12, height = 6.5, device = "pdf")



#################################
# Plot - USA Only
#################################
library(ggplot2)
library(viridis)
library(raster)
library(sf)

# Load the raster (species distribution)
species_distribution <- raster("~/R/Intermediate_Wheatgrass_collab/IWG_13_model_SDM_future_2050.tif")

# Load the U.S. shapefile
us_shape <- st_read("~/R/states_21basic/states.shp")

# Convert to Spatial object for raster operations
us_shape_sp <- as(us_shape, "Spatial")

# Crop and mask the raster to the USA's extent
species_us_crop <- crop(species_distribution, extent(us_shape_sp))
species_us_masked <- mask(species_us_crop, us_shape_sp)

# Save the cropped raster as a TIFF file
writeRaster(species_us_masked, 
            "~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_Future_2050_USA.tif", 
            format = "GTiff")

# Reduce resolution for plotting
species_us_reduced <- aggregate(species_us_masked, fact=5, fun=mean)

# Convert to data frame for ggplot
df_us <- as.data.frame(rasterToPoints(species_us_reduced))
names(df_us) <- c("Longitude", "Latitude", "Value")

# Plot for the USA
p_us <- ggplot(data = df_us, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Species Distribution") +
  labs(title = "Intermediate Wheatgrass - USA (2041-2060)") +
  theme_minimal()

# Save the plot as a PDF
ggsave("~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_Future_2050_USA.pdf", 
       plot = p_us, 
       width = 12, height = 6.5, device = "pdf")

# Display the plot
p_us


#################################
# Plot - Minnesota Only
#################################
library(ggplot2)
library(viridis)
library(raster)
library(sf)

# Load the raster (species distribution)
species_distribution <- raster("~/R/Intermediate_Wheatgrass_collab/IWG_13_model_SDM_future_2050.tif")

# Load the U.S. shapefile and extract Minnesota
us_shape <- st_read("~/R/states_21basic/states.shp")
minnesota_shape <- us_shape[us_shape$STATE_NAME == "Minnesota", ]

# Convert to Spatial object for raster operations
minnesota_shape_sp <- as(minnesota_shape, "Spatial")

# Crop and mask the raster to Minnesota's extent
species_mn_crop <- crop(species_distribution, extent(minnesota_shape_sp))
species_mn_masked <- mask(species_mn_crop, minnesota_shape_sp)

# Save the cropped raster as a TIFF file
writeRaster(species_mn_masked, 
            "~/R/Intermediate_Wheatgrass_collab/SDM_Future_2050/IWG_SDM_Future_2050_MN.tif", 
            format = "GTiff")

# Reduce resolution for plotting
species_mn_reduced <- aggregate(species_mn_masked, fact=5, fun=mean)

# Convert to data frame for ggplot
df_mn <- as.data.frame(rasterToPoints(species_mn_reduced))
names(df_mn) <- c("Longitude", "Latitude", "Value")

# Plot for Minnesota
p_mn <- ggplot(data = df_mn, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_raster() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Species Distribution") +
  labs(title = "Intermediate Wheatgrass - Minnesota (2041-2060)") +
  theme_minimal()

# Save the plot as a PDF
ggsave("~/R/Intermediate_Wheatgrass_collab/SDM_Future_2050/IWG_SDM_Future_2050_MN.pdf", 
       plot = p_mn, 
       width = 8, height = 6.5, device = "pdf")

# Display the plot
p_mn


###############################################################################################################################################
##########################################
# SDM Blankets 0.2 for USA
##########################################

# Load necessary libraries
library(raster)
library(terra)
library(sf)
library(maps)
library(maptools)

# Load your rasters using terra
r1 <- rast("~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_present_USA.tif")
r2_2050 <- rast("~/R/Intermediate_Wheatgrass_collab/SDM_Future_2050/IWG_SDM_Future_2050_USA.tif")

# Create a mask for values greater than 0.2
s <- app(r1, fun = function(x) { ifelse(x > 0.2, 1, NA) })
s2 <- app(r2_2050, fun = function(x) { ifelse(x > 0.2, 3, NA) })

# Load USA shapefile
usa <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))  # Convert US states map to sf object

# Reproject USA boundaries to match the CRS of the raster
usa <- st_transform(usa, crs = st_crs(r1))

# Crop rasters to USA boundaries
s_usa <- crop(s, usa)
s_usa <- mask(s_usa, usa)

s2_usa <- crop(s2, usa)
s2_usa <- mask(s2_usa, usa)

# Plotting
pdf("usa_overlay_with_state_borders_02.pdf", width = 11, height = 8.5)  # Save to PDF
plot(st_geometry(usa), col = "#f2f2f2", border = "#a0a0a0", lwd = 0.5, 
     main = "USA Overlay of Raster Datasets (Cutoff > 0.2)")
colors <- c("#00FF7F44", "#B2222244")  # RGBA values for Present and 2050
plot(s_usa, col = colors[1], legend = FALSE, add = TRUE)
plot(s2_usa, col = colors[2], legend = FALSE, add = TRUE)

# Add a custom legend
legend("topright", legend = c("Present", "2050"), fill = colors, bty = "n")

dev.off()  # Close the PDF device

#####################################################################################################################
# Pixel Count (LOCAL)
######################################
#USA
present_tiff <- "~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_present_USA.tif"
future_tiff <- "~/R/Intermediate_Wheatgrass_collab/SDM_Future_2050/IWG_SDM_Future_2050_USA.tif"

#WORLD
#present_tiff <- "/Users/annamccormick/R/barley_collab/Collab_Proposal_SDMs/Barley_SDM_presentday.tif"
#future_tiff <- "/Users/annamccormick/R/barley_collab/Collab_Proposal_SDMs/13_CMIP_tif/13_model_SDM_future_2070.tif"

# Function to calculate area in km² for a geographic CRS
calculate_area_km2 <- function(tiff_path, threshold, latitude = 0) {
  # Read the raster file
  raster_data <- raster(tiff_path)
  
  # Get the resolution in degrees
  res_x <- res(raster_data)[1]
  res_y <- res(raster_data)[2]
  
  # Convert resolution to km using latitude adjustment
  km_per_degree <- 111.32  # Approximate conversion for 1 degree to km at the equator
  pixel_width_km <- res_x * km_per_degree
  pixel_height_km <- res_y * km_per_degree * cos(latitude * pi / 180)  # Adjust for latitude
  
  # Calculate pixel area in km²
  pixel_area_km2 <- pixel_width_km * pixel_height_km
  
  # Count the number of pixels above the threshold
  pixel_count <- sum(raster_data[] > threshold, na.rm = TRUE)
  
  # Calculate the total area in km²
  total_area_km2 <- pixel_count * pixel_area_km2
  
  # Return results
  list(pixel_count = pixel_count, total_area_km2 = total_area_km2)
}

# Set the latitude of your study area (e.g., 40°N for the USA)
latitude <- 40

# Calculate area for present and future TIFFs
present_results <- calculate_area_km2(present_tiff, 0.2, latitude)
future_results <- calculate_area_km2(future_tiff, 0.2, latitude)

# Print the results
cat("Present Day:\n")
cat("Pixels above 0.2:", present_results$pixel_count, "\n")
cat("Total area above 0.2 (km²):", present_results$total_area_km2, "\n\n")

cat("Future:\n")
cat("Pixels above 0.2:", future_results$pixel_count, "\n")
cat("Total area above 0.2 (km²):", future_results$total_area_km2, "\n")

##########################################
# SDM Blankets 0.1 for USA
##########################################

# Load necessary libraries
library(raster)
library(terra)
library(sf)
library(maps)
library(maptools)

# Load your rasters using terra
r1 <- rast("~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_present_USA.tif")
r2_2050 <- rast("~/R/Intermediate_Wheatgrass_collab/SDM_Future_2050/IWG_SDM_Future_2050_USA.tif")

# Create a mask for values greater than 0.1
s <- app(r1, fun = function(x) { ifelse(x > 0.1, 1, NA) })
s2 <- app(r2_2050, fun = function(x) { ifelse(x > 0.1, 3, NA) })

# Load USA shapefile
usa <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))  # Convert US states map to sf object

# Reproject USA boundaries to match the CRS of the raster
usa <- st_transform(usa, crs = st_crs(r1))

# Crop rasters to USA boundaries
s_usa <- crop(s, usa)
s_usa <- mask(s_usa, usa)

s2_usa <- crop(s2, usa)
s2_usa <- mask(s2_usa, usa)

# Plotting
pdf("usa_overlay_with_state_borders_01.pdf", width = 11, height = 8.5)  # Save to PDF
plot(st_geometry(usa), col = "#f2f2f2", border = "#a0a0a0", lwd = 0.5, 
     main = "USA Overlay of Raster Datasets (Cutoff > 0.1)")
colors <- c("#00FF7F44", "#B2222244")  # RGBA values for Present and 2050
plot(s_usa, col = colors[1], legend = FALSE, add = TRUE)
plot(s2_usa, col = colors[2], legend = FALSE, add = TRUE)

# Add a custom legend
legend("topright", legend = c("Present", "2050"), fill = colors, bty = "n")

dev.off()  # Close the PDF device

#####################################################################################################################
# Pixel Count (LOCAL)
######################################
#USA
present_tiff <- "~/R/Intermediate_Wheatgrass_collab/SDM_Present_day/IWG_SDM_present_USA.tif"
future_tiff <- "~/R/Intermediate_Wheatgrass_collab/SDM_Future_2050/IWG_SDM_Future_2050_USA.tif"

#WORLD
#present_tiff <- "/Users/annamccormick/R/barley_collab/Collab_Proposal_SDMs/Barley_SDM_presentday.tif"
#future_tiff <- "/Users/annamccormick/R/barley_collab/Collab_Proposal_SDMs/13_CMIP_tif/13_model_SDM_future_2070.tif"

# Function to calculate area in km² for a geographic CRS
calculate_area_km2 <- function(tiff_path, threshold, latitude = 0) {
  # Read the raster file
  raster_data <- raster(tiff_path)
  
  # Get the resolution in degrees
  res_x <- res(raster_data)[1]
  res_y <- res(raster_data)[2]
  
  # Convert resolution to km using latitude adjustment
  km_per_degree <- 111.32  # Approximate conversion for 1 degree to km at the equator
  pixel_width_km <- res_x * km_per_degree
  pixel_height_km <- res_y * km_per_degree * cos(latitude * pi / 180)  # Adjust for latitude
  
  # Calculate pixel area in km²
  pixel_area_km2 <- pixel_width_km * pixel_height_km
  
  # Count the number of pixels above the threshold
  pixel_count <- sum(raster_data[] > threshold, na.rm = TRUE)
  
  # Calculate the total area in km²
  total_area_km2 <- pixel_count * pixel_area_km2
  
  # Return results
  list(pixel_count = pixel_count, total_area_km2 = total_area_km2)
}

# Set the latitude of your study area (e.g., 40°N for the USA)
latitude <- 40

# Calculate area for present and future TIFFs
present_results <- calculate_area_km2(present_tiff, 0.1, latitude)
future_results <- calculate_area_km2(future_tiff, 0.1, latitude)

# Print the results
cat("Present Day:\n")
cat("Pixels above 0.1:", present_results$pixel_count, "\n")
cat("Total area above 0.1 (km²):", present_results$total_area_km2, "\n\n")

cat("Future:\n")
cat("Pixels above 0.1:", future_results$pixel_count, "\n")
cat("Total area above 0.1 (km²):", future_results$total_area_km2, "\n")


