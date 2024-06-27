## 7 dec 2023
## katherine carbeck
## use shapefile range map created from 'eBirdstGetRangeMap.R' to clip climate rasters created in '00a.createRasters.R'

## we need range-wide climate data to fit trained gradient forest models to current/future climate
## transform  ClimateNA data gridded climate data from Lambert conformal to WGS84
##then extract climate data within the boundaries of range map shapefile (eBirdstGetRangeMap.R)

# steps:
# 1. download climateNA environmental vars rasters ('00a.createRasters.R')
# 2. clip climateNA rasters to shapefile mask
# 3. standardize values
# 4. create raster stack and extract values by lat long

######! bash
# copy raster files from home to workdir
cp -R /home/lc736_0001/song_sparrow/final_vcf/GF/climNAdat/2001_2011 /workdir/kcarbeck/data
cp -R /home/lc736_0001/song_sparrow/final_vcf/GF/climNAdat/2012_2022 /workdir/kcarbeck/data

# copy shapefiles from home to workdir - merged_range_shapefile.shp
cp -R /home/lc736_0001/song_sparrow/final_vcf/GF/eBird /workdir/kcarbeck/data

#########!

# install.packages("raster")
# install.packages("rgdal")
library(raster)
library(rgdal)
library(reshape2)

setwd("/workdir/kcarbeck/data")

# Load the shapefile
shapefile <- readOGR(dsn = "/workdir/kcarbeck/data/eBird/", layer = "merged_range_shapefile")
    # OGR data source with driver: ESRI Shapefile
    # Source: "/local/workdir/kcarbeck/data/eBird", layer: "merged_range_shapefile"
    # with 1 features
    # It has 1 fields
    # Integer64 fields read as strings:  FID
print(proj4string(shapefile))
    # "+proj=longlat +datum=WGS84 +no_defs"

# plot
png("shp.png", width = 800, height = 600, units = "px", res = 100)
print(plot(shapefile, main = ""))
dev.off()
#!bash:
#! imgcat /workdir/kcarbeck/data/shp.png



###*##################################*###
 ###*     clip 2001-2011 rasters     *###
###*###################################*###
setwd("/workdir/kcarbeck/data/2001_2011")
# List of raster files
raster_files <- c("DD_0.tif", "DD18.tif", "EMT.tif", "MAT.tif", "MSP.tif", "RH.tif", "SHM.tif", "TD.tif")

# Create an empty list to store the masked rasters
masked_rasters <- list()

# Loop through each raster file, clip, and mask to the shapefile
for (raster_file in raster_files) {
  raster_data <- raster(raster_file)
  
  # Ensure that the raster has the same CRS as the shapefile
  projection(raster_data) <- projection(shapefile)
  
  # Clip the raster to the extent of the shapefile
  clipped_raster <- crop(raster_data, extent(shapefile))
  
  # Mask the raster to the shapefile polygons
  masked_raster <- mask(clipped_raster, shapefile)
  
  # Add the masked raster to the list
  masked_rasters[[raster_file]] <- masked_raster
}

# save masked rasters
for (i in seq_along(raster_files)) {
  writeRaster(masked_rasters[[i]], filename = paste0("clipped_", raster_files[i]), format = "ascii")
}



#### plot #####
dd18<-masked_rasters[[2]]
NAvalue(dd18)
#-Inf

png("dd18.png", width = 800, height = 600, units = "px", res = 100)
print(plot(dd18, main = "dd18"))
dev.off()
#! imgcat /workdir/kcarbeck/data/2001_2011/dd18.png



# dd18 <- masked_rasters[["DD18.asc"]]
# png("dd18_w_shp.png", width = 800, height = 600, units = "px", res = 100)
# plot(dd18, main = "Clipped Raster DD18")
# plot(shapefile, add = TRUE, border = "red")  # Overlay shapefile in red
# dev.off()
#! imgcat /workdir/kcarbeck/data/Normal_1991_2020Y/dd18_w_shp.png





###*####################################*###
###*   standardize values (2001-2011)  *###
###*####################################*###

# Create an empty list to store the standardized rasters
standardized_rasters <- list()

# Loop through each clipped raster and apply scaling
for (i in seq_along(masked_rasters)) {
  # Extract values from the raster
  raster_values <- getValues(masked_rasters[[i]])
  
  # Scale the values
  scaled_values <- scale(raster_values)
  
  # Create a new raster with scaled values
  standardized_raster <- setValues(masked_rasters[[i]], scaled_values)
  
  # Add the standardized raster to the list
  standardized_rasters[[i]] <- standardized_raster
}

# save standardized clipped rasters
for (i in seq_along(raster_files)) {
  writeRaster(standardized_rasters[[i]], filename = paste0("standardized_", raster_files[i]), format = "ascii")
}

dd18<-standardized_rasters[[2]]
NAvalue(dd18)

#plot
png("dd18.png", width = 800, height = 600, units = "px", res = 100)
print(plot(dd18, main = "Standardized Raster dd18"))
dev.off()
#! imgcat /workdir/kcarbeck/data/2001_2011/dd18.png

crs(dd18)

###*##################################*###

 ###*   extract values to txt file   *###

###*###################################*###


raster_stack <- stack(standardized_rasters)
vals <- values(raster_stack)
coord <- xyFromCell(raster_stack,1:ncell(raster_stack))
combine <- cbind(coord,vals)
colnames(combine)[colnames(combine) == "x"] <- "lon"
colnames(combine)[colnames(combine) == "y"] <- "lat"
#            lon      lat DD_0 DD18 EMT MAT MSP RH SHM TD
# [1,] -178.6667 64.07783   NA   NA  NA  NA  NA NA  NA NA
# [2,] -178.6250 64.07783   NA   NA  NA  NA  NA NA  NA NA

# Write the combined data to a text file
write.table(combine, file = "combined_raster_data_2001_2011.txt", sep = "\t", row.names = FALSE)




#####################*#####################

#####################*#####################

#####################*#####################

#####################*#####################


###*##################################*###
 ###*     clip 2012-2022 rasters     *###
###*###################################*###
setwd("/workdir/kcarbeck/data/2012_2022")
# List of raster files
raster_files <- c("DD_0.tif", "DD18.tif", "EMT.tif", "MAT.tif", "MSP.tif", "RH.tif", "SHM.tif", "TD.tif")

# Create an empty list to store the masked rasters
masked_rasters <- list()

# Loop through each raster file, clip, and mask to the shapefile
for (raster_file in raster_files) {
  raster_data <- raster(raster_file)
  
  # Ensure that the raster has the same CRS as the shapefile
  projection(raster_data) <- projection(shapefile)
  
  # Clip the raster to the extent of the shapefile
  clipped_raster <- crop(raster_data, extent(shapefile))
  
  # Mask the raster to the shapefile polygons
  masked_raster <- mask(clipped_raster, shapefile)
  
  # Add the masked raster to the list
  masked_rasters[[raster_file]] <- masked_raster
}

# save masked rasters
for (i in seq_along(raster_files)) {
  writeRaster(masked_rasters[[i]], filename = paste0("clipped_", raster_files[i]), format = "ascii")
}



#### plot #####
dd18<-masked_rasters[[2]]
NAvalue(dd18)
#-Inf
crs(dd18)

png("dd18.png", width = 800, height = 600, units = "px", res = 100)
print(plot(dd18, main = "dd18"))
dev.off()
#! imgcat /workdir/kcarbeck/data/2012_2022/dd18.png


###*####################################*###
###*   standardize values (2012-2022)  *###
###*####################################*###

# Create an empty list to store the standardized rasters
standardized_rasters <- list()

# Loop through each clipped raster and apply scaling
for (i in seq_along(masked_rasters)) {
  # Extract values from the raster
  raster_values <- getValues(masked_rasters[[i]])
  
  # Scale the values
  scaled_values <- scale(raster_values)
  
  # Create a new raster with scaled values
  standardized_raster <- setValues(masked_rasters[[i]], scaled_values)
  
  # Add the standardized raster to the list
  standardized_rasters[[i]] <- standardized_raster
}

# save standardized clipped rasters
for (i in seq_along(raster_files)) {
  writeRaster(standardized_rasters[[i]], filename = paste0("standardized_", raster_files[i]), format = "ascii")
}

dd18<-standardized_rasters[[2]]
NAvalue(dd18)
#-Inf

#plot
png("dd18.png", width = 800, height = 600, units = "px", res = 100)
print(plot(dd18, main = "Standardized Raster dd18"))
dev.off()
#! imgcat /workdir/kcarbeck/data/2012_2022/dd18.png

crs(dd18)
# Deprecated Proj.4 representation: +proj=longlat +datum=WGS84 +no_defs


###*##################################*###

 ###*   extract values to txt file   *###

###*###################################*###


raster_stack <- stack(standardized_rasters)
vals <- values(raster_stack)
coord <- xyFromCell(raster_stack,1:ncell(raster_stack))
combine <- cbind(coord,vals)
colnames(combine)[colnames(combine) == "x"] <- "lon"
colnames(combine)[colnames(combine) == "y"] <- "lat"
head(combine)
#            lon      lat DD_0 DD18 EMT MAT MSP RH SHM TD
# [1,] -178.6667 64.07783   NA   NA  NA  NA  NA NA  NA NA
# [2,] -178.6250 64.07783   NA   NA  NA  NA  NA NA  NA NA

# Write the combined data to a text file
write.table(combine, file = "combined_raster_data_2012_2022.txt", sep = "\t", row.names = FALSE)









###*##################################*###
       
       ###*  WEST COAST    *###

       ###*  DEC 19 2023    *###

###*###################################*###

library(raster)
library(rgdal)
library(reshape2)

setwd("/workdir/kcarbeck/data/2001_2011")

# load maskfile
mask_raster <- raster("/workdir/kcarbeck/data/maskdir/west_states_mask.tif")


###*##################################################*###
  ###*  WEST clip & standardize 2001-2011 rasters   *###
###*###################################################*###
# List of raster files
raster_files <- c("DD_0.tif", "DD18.tif", "EMT.tif", "MAT.tif", "MSP.tif", "RH.tif", "SHM.tif", "TD.tif")

# Create an empty list to store the masked rasters
masked_rasters <- list()

# Loop through each raster file, clip, and mask to the shapefile
for (raster_file in raster_files) {
  raster_data <- raster(raster_file)
  
  # Ensure that the raster has the same CRS and extent
  crs(raster_data) <- crs(mask_raster)
  clipped_raster <- crop(raster_data, extent(mask_raster))
  
  # Mask the raster to the shapefile polygons
  masked_raster <- mask(clipped_raster, mask_raster)
  
  # Add the masked raster to the list
  masked_rasters[[raster_file]] <- masked_raster
}

dd0<-masked_rasters[[1]]
png("dd0.png", width = 800, height = 600, units = "px", res = 100)
print(plot(dd0, main = "dd0"))
dev.off()
#! imgcat dd0.png



# Create an empty list to store the standardized rasters
standardized_rasters <- list()

# Loop through each clipped raster and apply scaling
for (i in seq_along(masked_rasters)) {
  # Extract values from the raster
  raster_values <- getValues(masked_rasters[[i]])
  
  # Scale the values
  scaled_values <- scale(raster_values)
  
  # Create a new raster with scaled values
  standardized_raster <- setValues(masked_rasters[[i]], scaled_values)
  
  # Add the standardized raster to the list
  standardized_rasters[[i]] <- standardized_raster
}

# save standardized clipped rasters
for (i in seq_along(raster_files)) {
  writeRaster(standardized_rasters[[i]], filename = paste0("west_standardized_", raster_files[i]), format = "ascii", overwrite=TRUE)
}

dd0<-standardized_rasters[[1]]
png("dd0.png", width = 800, height = 600, units = "px", res = 100)
print(plot(dd0, main = "Standardized Raster dd0"))
dev.off()


###*   extract values to txt file   *###
raster_stack <- stack(standardized_rasters)
vals <- values(raster_stack)
coord <- xyFromCell(raster_stack,1:ncell(raster_stack))
combine <- cbind(coord,vals)
colnames(combine)[colnames(combine) == "x"] <- "lon"
colnames(combine)[colnames(combine) == "y"] <- "lat"
#            lon      lat DD_0 DD18 EMT MAT MSP RH SHM TD
# [1,] -178.6667 64.07783   NA   NA  NA  NA  NA NA  NA NA
# [2,] -178.6250 64.07783   NA   NA  NA  NA  NA NA  NA NA

# Write the combined data to a text file
write.table(combine, file = "west_combined_raster_data_2001_2011.txt", sep = "\t", row.names = FALSE)




###*##################################################*###
  ###*  WEST clip & standardize 2012-2022 rasters   *###
###*###################################################*###
setwd("/workdir/kcarbeck/data/2012_2022")
# List of raster files
raster_files <- c("DD_0.tif", "DD18.tif", "EMT.tif", "MAT.tif", "MSP.tif", "RH.tif", "SHM.tif", "TD.tif")

# Create an empty list to store the masked rasters
masked_rasters <- list()

# Loop through each raster file, clip, and mask to the shapefile
for (raster_file in raster_files) {
  raster_data <- raster(raster_file)
  
  # Ensure that the raster has the same CRS and extent
  crs(raster_data) <- crs(mask_raster)
  clipped_raster <- crop(raster_data, extent(mask_raster))
  
  # Mask the raster to the shapefile polygons
  masked_raster <- mask(clipped_raster, mask_raster)
  
  # Add the masked raster to the list
  masked_rasters[[raster_file]] <- masked_raster
}

dd0<-masked_rasters[[1]]
png("dd0.png", width = 800, height = 600, units = "px", res = 100)
print(plot(dd0, main = "dd0"))
dev.off()
#! imgcat dd0.png



# Create an empty list to store the standardized rasters
standardized_rasters <- list()

# Loop through each clipped raster and apply scaling
for (i in seq_along(masked_rasters)) {
  # Extract values from the raster
  raster_values <- getValues(masked_rasters[[i]])
  
  # Scale the values
  scaled_values <- scale(raster_values)
  
  # Create a new raster with scaled values
  standardized_raster <- setValues(masked_rasters[[i]], scaled_values)
  
  # Add the standardized raster to the list
  standardized_rasters[[i]] <- standardized_raster
}

# save standardized clipped rasters
for (i in seq_along(raster_files)) {
  writeRaster(standardized_rasters[[i]], filename = paste0("west_standardized_", raster_files[i]), format = "ascii", overwrite=TRUE)
}

dd0<-standardized_rasters[[1]]
png("dd0.png", width = 800, height = 600, units = "px", res = 100)
print(plot(dd0, main = "Standardized Raster dd0"))
dev.off()


###*   extract values to txt file   *###
raster_stack <- stack(standardized_rasters)
vals <- values(raster_stack)
coord <- xyFromCell(raster_stack,1:ncell(raster_stack))
combine <- cbind(coord,vals)
colnames(combine)[colnames(combine) == "x"] <- "lon"
colnames(combine)[colnames(combine) == "y"] <- "lat"
#            lon      lat DD_0 DD18 EMT MAT MSP RH SHM TD
# [1,] -178.6667 64.07783   NA   NA  NA  NA  NA NA  NA NA
# [2,] -178.6250 64.07783   NA   NA  NA  NA  NA NA  NA NA

# Write the combined data to a text file
write.table(combine, file = "west_combined_raster_data_2012_2022.txt", sep = "\t", row.names = FALSE)
