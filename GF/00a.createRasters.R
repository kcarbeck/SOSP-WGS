# create "present" and "future" rasters from downloaded rasters from climateNA
# katherine carbeck
# 12 dec 2023

library(raster)
library(rgdal)
library(reshape2)
library(tidyverse)


setwd("C:/Users/kcarbeck/OneDrive - UBC/Desktop/ClimNAdat/na4000")


folder_path <- getwd()
# list environmental variables
variables <- c('DD_0', 'DD18', 'EMT', 'MAT', 'MSP', 'SHM', 'TD', 'RH')

# read in by looping through each year
for (year in 2001:2022) {
  # create a list to store raster objects for each variable
  rasters <- list()
  
  # then loop through each environmental variable
  for (variable in variables) {
    # construct file name
    file_name <- paste0("Year_", year, "Y", "/", variable, ".asc")
    
    # read raster data
    raster_obj <- raster(file.path(folder_path, file_name))
    
    # assign name to the raster object
    assign(paste0(variable, "_", year), raster_obj, envir = .GlobalEnv)
    
    # store raster object in the list
    rasters[[variable]] <- raster_obj
  }
}



###*##################################*###
    ###*  average for 2001-2011  *###
###*###################################*###

# create empty list to store the average raster objects
average_rasters <- list()

# loop through each environmental variable
for (variable in variables) {
  # create an empty raster stack to store raster layers for the specified years
  raster_stack <- stack()
  
  # loop through each year from 2001 to 2011
  for (year in 2001:2011) {
    # construct raster object name
    raster_name <- paste0(variable, "_", year)
    
    # load raster object from the global environment
    raster_obj <- get(raster_name, envir = .GlobalEnv)
    
    # add raster layer to the raster stack
    raster_stack <- stack(raster_stack, raster_obj)
  }
  
  # calculate the mean of the raster stack
  average_raster <- mean(raster_stack, na.rm = TRUE)
  
  # assign a name to the average raster object
  assign(paste0(variable, "_average"), average_raster, envir = .GlobalEnv)
  
  # store the average raster object in the list
  average_rasters[[variable]] <- average_raster

  # construct file name for the average raster file
  output_file <- file.path("2001_2011", paste0(variable, ".tif"))
  
  # save the average raster as GeoTIFF
  writeRaster(average_raster, filename = output_file, format = "GTiff", overwrite = TRUE)
}

crs(MAT_2001)
crs(MAT_average)
plot(MAT_average)



###*##################################*###
    ###*  average for 2012-2022  *###
###*###################################*###

# create an empty list to store the average raster objects
average_rasters <- list()

# loop through each environmental variable
for (variable in variables) {
  # create an empty raster stack to store raster layers for the specified years
  raster_stack <- stack()
  
  # loop through each year from 2001 to 2011
  for (year in 2012:2022) {
    # construct raster object name
    raster_name <- paste0(variable, "_", year)
    
    # load raster object from the global environment
    raster_obj <- get(raster_name, envir = .GlobalEnv)
    
    # add raster layer to the raster stack
    raster_stack <- stack(raster_stack, raster_obj)
  }
  
  # calculate the mean of the raster stack
  average_raster <- mean(raster_stack, na.rm = TRUE)
  
  # assign a name to the average raster object
  assign(paste0(variable, "_average"), average_raster, envir = .GlobalEnv)
  
  # store the average raster object in the list
  average_rasters[[variable]] <- average_raster
  
  # construct the file name for the average raster file
  output_file <- file.path("2012_2022", paste0(variable, ".tif"))
  
  # save the average raster as GeoTIFF
  writeRaster(average_raster, filename = output_file, format = "GTiff", overwrite = TRUE)
}

crs(MAT_2012)
crs(MAT_average)
plot(MAT_average)


