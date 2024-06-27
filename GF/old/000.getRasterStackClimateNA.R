# download rasterstack for climate vars from climateNA
# 20 oct 2023
# katherine carbeck

setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS/GF")

#install climateNA
devtools::install_local("/Users/katherine/Downloads/ClimateNAr.zip", repos=NULL, type="source") # install climateNAr
library(ClimateNAr)

varList <- c('mat', 'td','msp', 'dd_0', 'dd18', 'emt', 'rh')  


clm <- ClimateNA_API2(ClimateBC_NA='NA', inputFile=env, period='Normal_1961_1990.nrm',MSY='SY') 
head(clm);dim(clm) 
write.csv(clm, file="data/climateNA_output.csv", row.names = FALSE)  



rasterDownload (region='BC',res='800m',period='Normal_1961_1990',varList=varList,sDir='C:/temp')  






#define projection 
wkt <- sf::st_crs(4326)[[2]] #projection for lat-long coordinates - WGS84 ellipsoid
#data source 'src' is Copernicus Digital Elevation Model (DEM) - an AMAZON Sustainability Data Initiative
#'z' aka zoom, determines ground resolution, at 45 degrees lat z = 9 is ~216 meters resolution (https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-the-ground-resolution)
occurrence.elevation <- get_elev_point(locations = occurrence.location, units = 'meters', prj = wkt, src = 'aws', z = 9) #~3 minute processing 




# elevation file from adaptwest - elevation.tif

Coordinate Systems
------------------
  
  X,Y: 
  Projection: Lambert_Azimuthal_Equal_Area
False_Easting: 0.0
False_Northing: 0.0
Central_Meridian: -100.0
Latitude_Of_Origin: 45.0
Linear Unit: Meter 
Geographic Coordinate System: GCS_WGS_1984
Datum: D_WGS_1984

ELEV:
  Elevation in meters above sea level Source: GMTED2010 data at 15 arc-second resolution. Danielson, J.J., and Gesch, D.B., 2011, Global multi-resolution terrain elevation data 2010 (GMTED2010): USGS Open-File Report 2011–1073. https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-multi-resolution-terrain-elevation





































setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS")
library(dplyr)
library(sf)
library(tidyr)
library(stringr)
library(rgbif) # Interface to the Global 'Biodiversity' Information Facility API




###### ---------------- Create input file -------------- #######
# read in initial df
# ClimateNAr requires a properly formatted .csv input file that includes the following headers: ID1, ID2, lat, long, el
env<-read.csv("data/climateNA_input.csv", stringsAsFactors = F)
env<- env %>%
  rename(c(ID1=sample_id, 
           ID2=subspecies))
head(env)
# ID1     ID2      Lat      Long
# 1 AK_caurina_S11 caurina 60.41750 -145.4039
# 2 AK_caurina_S12 caurina 60.41750 -145.4039
# 3 AK_caurina_S13 caurina 60.35100 -145.4200
# 4 AK_caurina_S14 caurina 60.35100 -145.4200
# 5 AK_caurina_S15 caurina 62.60757 -144.5637
# 6 AK_caurina_S16 caurina 60.41750 -145.4039
str(env)
env$lat <- as.numeric(env$lat)
env$long <- as.numeric(env$long)


# extract elevation using rgbif
# To get a GeoNames user name, register for an account at http://www.geonames.org/login- then you can enable your account for the GeoNames webservice on your account page (http://www.geonames.org/manageaccount). Once you are enabled to use the webservice, you can pass in your username to the username parameter. Better yet, store your username in your .Renviron file, or similar (e.g., .zshrc or .bash_profile files) and read it in via Sys.getenv() as in the examples below. By default we do Sys.getenv("GEONAMES_USER") for the username parameter.
user<-Sys.getenv("kcarbeck") 
el <- elevation(latitude = env$lat, longitude = env$long, elevation_model = "srtm1", 
                username="kcarbeck")
env<- env %>%
  bind_cols(el=el$elevation_geonames) 
dim(env)
write.csv(env, file="data/climateNA_input_el.csv", row.names = FALSE)  




### Extract ClimateNA monthly summaries directly with R 
devtools::install_local("/Users/katherine/Downloads/ClimateNAr.zip", repos=NULL, type="source") # install climateNAr
library(ClimateNAr)

clm <- ClimateNA_API2(ClimateBC_NA='NA', inputFile=env, period='Normal_1961_1990.nrm',MSY='SY') 
head(clm);dim(clm) 
write.csv(clm, file="data/climateNA_output.csv", row.names = FALSE)  

