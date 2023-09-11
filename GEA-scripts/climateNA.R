# create environmental variables data file for GEA/GF 
# 22 july 2023
# katherine carbeck

setwd("/Users/katherine/Desktop/PhD/github/song-sparrow-WGS")
library(dplyr)
library(sf)
library(tidyr)
library(stringr)
library(rgbif) # Interface to the Global 'Biodiversity' Information Facility API




###### ---------------- Create input file -------------- #######
# read in initial df
# ClimateNAr requires a properly formatted .csv input file that includes the following headers: ID1, ID2, lat, long, el
env<-read.csv("GEA-data/ClimateNA_input.csv", stringsAsFactors = F)
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



install.packages("elevatr")
library(elevatr)
library(sf)
library(sp)

#subset for just coordinates
head(env)
sub<-subset(env, select=c(long,lat))
head(sub)

# elevatr package to get point elevation
?get_elev_point
prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
elev <- elevatr::get_elev_point(sub, prj = prj_dd, src="aws",z=7)
elev$elevation

#bind env data with elevation
env_el <- cbind(env,elev$elevation)


write.csv(env_el, file="GEA-data/ClimateNA_input_el2.csv", row.names = FALSE)  

env_el<- read.csv("GEA-data/ClimateNA_input_el2.csv")
head(env_el)


### Extract ClimateNA monthly summaries directly with R 
# can only do 100 lines at a time so split df into 100 lines each
t<-split(env_el, (seq(nrow(env_el))-1) %/% 100) 
lapply(t,dim)
env_el1<- t[["0"]]
env_el2<- t[["1"]]
env_el3<- t[["2"]]
env_el4<- t[["3"]]

devtools::install_local("/Users/katherine/Downloads/ClimateNAr_1.2.0.zip", repos=NULL, type="source") # install climateNAr
library(ClimateNAr)
?ClimateNA_API2

clm1 <- ClimateNA_API2(ClimateBC_NA='NA', inputFile=env_el1, period='Normal_1991_2020.nrm',MSY='SY') 
clm2 <- ClimateNA_API2(ClimateBC_NA='NA', inputFile=env_el2, period='Normal_1991_2020.nrm',MSY='SY') 
clm3 <- ClimateNA_API2(ClimateBC_NA='NA', inputFile=env_el3, period='Normal_1991_2020.nrm',MSY='SY') 
clm4 <- ClimateNA_API2(ClimateBC_NA='NA', inputFile=env_el4, period='Normal_1991_2020.nrm',MSY='SY') 

clm<-do.call(rbind, list(clm1, clm2, clm3, clm4)) #rbind all climatena output
write.csv(clm, file="GEA-data/ClimateNA_output_1991-2020.csv", row.names = FALSE)  








####################      old       ######################

# old method of extracting elevation:
# extract elevation using rgbif
# To get a GeoNames user name, register for an account at http://www.geonames.org/login- then you can enable your account for the GeoNames webservice on your account page (http://www.geonames.org/manageaccount). Once you are enabled to use the webservice, you can pass in your username to the username parameter. Better yet, store your username in your .Renviron file, or similar (e.g., .zshrc or .bash_profile files) and read it in via Sys.getenv() as in the examples below. By default we do Sys.getenv("GEONAMES_USER") for the username parameter.
user<-Sys.getenv("kcarbeck") 
el <- elevation(latitude = env$lat, longitude = env$long, elevation_model = "srtm1", 
                username="kcarbeck")
#bind elevation to env df
env<- env %>%
  bind_cols(el=el$elevation_geonames) 
dim(env)
write.csv(env, file="GEA-data/ClimateNA_input_el.csv", row.names = FALSE)  


