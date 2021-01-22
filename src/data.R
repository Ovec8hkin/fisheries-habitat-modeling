library(dplyr)
library(rgdal)
source("https://gist.githubusercontent.com/TonyLadson/f37aab3e2ef517188a7f27166307c985/raw/0822970769bc90fcc28052a91b375399d665286e/UTM2deg.R")

get_ecomon_catch_data <- function(species, fname){
  catch_data = read.csv(fname)
  
  catch_data = dplyr::select(catch_data, -c("X", "year")) #drop unnecsarry columns
  catch_data = dplyr::rename(catch_data, c("catch"=species)) # restrict to species of interest
  catch_data$btm_sub = factor(catch_data$btm_sub)
  return(catch_data)
}

convert_to_northing_easting <- function(lat, lon){
  coords.dec = SpatialPoints(cbind(lon, lat), proj4string=CRS("+proj=longlat"))
  coords.UTM = spTransform(coords.dec, CRS("+init=epsg:26983"))
  cs = coordinates(coords.UTM)
  colnames(cs)=c("easting", "northing")
  return(cs)
}

load_environmental_data <- function(y, m, s, rivers=FALSE){
  
  month = str_pad(m, 2, pad="0")
  fname=paste("auxdata/environmental-data/environmental_data_", y, month, ".csv", sep="")
  prediction.data = read.csv(fname)
  colnames(prediction.data) = c("east", "north", "sfc_temp", "btm_temp", "sfc_salt", "btm_salt", "depth", "area")
  prediction.data$area_km2 = prediction.data$area/1000/1000
  edata = prediction.data
  
  if(!rivers){
    edata = remove_rivers(prediction.data)
  }
  
  len = length(edata[,1])
  edata$year = rep(y, len)
  edata$month = rep(m, len)
  edata$species = rep(s, len)
  
  return(edata)
  
}

remove_rivers <- function(edata){
  sub = (edata$sfc_salt > 30 & edata$sfc_salt < 35)
  edata = edata[sub,]
  return(edata)
}