
message(noquote("Loading libraries ..."))
require(rgdal)
require(raster)
require(caret)
require(dplyr)
require(sf)
require(sp)
require(ranger)
require(randomForest)
require(hydroGOF)
require(Metrics)
require(ggplot2)
require(quantregForest)

rm(list = ls())

setwd(path_to_input)
message(noquote("reading shape extent..."))
aoi <- readOGR("wheat_mask.shp")

message(noquote("reading covariate and masking..."))
stack_normal <- stack("stack_eth_normal.tif") %>% crop(aoi) %>% crop(aoi)

names(stack_normal) 
  c(
    "bdod",
    "cec",
    "cfvo",
    "clay",
    "dem",
    "landform",
    "nitrogen",
    "ocd",
    "ocs",
    "phh2o",
    "sand",
    "silt",
    "slope",
    "soc",
    "tpi",
    "tri", 
    "prcp1",
    "prcp2",
    "prcp3",       
    "srad1",
    "srad2",        
    "srad3",
    "tmax1",        
    "tmax2",       
    "tmax3",       
    "tmin1",      
    "tmin2",       
    "tmin3"
  )

# ------------------------------------------------------------------------------
# read csv data and convert to spatial
  points <- read.csv("eth_csv.csv", header = T, sep = ",") %>% 
    select(points, -c("lon","lat"))
  new_crs <- proj4string(stack_normal)
  coordinates(points) <- ~lon+lat
  proj4string(points) <- new_crs