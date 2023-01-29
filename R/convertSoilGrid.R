# ------------------------------------------------------------------------------
# converting list of SoilGrid geo-tiff files using conversion factors
# and calculating weighted average
# ------------------------------------------------------------------------------
convertSoilGrids <- function(sou_path, dest_path = NULL){
  require(raster)
  require(rgdal)
  library(dplyr)
  
  setwd(sou_path) #soil grid work space
  if(is.null(dest_path)){
    dest_path = getwd()
  }
  else{
    dest_path = dest_path
  }
  rasterlist <- list.files(path = ".", pattern = ".tif$")
  rasters <- lapply(rasterlist, raster)  %>% stack()
  soil_stack <- stack()
  
  #multiple rasters by factor
  for(i in 1:nlayers(rasters) ) { 
    if("bdod" %in% names(rasters[[i]]) ||
       "nitrogen" %in% names(rasters[[i]]))
    {
      rout <- rasters[[i]] * 0.01
    }else{
      rout <- rasters[[i]] * 0.1
    }
    soil_stack <- addLayer(soil_stack, rout)
  }
  
  #create a weighted average for each parameter
  stack_weigh <- stack()
  for(i in seq(from = 1, to = nlayers(soil_stack), by = 3)) { 
    ras <- soil_stack[[i]] * (5 / 30) + soil_stack[[i+1]] * (10 / 30) + soil_stack[[i+2]] * (15 / 30)
    stack_weigh <- addLayer(stack_weigh, ras)
  }
  
  #write the raster stack
  writeRaster(
    stack_weigh,
    filename = paste(dest_path, "soil_stack.tif", sep = "/"), 
    format = "GTiff",
    overwrite = T)
}

