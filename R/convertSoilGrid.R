# ------------------------------------------------------------------------------
# converting list of SoilGrid geo-tiff files using conversion factors
# ------------------------------------------------------------------------------
convertSoilGrids <- function(sou_path, dest_path = NULL){
  require(raster)
  require(rgdal)
  
  setwd(sou_path)
  
  if(is.null(dest_path)){
    dest_path = getwd()
  }
  else{
    dest_path = dest_path
  }
  rasterlist <- list.files(path = ".", pattern = ".tif$")
  rasters <- lapply(rasterlist, raster)
  for(i in 1:length(rasters) ) { 
    if("bdod" %in% names(rasters[[i]]) ||
       "nitrogen" %in% names(rasters[[i]]))
    {
      rout <- rasters[[i]] * 0.01
    }else{
      rout <- rasters[[i]] * 0.1
    }
    rout_name <- names(rasters[[i]])
    writeRaster(
      rout,
      filename = paste(dest_path, rout_name, sep = "/"), 
      format = "GTiff",
      overwrite = T
    )
  }
}

