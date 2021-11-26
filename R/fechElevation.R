# ------------------------------------------------------------------------------
# Fetches elevation data from Amazon Web Services using and area of interest
# ------------------------------------------------------------------------------ 
fetchElevation <- function(aoi, aoi_path = NULL){
  require(raster)
  require(sf)
  require(rgl)
  require(rgdal)
  require(elevatr)
  
# ------------------------------------------------------------------------------
# read area of interest 
  if (is.null (aoi_path)){
    path <- getwd()
    aoi_path <- setwd(path)
    aoi <- st_read(dsn = aoi_path, layer = aoi)
  }
  else {
    aoi_path <- aoi_path
    aoi <- st_read(dsn = aoi_path, layer = aoi)
  }
  
# ------------------------------------------------------------------------------
# fetches data using elevatr package 
  dem <- get_elev_raster(aoi, z = 8)
  
# ------------------------------------------------------------------------------
# derive slope and aspect
  slope <- terrain(dem,opt='slope', unit='degrees') 
  aspect <- terrain(dem,opt='aspect',unit='degrees') 
  tpi <- terrain(dem,opt='TPI') 
  tri <- terrain(dem,opt='TRI')  
  SD <- sd(tpi[],na.rm=T)
  landform <- reclassify(tpi, matrix(c(-Inf, -SD, 1,
                          -SD, -SD/2, 2,
                          -SD/2, 0, 3,
                          0, SD/2, 4,
                          SD/2, SD, 5,
                          SD, Inf, 6),
                          ncol = 3, byrow = T),
                          right = T)
  
  landform <- as.factor(landform) 
  rat <- levels(landform)[[1]]
  rat[["landform"]] <- c('Valley', 'Lower Slope', 
                         'Flat Area','Middle Slope', 
                         'Upper Slope', 'Ridge')
  levels(landform) <- rat 
  
# ------------------------------------------------------------------------------
# crop using aoi
  dem <- crop(dem, aoi)
  slope <- crop(slope, aoi)
  aspect <- crop(aspect, aoi)
  tpi <- crop(tpi, aoi)
  tri <- crop(tri, aoi)
  landform <- crop(landform, aoi)
# ------------------------------------------------------------------------------
# export raster 
  writeRaster(
    dem,
    filename = "dem",
    bylayer = TRUE,
    format = "GTiff",
    overwrite = T
  )
  writeRaster(
    slope,
    filename = "slope",
    bylayer = TRUE,
    format = "GTiff",
    overwrite = T
  )
  writeRaster(
    aspect,
    filename = "aspect",
    bylayer = TRUE,
    format = "GTiff",
    overwrite = T
  )
  writeRaster(
    tpi,
    filename = "tpi",
    bylayer = TRUE,
    format = "GTiff",
    overwrite = T
  ) 
  writeRaster(
    tri,
    filename = "tri",
    bylayer = TRUE,
    format = "GTiff",
    overwrite = T
  )
  writeRaster(
    landform,
    filename = "landform",
    bylayer = TRUE,
    format = "GTiff",
    overwrite = T
  )
}

