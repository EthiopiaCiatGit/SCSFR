# ------------------------------------------------------------------------------
# Fetches elevation data from Amazon Web Services using and area of interest
# Projects to (WGS 1984 UTM Zone 37N) calculates slope, aspect, tpi, tri and 
# degree and exports the layers as a tiff
# landform. The projection of the aoi should be in gcs_wgs_84
# ------------------------------------------------------------------------------ 
fetchElevation <- function(aoi, aoi_path = NULL){
  library(raster)
  library(sf)
  library(rgl)
  library(rgdal)
  library(elevatr)
  
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
  # utm <- "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs"
  # dem <- projectRaster(dem, crs = utm, res = 250)
  
# ------------------------------------------------------------------------------
# derive slope and aspect
  slope <- terrain(dem,opt='slope', unit='degrees') #slop
  aspect <- terrain(dem,opt='aspect',unit='degrees') #aspect
  tpi <- terrain(dem,opt='TPI') #topographic position index
  tri <- terrain(dem,opt='TRI') #topographic ruggedness index 
  # aoi <- st_transform(aoi, utm)
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

#Example
fetchElevation(aoi = "tigray_prj", "D:/tigray")
fetchElevation(aoi = "bg_prj")
