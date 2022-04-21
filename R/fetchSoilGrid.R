# ------------------------------------------------------------------------------   
# Function to access and download ISRIC Soil grids 
# ------------------------------------------------------------------------------

fetchSoilGrid <-
  function(voi, depth, quantile, aoi){
    require(XML)
    require(sf)
    require(sp)
    require(raster)
    require(rgdal)
    require(gdalUtils)
    
    aoi <- st_read(aoi)
    igh <-'+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs'
    aoi_igh <- st_transform(aoi, igh)
    bbox <- st_bbox(aoi_igh)
    ulx <- bbox$xmin
    uly <- bbox$ymax
    lrx <- bbox$xmax
    lry <- bbox$ymin
    bb <- c(ulx, uly, lrx, lry)
    
    url <- "/vsicurl/https://files.isric.org/soilgrids/latest/data/"
    
    cov <- voi 
    soil_dep <- depth
    quant <- quantile
    lfile <- paste0(as.character(cov),".tif")
    
    voi_path <- paste(cov, paste0(paste(cov,soil_dep, quant, sep = "_"),".vrt"), sep = "/")
    voi_url <- paste0(url,voi_path)
    
    gdal_translate(
      voi_url,
      lfile ,
      tr = c(250, 250),
      of = "GTiff",
      projwin = bb,
      projwin_srs = igh,
      verbose = TRUE
    )
  }
