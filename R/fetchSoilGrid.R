# ------------------------------------------------------------------------------   
# Function to access and download ISRIC Soil grids 
# ------------------------------------------------------------------------------

fetchSoilGrid <-
  function(voi, depth, quantile, aoi_path, aoi){
    require(XML)
    require(sf)
    require(rgdal)
    require(gdalUtils)
    
    if(.Platform$OS.type == "windows"){
      gdal.dir <- shortPathName("C:/Program Files/GDAL")
      gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
      gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
      gdalinfo <- paste0(gdal.dir, "/gdalinfo.exe")
    } else {
      gdal_translate = "gdal_translate"
      gdalwarp = "gdalwarp"
      gdalinfo = "gdalinfo"
    }
    
    igh <- '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs'
    aoi <- st_read(paste(aoi_path, aoi, sep = "/"))
    aoi_igh <- st_transform(aoi, igh)
    bbox <- st_bbox(aoi_igh)
    ulx <- bbox$xmin
    uly <- bbox$ymax
    lrx <- bbox$xmax
    lry <- bbox$ymin
    bb <- c(ulx, uly, lrx, lry)
    # ------------------------------------------------------------------------------   
    # layer of interest
    voi_layer <-
      paste(voi, depth, quantile, sep = "_") 
    
    # ------------------------------------------------------------------------------  
    # Path to the WCS. See maps.isric.org
    wcs_path <-
      paste0("https://maps.isric.org/mapserv?map=/map/", voi, ".map") 
    wcs_service <- "SERVICE=WCS"
    
    # ------------------------------------------------------------------------------
    # This works for gdal >=2.3; "VERSION=1.1.1" works with gdal < 2.3.
    wcs_version <-
      "VERSION=2.0.1" 
    wcs_request <- "DescribeCoverage"
    wcs <-
      paste(wcs_path, wcs_service, wcs_version, wcs_request, sep = "&")
    l1 <- newXMLNode("WCS_GDAL")
    l1.s <- newXMLNode("ServiceURL", wcs, parent = l1)
    l1.l <- newXMLNode("CoverageName", voi_layer, parent = l1)
    xml.out <- paste(voi,".xml")
    saveXML(l1, file = xml.out)
    file.out <- paste(voi, "tif", sep = ".")
    gdal_translate(
      xml.out,
      file.out,
      tr = c(250, 250),
      of = "GTiff",
      projwin = bb,
      projwin_srs = igh,
      co = c(
        "TILED=YES",
        "COMPRESS=DEFLATE",
        "PREDICTOR=2",
        "BIGTIFF=YES",
        "GDAL_HTTP_UNSAFESSL=YES"
      ),
      verbose = TRUE
    )
  }


