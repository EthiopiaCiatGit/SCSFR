# ------------------------------------------------------------------------------   
# Function to access and download ISRIC Soil grids 
# ------------------------------------------------------------------------------

fetchSoilGrid <-
  function(voi, depth, quantile) {
    library(XML)
    library(sf)
    library(rgdal)
    library(gdalUtils)
    
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
    
# ------------------------------------------------------------------------------   
# bbox and projection
    bbox <- c(32.999939, 14.899958, 47.986179, 3.322099)
    wgs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  

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
      projwin = bbox,
      projwin_srs = wgs,
      co = c(
        "TILED=YES",
        "COMPRESS=DEFLATE",
        "PREDICTOR=2",
        "BIGTIFF=YES"
      ),
      verbose = TRUE
    )
  }


