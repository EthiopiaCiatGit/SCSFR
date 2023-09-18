
# --------------------------------------------------------------------
# A function for downloading chirps precipitation data from GEE
# --------------------------------------------------------------------

calcfactorR <- function(inp_wd, aoi, st_date, en_date){
  
  require(rgee)
  require(reticulate)
  require(sf)
  require(terra)
  require(rts)
  
  # earthengine_python <- Sys.getenv("EARTHENGINE_PYTHON", unset = NA)
  # Sys.setenv(RETICULATE_PYTHON = earthengine_python)
  # reticulate::py_config()
  #ee_clean_credentials() # Remove credentials if you want to change the user
  ee_Initialize(drive = TRUE) # initialize GEE, this will have you log in to Google Drive
  # ee_check() # Check non-R dependencies
  setwd(inp_wd)
  aoi <- st_read(aoi)

  bbox_aoi <- ee$Geometry$Rectangle(
    coords = st_bbox(aoi),
    proj = "EPSG:4326",
    geodesic = FALSE
  )

  #chirps
  ch_data <- ee$ImageCollection('UCSB-CHG/CHIRPS/DAILY')$filter(ee$Filter$date(st_date, en_date))
  data <- ch_data$select('precipitation')$toBands()
  ch_data <- ee_as_raster(data,
                          region = bbox_aoi,
                          via = "drive",
                          scale = 1000)

  time <- datevec <-
    seq(as.Date(st_date),
        by = "day",
        length.out = difftime(en_date, st_date, 'days', tz = "Africa/Addis_Ababa"))
  ch_rts <- rts(ch_data, time)
  new_crs <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
  mon_mean <- ch_rts %>% apply.monthly(mean) %>% apply.yearly(sum)
  mon_mean_ras <- mon_mean@raster %>% projectRaster(crs = new_crs, method = 'bilinear') %>%
    crop(aoi) %>% mask(aoi)

  writeRaster(mon_mean_ras, filename = "rainfall.tif", format = "GTiff", overwrite = TRUE)
}

inp_wd = "C:/Users/ATilaye/Documents/01My_Docs/01CIAT/01Landscape_Doctor/RScrips/geda_wshed"
st_date = "2010-01-01"
en_date = "2010-12-31"
aoi = "geda_kebele.shp"

calcfactorR(inp_wd, aoi, st_date, en_date)
