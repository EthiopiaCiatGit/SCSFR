#------------------------------------------------------------------------------   
# Function to access and download ISRIC Soil grids. Currently we have access to 
# EthioSIS soil parameter layers and we access them locally since they are not 
# publicly available freely (check https://nsis.moa.gov.et/geoportal/).
# ------------------------------------------------------------------------------

#' @param voi - variable of interest to be downloaded (bdod, clay, sand, ...)
#' @param depth - depth of soil layer (0-5cm, 5-15cm ...)
#' @param quantile - the quantiles (mean, q0.5 ... )
#' @param aoi - a shapefile for the study area 


fetchSoilGrid <-
  function(voi, depth, quantile, aoi){
   
    packages_required <- c("XML", "sf", "rgdal", "sp", "raster", "gdalUtils")
    
    # check and install packages that are not yet installed
    installed_packages <- packages_required %in% rownames(installed.packages())
    if(any(installed_packages == FALSE)){
      install.packages(packages_required[!installed_packages])}
    
    # load required packages
    invisible(lapply(packages_required, library, character.only = TRUE))
    
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
    
    # ------------------------------------------------------------------------------   
    # layer of interest
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


# ------------------------------------------------------------------------------
# converting list of SoilGrid geo-tiff files using conversion factors
# and calculating weighted average
# ------------------------------------------------------------------------------

#' @param sou_path - the path of downloaded ISRIC data
#' @param dest_path - path for the final output

convertSoilGrids <- function(sou_path, dest_path = NULL){
  
  packages_required <- c("raster", "rgdal", "dplyr")
  
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  # load required packages
  invisible(lapply(packages_required, library, character.only = TRUE))
  
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


# ------------------------------------------------------------------------------
# Fetches elevation data from Amazon Web Services using and area of interest. 
# Currently we are using SRTM DEM downloaded and stored in CGLabs.
# ------------------------------------------------------------------------------ 

#' @param aoi - a shapefile for the study area
#' @param dest_path - path for the final output

fetchElevation <- function(aoi, aoi_path = NULL){
  
  packages_required <- c("raster", "rgdal", "elevatr", "sf", "rgl")
  
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  # load required packages
  invisible(lapply(packages_required, library, character.only = TRUE))
  
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
  dem <- get_elev_raster(aoi, z = 9)
  
  # ------------------------------------------------------------------------------
  # derive slope and aspect
  slope <- terrain(dem, opt = 'slope', unit = 'degrees')
  aspect <- terrain(dem, opt = 'aspect', unit = 'degrees')
  tpi <- terrain(dem, opt = 'TPI')
  tri <- terrain(dem, opt = 'TRI')
  SD <- sd(tpi[], na.rm = T)
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

# --------------------------------------------------------------------
# A function for downloading chirps precipitation data from GEE
# --------------------------------------------------------------------

#' @param inp_wd - input working directory
#' @param aoi - a shapefile for the study area
#' @param st_date - start date
#' @param en_date - end date

fetchClimData <- function(inp_wd, aoi, st_date, en_date){
  
  packages_required <- c("rgee", "reticulate", "sf", "terra", "rts")
  
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  # load required packages
  invisible(lapply(packages_required, library, character.only = TRUE))
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
  mon_mean <- ch_rts %>% apply.monthly(mean)
  mon_mean_ras <- mon_mean@raster %>% projectRaster(crs = new_crs, method = 'bilinear') %>%
    crop(aoi) %>% mask(aoi)
  
  writeRaster(mon_mean_ras, filename = "rainfall.tif", format = "GTiff", overwrite = TRUE)
}

# ------------------------------------------------------------------------------
#' Extraction of climate monthly data (precipitaion - CHIRPS & others from AgEra5) 
#' using field trial data based on different crops. 
#' The climate variables in this extraction are precipitation, temp max,
#' temp min, solar radiation and relative humidity. 
#' @description - this function selects a crop from trial data, extracts climate 
#' data and finally writes the results for the crop specified. Besides it 
#' generates a raster stack for three different scenarios(above, below & normal)
#' for the ML prediction   
#' @param climate the climate variable to be extracted should be one of 
#' precipitation, tempMax, tempMin, solaRadiation & relativeHumidity.
#' @param crops the name of the crop is one of  maize, sorghum teff or wheat
#' @example - get_geoSpatialClimate(climate = "precipitation", crops = "Maize")
# ------------------------------------------------------------------------------

get_geoSpatialClimate <- function(climate, crops){
  # install & load packages
  packages_required <- c("terra", "tidyverse", "sf", "lubridate", "plyr")
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  # load required packages
  invisible(lapply(packages_required, library, character.only = TRUE))
  
  # read trial sample data
  pathIn_trial <- "~/Eth_DST_Harmonization/data/national_Data/raw/field_Data"
  trial <- readxl::read_xlsx(paste(pathIn_trial, "trial_data.xlsx", sep = "/"), 
                             col_names = T) |> dplyr::mutate(crop = tolower(crop)) |>
    dplyr::filter(crop == tolower(crops)) |> na.omit()
  
  # change trial data to spatial points
  points <- vect(trial, geom = c("lon", "lat"), crs = "epsg:4326")
  
  # read sowing date and growing length decade data
  pathIn_sow <- "~/Eth_DST_Harmonization/data/national_Data/raw/geospatial/sowing_Date"
  pathIn_gleng <- "~/Eth_DST_Harmonization/data/national_Data/raw/geospatial/growing_Length"
  sow_start <- rast(paste(pathIn_sow, paste0(tolower(crops), ".tif"), sep = "/")) * 10
  sow_start[sow_start <= 0 | sow_start > 366] <- NA
  grow_length <- rast(paste(pathIn_gleng, paste0(tolower(crops), ".tif"), sep = "/")) * 10
  grow_length[grow_length <= 0 | grow_length > 365] <- NA
  
  #mechanically set the sowing date and harvesting date 
  #values(sow_start) <- 91 #start of April 
  #values(grow_length) <- 210 #end of October for maize physiologic maturity
  
  # extract the sowing start and growing length by the trial points data
  start <- sow_start |> terra::extract(points)
  grow_len <- grow_length |> terra::extract(points)
  
  grow_period <- cbind(start[, -1], grow_len[, -1]) 
  colnames(grow_period) <- c("sowing_start", "grow_length")
  
  # bind the trial data with sowing date and growing length
  trial_with_grow_len <- cbind(trial, grow_period)
  trial_with_grow_len <- trial_with_grow_len |> 
    dplyr::mutate(full_date = as.Date(paste(year, "01", "01", sep = "-"))) |>
    dplyr::mutate(pl_date = ymd(full_date) + days(sowing_start)) |>
    dplyr::mutate(hv_date = ymd(full_date) + days(sowing_start+grow_length)) |> 
    na.omit()
  
  # read the climate monthly data from different directories
  if(climate == "precipitation"){
    path <- "~/Eth_DST_Harmonization/data/national_Data/raw/geospatial/climate/rainfall/monthly"
    path_daily <- "~/Eth_DST_Harmonization/data/national_Data/raw/geospatial/climate/rainfall/daily"
  }else if(climate == "tempMax"){
    path <- "~/Eth_DST_Harmonization/data/national_Data/raw/geospatial/climate/temp_Max/monthly"
  }else if(climate == "tempMin"){
    path <- "~/Eth_DST_Harmonization/data/national_Data/raw/geospatial/climate/temp_Min/monthly"
  }else if(climate == "solarRadiation"){
    path <- "~/Eth_DST_Harmonization/data/national_Data/raw/geospatial/climate/solar_Radiation/monthly"
  }else if(climate == "relativeHumidity"){
    path <- "~/Eth_DST_Harmonization/data/national_Data/raw/geospatial/climate/relative_Humidity/monthly/"
  }else{
    print("Please enter a valid climate variable name")
  }
  
  unique_grow_len <- trial_with_grow_len |> dplyr::select(year, pl_date, hv_date) |>
    unique()
  
  # loop unique planting and harvesting date to extract the climate data
  f_df <- c()
  if(climate != "precipitation"){
    for(i in 1:nrow(unique_grow_len)){
      print(i)
      d <- trial_with_grow_len |> dplyr::filter(pl_date == unique_grow_len$pl_date[i] & hv_date == unique_grow_len$hv_date[i])
      pl_date <- d$pl_date[1]
      hv_date <- d$hv_date[1]
      year <- d$year[1]
      d <- vect(d, geom = c("lon", "lat"), crs = "epsg:4326")
      if(pl_date < hv_date){
        r <- rast(paste(path, paste0(year,".tif"), sep = "/"))
        month_start <- as.numeric(format(pl_date, "%m"))
        month_end <- as.numeric(format(hv_date, "%m"))
        grow_len_r <- r[[month_start:month_end]]
      }else if(pl_date > hv_date){
        month_start <- as.numeric(format(pl_date, "%m"))
        month_end <- as.numeric(format(hv_date, "%m"))
        r1 <- rast(paste(path, paste0(year,".tif"), sep = "/"))
        r2 <- rast(paste(path, paste0((year)+1, ".tif"), sep = "/"))
        r <- c(r1,r2)
        grow_len_r1 <- r[[month_start:12]]
        grow_len_r2 <- r[[13:month_end + 12]]
        grow_len_r <- c(grow_len_r1, grow_len_r2)
      }
      r_df <- terra::extract(grow_len_r, d) |> select(-ID)
      colnames(r_df) <- paste0(climate, "_g_len_", seq(1, ncol(r_df),1))
      f_df <- rbind.fill(f_df, r_df) # row binds df with different column numbers
    }
  }else if(climate == "precipitation"){
    for(i in 1:nrow(unique_grow_len)){
      print(i)
      d <- trial_with_grow_len |> dplyr::filter(pl_date == unique_grow_len$pl_date[i] & hv_date == unique_grow_len$hv_date[i])
      pl_date <- d$pl_date[1]
      hv_date <- d$hv_date[1]
      pl_date_day <- yday(d$pl_date[1]) 
      hv_date_day <- yday(d$hv_date[1])
      year <- d$year[1]
      d <- vect(d, geom = c("lon", "lat"), crs = "epsg:4326")
      if(pl_date < hv_date){
        r_mon <- rast(paste(path, paste0(year,".tif"), sep = "/"))
        r_daily <- rast(paste(path_daily, paste0(year,".tif"), sep = "/"))
        month_start <- as.numeric(format(pl_date, "%m"))
        month_end <- as.numeric(format(hv_date, "%m"))
        grow_len_r_mon <- r_mon[[month_start:month_end]]
        grow_len_r_daily <- r_daily[[pl_date_day:hv_date_day]]
      }else if(pl_date > hv_date){
        month_start <- as.numeric(format(pl_date, "%m"))
        month_end <- as.numeric(format(hv_date, "%m"))
        r1 <- rast(paste(path, paste0(year,".tif"), sep = "/"))
        r2 <- rast(paste(path, paste0((year)+1, ".tif"), sep = "/"))
        r <- c(r1,r2)
        grow_len_r1 <- r[[month_start:12]]
        grow_len_r2 <- r[[13:month_end + 12]]
        grow_len_r <- c(grow_len_r1, grow_len_r2)
        r1_day <- rast(paste(path_mon, paste0(year,".tif"), sep = "/"))
        r2_day <- rast(paste(path_mon, paste0((year)+1, ".tif"), sep = "/"))
        r_daily <- c(r1_day, r2_day)
        grow_len_r1_daily <- r_daily[[pl_date_day:365]]
        grow_len_r2_daily <- r_daily[[366:365+hv_date_day]]
        grow_len_r_daily <- c(grow_len_r1_daily, grow_len_r2_daily)
      }
      #total rainfall for the growing period
      totalRF <- grow_len_r_daily |> terra::app(fun = sum, na.rm = T)
      names(totalRF) <- "totalRF"
      totalRF_d <- terra::extract(totalRF, d) |> 
        dplyr::select(-ID)
      #number of rainy days for the growing period
      nRainyDays <- terra::extract(grow_len_r_daily, d) |> 
        dplyr::select(-ID)
      nRainyDays_t <- t(nRainyDays)
      nRainyDays_t[nRainyDays_t < 2] <- 0
      nRainyDays_t[nRainyDays_t >= 2] <- 1
      if(nrow(d) <= 1){
        nRainyDays_t2 <- sum(nRainyDays_t[, 1], na.rm = T)
      }else{
        nRainyDays_t2 <- colSums(nRainyDays_t[1:nrow(nRainyDays_t),], na.rm = T)
      }
      nrRainyDays <- as.data.frame(nRainyDays_t2)
      colnames(nrRainyDays) <- "nrRainyDays"
      r_df <- terra::extract(grow_len_r_mon, d) |> select(-ID)
      colnames(r_df) <- paste0(climate, "_g_len_", seq(1,ncol(r_df),1))
      final_df <- cbind(r_df, nrRainyDays, totalRF_d)
      f_df <- rbind.fill(f_df, final_df) # row binds df with different column numbers
    }
  }
  
  # convert to appropriate units
  if(climate %in% c("tempMax","tempMin")){
    f_df <- f_df - 274
  }else if(climate == "solarRadiation"){
    f_df <- f_df / 1000000
  }
  
  # fill na values with row means - with different growing lengths some areas have 3,
  # others have 4, others have 5. The shorter ones have NA values to match the larger
  # ones. The NA values will be filled by row means
  
  # k <- which(is.na(f_df), arr.ind = TRUE)
  # f_df[k] <- rowMeans(f_df, na.rm = TRUE)[k[, 1]]
  
  final_df <- trial_with_grow_len |> dplyr::select(id, lon, lat, year, pl_date, hv_date) |>
    cbind(f_df)
  
  # ----------------------------------------------------------------------------
  # create a prediction scenario raster
  min_pl_date <- min(format(as.Date(final_df$pl_date,format="%d/%m/%Y"),"%m")) |>
    as.numeric()
  max_hv_date <- max(format(as.Date(final_df$hv_date,format="%d/%m/%Y"),"%m")) |>
    as.numeric()
  clim_rast <- list.files(path = path, pattern = "*.tif$", full.names = T) |> 
    rast()
  
  q1_clim <- rast()
  q2_clim <- rast()
  q3_clim <- rast()
  
  for (i in 1:12){
    clim_year <- rast()
    k <- i
    for(j in 1:(nlyr(clim_rast)/12)){
      print(k)
      add(clim_year) <- clim_rast[[k]]
      k <- k + 12
    }
    q1 <- terra::quantile(clim_year, probs = 0.25)
    add(q1_clim) <- q1
    q2 <- terra::quantile(clim_year, probs = 0.5)
    add(q2_clim) <- q2
    q3 <- terra::quantile(clim_year, probs = 0.75)
    add(q3_clim) <- q3
  }
  
  clim_below <- q1_clim[[min_pl_date:max_hv_date]]
  clim_normal <- q2_clim[[min_pl_date:max_hv_date]]
  clim_above <- q3_clim[[min_pl_date:max_hv_date]]
  names(clim_above) <- names(clim_normal) <- names(clim_above) <-
    paste0(climate, "_g_len_", seq(1,nlyr(clim_below),1))
  
  # convert to appropriate units
  if(climate %in% c("tempMax","tempMin")){
    clim_above <- clim_above - 274
    clim_normal <- clim_normal - 274
    clim_below <- clim_below - 274
  }else if (climate == "solarRadiation"){
    clim_above <- clim_above / 1000000
    clim_normal <- clim_normal / 1000000
    clim_below <- clim_below / 1000000
  }
  
  # mask using the growing area of the crop
  path_CropMask <- "~/Eth_DST_Harmonization/data/national_Data/raw/geospatial/crop_Mask/"
  crop_mask <- terra::rast(paste(path_CropMask, paste0(tolower(crops), ".tif"), sep = "/"))
  crop_mask <- crop_mask |> terra::resample(clim_above)
  clim_above <- clim_above |> terra::crop(crop_mask) |> terra::mask(crop_mask)
  clim_normal <- clim_normal |> terra::crop(crop_mask) |> terra::mask(crop_mask)
  clim_below <- clim_below |> terra::crop(crop_mask) |> terra::mask(crop_mask)
  
  # ----------------------------------------------------------------------------
  # export the results to the crop path  
  pathOut_gps <- paste("~/Eth_DST_Harmonization/data/data_Processing", 
                       crops, "geospatial/data_Driven/gps", sep = "/")
  if(!exists(pathOut_gps)){
    suppressWarnings(dir.create(pathOut_gps, recursive = T))
  }
  
  pathOut_layers <- paste("~/Eth_DST_Harmonization/data/data_Processing", 
                          crops, "geospatial/data_Driven/climate", sep = "/")
  if(!exists(pathOut_layers)){
    suppressWarnings(dir.create(pathOut_layers, recursive = T))
  }
  
  final_df <- final_df |> select(-c(hv_date, pl_date))
  saveRDS(final_df, paste(pathOut_gps, paste0(climate, ".rds"), sep = "/"))
  
  terra::writeRaster(
    clim_above,
    filename = paste(pathOut_layers,paste0(climate, "_above", ".tif"), sep = "/"),
    filetype = "GTiff",
    overwrite = T
  )
  terra::writeRaster(
    clim_normal,
    filename = paste(pathOut_layers,paste0(climate, "_normal", ".tif"), sep = "/"),
    filetype = "GTiff",
    overwrite = T
  )
  terra::writeRaster(
    clim_below,
    filename = paste(pathOut_layers,paste0(climate, "_below", ".tif"), sep = "/"),
    filetype = "GTiff",
    overwrite = T
  )
}


# --------------------------------------------------------------------
# A function for accessing trial data
# --------------------------------------------------------------------

#' @param inp_wd - input working directory
#' @param aoi - a shapefile for the study area
#' @param st_date - start date
#' @param en_date - end date
#' 
get_trialData <- function(dest_path = NULL){
  
  remotes::install_github("iqss/dataverse-client-r")
  packages_required <- c("dataverse", "dplyr", "stringr", "ggplot2", "ggpubr", "tibble",
                         "reshape", "sp", "terra", "readxl")
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  suppressWarnings(suppressPackageStartupMessages(invisible(lapply(packages_required, library, character.only = TRUE))))
  
  #Read the data from dataverse
  data <- 
    get_dataframe_by_name(
      filename   = "Fertilizer data.tab",
      dataset    = "doi:10.7910/DVN/M8FCSL",
      server     = "dataverse.harvard.edu",
      original   = TRUE,
      .f         = readxl::read_xlsx) %>% as.data.frame() %>%
    select(id, long, lat, crop, year, n, p, k, grain_yield)
  
  
  setwd(dest_path)
  write.csv(data, "trial_data.csv", col.names = T, row.names = F)
}

# --------------------------------------------------------------------
# A function for extracting monthly clim data based on trial data
# --------------------------------------------------------------------

#' @param inp_wd - input working directory
#' @param clim_data - climate data
#' @param trial_data - trial_data
#' @param inp_wd - destination working directory


extract_Climate <- function(inp_wd, clim_data, dest_wd, trial_data){
  
  packages_required <- c("sp", "rgdal", "sf", "terra", "dplyr", "doParallel")
  
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  #parallelize
  rm(list = ls())
  cl <- makeCluster(detectCores()/3, type='PSOCK')
  registerDoParallel(cl)
  
  setwd(inp_wd)
  clim <- rast(clim_data)
  
  df <- read.csv(trial_data, header = T, sep = ",") %>% 
    select(id, lon,lat, year) %>% filter(year >= 1981)
  
  f_df <- NULL
  for(i in 1:nrow(df)){
    print(paste0("row_",i))
    r_df <- NULL
    for(j in 1:nlyr(clim)){
      z <- time(clim[[j]])
      x <- as.numeric(format(z, "%Y"))
      #df_t <- as.numeric(df$year[i])
      if( x == as.numeric(df$year[i])){
        r <- subset(clim[[j]], 1)
        pt <- vect(df[i,], geom=c("lon", "lat"), crs="epsg:4326")
        r_df1 <- extract(r, pt, xy = T)
        colnames(r_df1) <- c("ID",paste0("clim", i))
        r_df <- rbind(r_df, r_df1)
      }
    }
    #selects the layers 6:8 and transposes for the months
    t <- as.data.frame(t(r_df[5:10,2])) %>% cbind(df[i,])
    f_df <- rbind(t,f_df)
  }
  
  colnames(f_df)[1:6] <- paste0("clim", seq(1, 6, 1))
  setwd(dest_wd)
  write.csv(f_df, "clim_input.csv", col.names = T, row.names = F)
  stopCluster(cl)
}

