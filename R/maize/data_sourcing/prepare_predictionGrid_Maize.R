#' -----------------------------------------------------------------------------
#' source all the generic data sourcing function and get soil, topography and 
#' climate data for the field trial points and ML prediction grid
#' 1. get point based climate + Topography + soil for the crop specified
#' 2. bind all the data and prepare the output for the data driven approach
#' 3. synchronize all the raster layers for the prediction grid in three 
#' different scenarios.
#' @param crops - the crop in which the data is generated
#' example - generate_AnalysisReadyData(crops = "Maize", useCaseName = NULL)
#' -----------------------------------------------------------------------------

generate_AnalysisReadyData <- function(crops, useCaseName){
  
  packages_required <- c("terra", "sf", "tidyverse", "readxl")
  
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  # load required packages
  invisible(lapply(packages_required, library, character.only = TRUE))
  
  pathIn <- paste("~/shared-data/Data", crops, "Intermediate/ML", sep = "/")
  
  
  #soil
  n <- rast(paste(path = pathIn, pattern = "N.tif", sep = "/"))
  # s <- rast(paste(path = pathIn, pattern = "S.tif", sep = "/")) |>
  #   crop(n) |> mask(n)
  # b <- rast(paste(path = pathIn, pattern = "B.tif", sep = "/")) |>
  #   crop(n) |> mask(n)
  # ca <- rast(paste(path = pathIn, pattern = "Ca_sat.tif", sep = "/")) |>
  #   crop(n) |> mask(n)
  # om <- rast(paste(path = pathIn, pattern = "OM.tif", sep = "/")) |> 
  #   terra::project(y = n)|>
  #   crop(n) |> mask(n)
  # 
  # soil_sp <- c(n, s, b, ca, om)
  
  tpi <- rast(paste(pathIn, "tpi.tif", sep = "/"))|> 
    terra::project(y = n) |>
    crop(n) |> mask(n)
  
  
  clim_norm_f <- rast()
  clim_norm <- list.files(path = pathIn, pattern = "_normal.tif", full.names = T)
  for(i in 1:length(clim_norm)){
    print(i)
    r <- rast(clim_norm[i]) |> terra::resample(n) |>
        terra::crop(n) |> terra::mask(n)
    terra::add(clim_norm_f) <- r
  }
  
  #topography
  # dem <- rast("~/shared-data/Data/General/Global_GeoData/Landing/ETH_dem.tif") |>
  #   terra::project(y = n) |> 
  #   terra::resample(n) |>
  #   crop(n) |> mask(n)
  
  sc_normal <- c(clim_norm_f[[c( 2, 5, 7, 8, 10, 11, 14, 15, 16, 17)]], tpi)
  nlyr(sc_normal)
  names(sc_normal)
  plot(sc_normal[[1]])
  
  names(sc_normal) <- c("Rainfall_month2", "totalRF", 
                        'relativeHumid_month2', 'relativeHumid_month3',
                        'solarRad_month2', 'solarRad_month3', 'Tmax_month3',
                        'Tmin_month1','Tmin_month2', 'Tmin_month3', "TPI")
  
  pathOut <- "~/my_data/data/scenario_raster/"
  writeRaster(sc_normal, paste(pathOut, "scenario_normal.tif", sep = "/"), 
              filetype = "GTiff", overwrite = T)
  
  clim_above_f <- rast()
  clim_above <- list.files(path = pathIn, pattern = "_above.tif", full.names = T)
  for(i in 1:length(clim_above)){
    print(i)
    r <- rast(clim_above[i]) |> terra::resample(n) |>
      terra::crop(n) |> terra::mask(n)
    terra::add(clim_above_f) <- r
  }
  
  sc_above <- c(clim_above_f[[c( 2, 5, 7, 8, 10, 11, 14, 15, 16, 17)]], tpi)
  nlyr(sc_above)
  plot(sc_above[[1]])
  
  names(sc_above) <- c("Rainfall_month2", "totalRF", 
                       'relativeHumid_month2', 'relativeHumid_month3',
                       'solarRad_month2', 'solarRad_month3', 'Tmax_month3',
                       'Tmin_month1','Tmin_month2', 'Tmin_month3', "TPI")
  writeRaster(sc_above, paste(pathOut, "scenario_above.tif", sep = "/"), 
              filetype = "Gtiff", overwrite = T)
  
  
  clim_below_f <- rast()
  clim_below <- list.files(path = pathIn, pattern = "_below.tif", full.names = T)
  for(i in 1:length(clim_below)){
    print(i)
    r <- rast(clim_below[i]) |> terra::resample(n) |>
      terra::crop(n) |> terra::mask(n)
    terra::add(clim_below_f) <- r
  }
  sc_below <- c(clim_below_f[[c( 2, 5, 7, 8, 10, 11, 14, 15, 16, 17)]], tpi)
  nlyr(sc_below)
  plot(sc_below[[1]])
  
  names(sc_below) <- c("Rainfall_month2", "totalRF", 
                       'relativeHumid_month2', 'relativeHumid_month3',
                       'solarRad_month2', 'solarRad_month3', 'Tmax_month3',
                       'Tmin_month1','Tmin_month2', 'Tmin_month3', "TPI")
  writeRaster(sc_below, paste(pathOut, "scenario_below.tif", sep = "/"), 
              filetype = "Gtiff", overwrite = T)
}

