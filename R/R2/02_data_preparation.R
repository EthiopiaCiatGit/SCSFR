#------------------------------------------------------------------------------   
# Function to synchronize the geo-spatial layers and create probability layers
# ------------------------------------------------------------------------------

#' @param clim_path - path for the climate data
#' @param soil_path - path for soil data
#' @param topo_path - path for topographic data
#' @param aoi_path - path for a shapefile 
#' @param aoi- a shapefile for the study area 
#' @param dest_path - path for the output

data_Prep <- function(clim_path, soil_path, topo_path, aoi_path, aoi, dest_path)
{
  
  packages_required <- c("raster", "sf", "rgdal", "sp", "dplyr", "spatial.tools")
  
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  # load required packages
  invisible(lapply(packages_required, library, character.only = TRUE))
  
  setwd(aoi_path)
  aoi <- st_read(aoi)
  
  setwd(soil_path)
  soil_files <- list.files(path = ".", pattern = ".tif$", all.files = T)
  soil <- lapply(soil_files, raster)
  ocs <- soil[[4]]
  
  # check projection of the raster layers??
  
  setwd(clim_path)
  clim_files <- list.files(path = ".", pattern = ".tif$", all.files = T)
  climate <- lapply(clim_files, raster)
  
  for(i in 1:nlayers(climate)){
    climate[[i]] <- spatial_sync_raster(climate[[i]],ocs, method = "bilinear")%>%
                    crop(aoi) %>% mask(aoi)
  }
  
  setwd(topo_path)
  topo_files <- list.files(path = ".", pattern = ".tif$", all.files = T)
  topography <- lapply(topo_files, raster)

  for(i in 1:nlayers(topography)){
    topography[[i]] <- spatial_sync_raster(topography[[i]],ocs, method = "bilinear")%>%
      crop(aoi) %>% mask(aoi)
  }
  
  clim_below <- quantile(climate, probs = c(0.25))
  clim_normal <- quantile(climate, probs = c(0.5))
  clim_above <- quantile(climate, probs = c(0.75))
  
  stack_below <- stack(soil, topography, clim_below)
  stack_normal <- stack(soil, topography, clim_normal)
  stack_above <- stack(soil, topography, clim_above)
  
  setwd(dest_path)
  writeRaster(stack_below,
              filename = "stack_below.tif",
              format = "GTiff", overwrite = T)
 

  writeRaster(stack_normal,
              filename = "stack_normal.tif",
              format = "GTiff", overwrite = T)
 
  writeRaster(stack_above,
              filename = "stack_above.tif",
              format = "GTiff", overwrite = T)
  }

#------------------------------------------------------------------------------   
# Function to extract soil and topography data
# ------------------------------------------------------------------------------

#' @param stack_path - path for the topographic and raster paths
#' @param trial_path - path for trial/agronomic data
#' @param trial_data - trial/agronomic data
#' @param dest_path - path for the output

extract_SoilTopo <- function(stack_path, r_stack, trial_path, trial_data, dest_path)
{
  
  packages_required <- c("dplyr", "terra", "rgdal", "sp")
  
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  # load required packages
  invisible(lapply(packages_required, library, character.only = TRUE))
  
  setwd(stack_path)
  stack <- rast(stack)
  
  #select the layers that contain the soil and topographic variables
  soil_topo <- rast[[1:17]]
  
  #read trial data and change to spatial object
  trial <- read.csv(trial_data, header = T, sep = ",")
  pts <- vect(trial, geom = c("long", "lat"), crs = "epsg:4326")
  
  #extract the topographic and soil data using points
  soil_topo_df <- terra::extract(soil_topo, pts)
  soil_topo_df2 <- trial %>% select(id, long, lat) %>% cbind(soil_topo_df)
  
  setwd(dest_path)
  write.csv(clim_trial, "soil_topo.csv", col.names = T, row.names = F)
}

#-------------------------------------------------------------------------------   
# Function to bind the soil, topographic & climate data with trial data
# ------------------------------------------------------------------------------

#' @param clim_csv_path - path for the climate extracted data
#' @param trial_path - path for trial/agronomic data
#' @param clim_csv - path for topographic data
#' @param trial_data - trial/agronomic data
#' @param dest_path - path for the output

bind_csv_Data <- function(clim_csv_path, trial_path, clim_csv, trial_data, 
                          soil_topo_path, soil_topo_csv, dest_path)
{
  
  packages_required <- c("dplyr")
  
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  # load required packages
  invisible(lapply(packages_required, library, character.only = TRUE))
  
  setwd(clim_csv_path)
  clim_extr <- read.csv(clim_csv, header = T, sep = ",")
  
  setwd(soil_topo_path)
  soil_topo <- read.csv(soil_topo_csv, header = T, sep = ",")
  
  setwd(trial_path)
  trial <- read.csv(trial_data, header = T, sep = ",")
  clim_trial <- dplyr::inner_join(trial, soil_topo, clim_extr, by = c(id, long, lat))
  
  setwd(dest_path)
  write.csv(clim_trial, "covariates_trial.csv", col.names = T, row.names = F)
}
