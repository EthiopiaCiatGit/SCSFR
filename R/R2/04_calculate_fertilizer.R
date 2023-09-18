#-------------------------------------------------------------------------------   
# Function to generate optimal nuitrients from the yield
# ------------------------------------------------------------------------------

#' @param layer_path - path for the predicted layers
#' @param dest_path - path for the output

generate_optimalNutrient <- function(layer_path, dest_path){
  
  packages_required <- c("raster","dplyr", "sf",  "rgdal")
  
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  # load required packages
  invisible(lapply(packages_required, library, character.only = TRUE))
  
  # read predicted layers and calculate optimal yield
  setwd(layer_path)
  rfiles <- list.files(path = ".", pattern = ".tif$", all.files = T)
  yield_ras <- lapply(rfiles, raster::raster)
  yield_stack <- stack(yield_ras)
  optimal_yield <- max(yield_stack)
  
  setwd(dest_path)
  writeRaster(
    optimal_yield,
    filename = "wheat_optimal_yield.tif",
    format = "GTiff",
    overwrite = TRUE)
  
  # generate n, p, k layers that generates the optimal yield
  grid <- as(optimal_yield, "SpatialGridDataFrame")
  grid2 <- as.data.frame(grid)
  vals <- values(yield_stack) # creates a matrix
  vals[is.na(vals)] <- -1 # changing NA values to -1
  col_vals <- colnames(vals)[max.col(vals, ties.method = "first")]
  
  # capture layer names and change to data frame
  col_vals <- as.data.frame(col_vals)
  names(col_vals) <- "yield"
  nps_df <- col_vals %>% 
    separate(yield, into = paste0('yield', 1:3), sep = '[+.]') # yield - name of the column to separate
  names(nps_df) <- c("yield", "n", "p")
  nps_df$n <- as.integer(nps_df$n)
  nps_df$p <- as.integer(nps_df$p)
  #nps_df$k <- as.integer(nps_df$k)

  # create a csv file for the layers
  n_layer <- p_layer <- raster(optimal_yield)
  values(n_layer) <- nps_df$n
  values(p_layer) <- nps_df$p
  #values(k_layer) <- nps_df$k
  
  n_layer2 <- mask(n_layer, optimal_yield)
  p_layer2 <- mask(p_layer, optimal_yield)
  #k_layer2 <- mask(k_layer, optimal_yield)
  
  writeRaster(n_layer2, filename = "wheat_n_normal", format = "GTiff", overwrite = T)
  writeRaster(p_layer2, filename = "wheat_p_normal", format = "GTiff", overwrite = T)
}


#-------------------------------------------------------------------------------   
# Function to generate dominant nutrients scenario from EDACaP layer
# ------------------------------------------------------------------------------

#' @param nuitrient_path - path for the nutrients
#' @param dominant_prob_path - path for the climate scenario
#' @param dominant_prob - dominant probability raster
#' @param dest_path - destination path for output

generate_dominantNutrient <- function(nuitrient_path, dest_path,
                                      dominant_prob_path, dominant_prob){
  
  packages_required <- c("raster","dplyr", "sf",  "rgdal")
  
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  # load required packages
  invisible(lapply(packages_required, library, character.only = TRUE))
  
  setwd(nuitrient_path)
  nutrient <- list.files(path = ".", pattern = ".tif$", all.files = T)
  nutrient <- lapply(nutrient, raster)
  names(nutrient) <- c("below", "above", "normal")
  
  setwd(dominant_prob_path)
  dominant <- raster(dominant_prob)
  
  #create a stack of all layers
  rstack <- stack(nutrient, dominant)
  
  df <- as.data.frame(rstack, xy = T)
  
  df2 <- df %>% mutate(dmnt = case_when(dominant >= 0 & dominant < 100 & !is.na(dominant) ~ above,
                                     dominant >= 100 & dominant < 200 & !is.na(dominant) ~ normal,
                                     dominant >= 200 & !is.na(dominant) ~ below))
  
  dominant <- raster(nutrient[[1]])
  values(p_dominant) <- df2$dmnt
  setwd(dest_path)
  writeRaster(p_dominant, filename = "wheat_p_dominant_ngb", format = "GTiff", overwrite = T)
}


#-------------------------------------------------------------------------------
#Script to calculate nps, urea, compost and vermi compost
#-------------------------------------------------------------------------------

#' @param sou_path - path for the nutrients
#' @param dominant_n - path for the climate scenario
#' @param dominant_p - dominant probability raster
#' @param dest_path - destination path for output

calculate_fertilizer <- function(sou_path, dominant_n, dominant_p, 
                            dest_path){
  
  require(raster)
  
  setwd(sou_path)
  n <- raster(dominant_n)
  p <- raster(dominant_p)
  comp <- (n * 100 / 0.8) / 1000
  v_comp <- (n * 100 / 1) / 1000
  nps <- p * 100 / (16.59)
  urea <- (n - (nps * 19 / 100)) * 100 / 46
  
  setwd(dest_path)
  writeRaster(comp, file = "wheat_compost.tif", format = "GTiff")
  writeRaster(v_comp, file = "wheat_v_compost.tif", format = "GTiff")
  writeRaster(nps, file = "wheat_nps.tif", format = "GTiff")
  writeRaster(urea, file = "wheat_urea.tif", format = "GTiff")
}


#-------------------------------------------------------------------------------
#Script to to extract the fertilizer rates based on point data
#-------------------------------------------------------------------------------

#' @param sou_path - path for the nutrients
#' @param csv_path - path for points data
#' @param csv - points data
#' @param dest_path - destination path for output

extract_fertilizerPoint <- function(sou_path, csv_path, csv, 
                                 dest_path){
  require(terra)
  require(dplyr)
  
  setwd(sou_path)
  r <- list.files(pattern = ".tif$", all.files = T) %>% lapply(rast)
  stack <- c(r) %>% rast()
  
  setwd(csv_path)
  data <- read.csv(csv, header = T, sep = ",")
  pts <- vect(data, geom = c("lon", "lat"), crs = "epsg:4326")
  
  pts_value <- extract(stack, pts) %>% dplyr::select(-c("ID"))
  
  final_data <- cbind(data, pts_value)
  write.csv(final_data, "point_based_recomm.csv", col.names = T)
}
  

#-------------------------------------------------------------------------------   
# Function to generate advisory based on shapefile (kebele or woreda)
# ------------------------------------------------------------------------------

#' @param sou_path - path for the fertilizers
#' @param shp_path - path for the shapefile
#' @param shp - shapefile
#' @param dest_path - destination path for output

generate_recommShape <- function(sou_path, shp_path, shp,
                                      dest_path){
  
  packages_required <- c("terra","dplyr", "sf",  "rgdal")
  
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  # load required packages
  invisible(lapply(packages_required, library, character.only = TRUE))
  
  setwd(sou_path)
  r <- list.files(pattern = ".tif$", all.files = T) %>% lapply(rast)
  stack <- c(r) %>% rast()
  names(stack) <- c("wheat_compost", "wheat_n", "wheat_nps", "wheat_optimal_yield",
                    "wheat_p", "wheat_urea", "wheat_vcompost")
  setwd(shp_path)
  aoi <- vect(shp)
  aoi_val <- extract(stack, kebele, fun = mean, na.rm = T)
  
  for(j in 2:ncol(keb_val)){
    final_val <- cbind(kebele, as.data.frame(keb_val[,j]))
    names(final_val)[ncol(final_val)] <- "value"
    writeVector(final_val, filename = paste0(names(stack[[j-1]]), ".shp"), filetype = "ESRI Shapefile", overwrite = T)
  }
}
