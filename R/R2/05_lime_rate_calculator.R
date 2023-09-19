#------------------------------------------------------------------------------   
# Function to calculate lime requirement using different analytical methods
# ------------------------------------------------------------------------------

#' @param sou_path - the path of required raster files
#' @param dest_path - path for the final output


calculate_limeRate <- function(sou_path, dest_path){
    
    packages_required <- c("sf", "rgdal", "sp", "raster", "gdalUtils")
    
    # check and install packages that are not yet installed
    installed_packages <- packages_required %in% rownames(installed.packages())
    if(any(installed_packages == FALSE)){
      install.packages(packages_required[!installed_packages])}
    
    # load required packages
    invisible(lapply(packages_required, library, character.only = TRUE))
    
    #PH method
    cec <- raster("CEC_rf.tif")
    ph <- raster("pH_rf.tif")
    
    ph[ph > 6] <- NA
    
    #create a raster with target pH value
    target_ph <- raster(ph)
    values(target_ph) <- 6
    
    #use the formula lime = (target_ph - current_ph)*cec
    lime <- (target_ph - ph) * cec
    lime[lime <= 0] <- 0
  
    #acid saturation method
    al <- raster("al_0_20_isda.tif")
    al_mol <- (al)/10 #aluminium centimol per kg
    al <- al_mol %>% projectRaster(crs = crs(cec)) %>% 
      resample(cec, method = 'bilinear') %>% 
      crop(cec) %>% mask(cec)
    
    ecec <- raster("ecec_0-20cm_isda.tif") %>% projectRaster(crs = crs(cec)) %>% 
      resample(cec, method = 'bilinear') %>% 
      crop(cec) %>% mask(cec)
    
    lr <- 1160 * (al - (0.1*ecec)) * 0.001 # 0.001 is to convert ton/ha
    lr[lr <=0 ] <- 0
    
    #LiTAS model equation
    lr_litas <- (al - (0.1 * ecec)) / ((0.6 + 0.1) * (0.92 - 0.6))
    lr_litas[lr <= 0] <- 0
    plot(lr_litas)
    
    setwd(dest_path)
    writeRaster(lr, filename = "acid_sat_eth.tif", format = "GTiff") # kg/ha
    writeRaster(lime, filename = "ph_method.tif", format = "GTiff") # kg/ha
    writeRaster(lr_litas, filename = "lr_litas2.tif", format = "GTiff") # ton/ha
    
  }