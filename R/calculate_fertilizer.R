#-------------------------------------------------------------------------------
#Script to calculate nps, urea, compost and vermi compost
#-------------------------------------------------------------------------------

calc_nps_urea <- function(s_path, n, p, d_path){
  
    setwd(s_path)
    n <- raster(n)
    p <- raster(p)
    comp <- (n*100/0.8)/1000
    v_comp <- (n*100/1)/1000
    nps <- p*100/(16.59)
    urea <- (n-(nps*19/100))*100/46
    setwd(d_path)
    writeRaster(comp, file = "wheat_compost.tif", format = "GTiff")
    writeRaster(v_comp, file = "wheat_v_compost.tif", format = "GTiff")
    writeRaster(nps, file = "wheat_nps.tif", format = "GTiff")
    writeRaster(urea, file = "wheat_urea.tif", format = "GTiff")
    
}


