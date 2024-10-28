library(terra)
library(tidyverse)

setwd("~/my_data/data/predicted_raster/optimal/")
n <- rast("maize_n_dominant.tif")
p <- rast("maize_p_dominant.tif")

comp <- (n * 100 / 0.8) / 1000
v_comp <- (n * 100 / 1) / 1000
nps <- p * 100 / (16.59)
urea <- (n - (nps * 19 / 100)) * 100 / 46

setwd("~/my_data/data/fertilizer/")
writeRaster(comp, file = "maize_compost.tif", filetype = "GTiff", overwrite = T)
writeRaster(v_comp, file = "maize_v_compost.tif", filetype = "GTiff", overwrite = T)
writeRaster(nps, file = "maize_nps.tif", filetype = "GTiff", overwrite = T)
writeRaster(urea, file = "maize_urea.tif", filetype = "GTiff", overwrite = T)
