
#generate dominant fertilizer based on climate dominant probability
library(raster)
library(rgdal)
library(maptools)
library(dplyr)

rm(list = ls())
setwd("C:\\Users\\ATilaye\\Documents\\01My_Docs\\01CIAT\\02Fertilizer\\prediction_2023/data/final/below")
list.files()
p_below <- raster("wheat_p_below.tif")
setwd("../above/")
list.files()
p_above <- raster("wheat_p_above.tif")
setwd("../normal/")
p_normal <- raster("wheat_p_normal.tif")
setwd("../dominant/input/")
dominant <- raster("aclimate_et_seasonal_country_et_dominant.tif") %>% resample(p_normal, method = 'ngb') %>%
  crop(p_normal) %>% mask(p_normal)

rstack <- stack(p_above,p_below,p_normal, dominant)

names(rstack) <- c("p_above","p_below","p_normal", "dominant")
df <- as.data.frame(rstack, xy = T)

df2 <- df %>% mutate(g = case_when(dominant >= 0 & dominant < 100 & !is.na(dominant) ~ p_above,
                                   dominant >= 100 & dominant < 200 & !is.na(dominant) ~ p_normal,
                                   dominant >= 200 & !is.na(dominant) ~ p_below))

p_dominant <- raster(p_normal)
values(p_dominant) <- df2$g
setwd("../final/")
writeRaster(p_dominant, filename = "wheat_p_dominant_ngb", format = "GTiff", overwrite = T)
#values(k_layer) <- nps_df$k

