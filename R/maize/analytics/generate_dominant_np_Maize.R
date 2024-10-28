packages_required <- c("raster","dplyr", "sf",  "rgdal", "terra")

# check and install packages that are not yet installed
installed_packages <- packages_required %in% rownames(installed.packages())
if(any(installed_packages == FALSE)){
  install.packages(packages_required[!installed_packages])}

# load required packages
invisible(lapply(packages_required, library, character.only = TRUE))

rm(list = ls())
setwd("~/my_data/data/predicted_raster/optimal/")
nutrient <- list.files(path = ".", pattern = "maize_p_", all.files = T)
nutrient <- terra::rast(nutrient)
names(nutrient) <- c("above", "below", "normal")

setwd("~/my_data/data/climate_dominant/")
dominant <- terra::rast("dominant.tif") |> 
  terra::resample(nutrient[[1]], method = "near") |>
  terra::crop(nutrient[[1]]) |>
  terra::mask(nutrient[[1]])

names(dominant) <- "dominant" 
#create a stack of all layers
rstack <- c(nutrient, dominant)
class(rstack)
nlyr(rstack)
plot(rstack)
df <- as.data.frame(rstack, xy = T, na.rm = F)
colnames(df)

df2 <- df %>% mutate(domnt = case_when(dominant == 1 & !is.na(dominant) ~ below,
                                      dominant == 2 & !is.na(dominant) ~ normal,
                                      dominant == 3 & !is.na(dominant) ~ above))

p_dominant <- rast(nutrient[[1]])
values(p_dominant) <- df2$domnt
setwd("~/my_data/data/predicted_raster/optimal")
writeRaster(p_dominant, filename = "maize_p_dominant.tif", filetype = "GTiff", overwrite = T)
