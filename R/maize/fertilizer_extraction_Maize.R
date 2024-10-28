rm(list = ls())
library(terra)
library(dplyr)

setwd("~/my_data/data/fertilizer/")
nps <- rast("maize_nps.tif")
urea <- rast("maize_urea.tif")
comp <- rast("maize_compost.tif")
vcomp <- rast("maize_v_compost.tif")

stack <- c(nps, urea)
names(stack) <- c("nps", "urea")
plot(stack)

setwd("~/my_data/data/farmer_gps/")
data <- read.csv("gps_reading_mareko_meskan.csv", header = T, sep = ",")
pts <- vect(data, geom=c("X_GPS.coordinates_longitude",
                         "X_GPS.coordinates_latitude"), 
            crs="epsg:4326")
#pts2 <- terra::project(pts,y="epsg:4326")

pts_value <- extract(stack, pts) %>% 
  dplyr::select(-c("ID")) |> round(0)

colnames(data)
final_data <- cbind(data[,c(2:5, 9, 10)], pts_value)
colnames(final_data) <- c("woreda", "kebele", "got", "farmer" ,"latitude", "longitude",
                         "nps_kgpha", "urea_kgpha")
write.csv(final_data, "mareko_meskan_advisory.csv", col.names = T, row.names = F)
