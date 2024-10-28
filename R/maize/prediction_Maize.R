#prediction
library(terra)
library(rgdal)
library(sp)
library(caret)
library(ranger)

rm(list = ls())
woreda <- vect("~/my_data/data/woreda/woreda.shp")
setwd("~/my_data/data/scenario_raster/")
above <- rast("scenario_above.tif") |> 
  terra::crop(woreda) |>
  terra::mask(woreda)

#read final fitted model
setwd("../model/")
model <- readRDS("model_rf_DST_14vars.rds")

N <- c(seq(0, 140, 23)) 
P <- c(seq(0, 30, 10))

np <- expand.grid(N = N, P = P)
dim(np)

#create path and assocaited filename
pathOut <- "~/my_data/data/predicted_raster/above"

if(!dir.exists(pathOut)){
  suppressWarnings(dir.create(pathOut, recursive = T))
}

a <- apply(np, 1, function(i) paste(i, collapse = "."))
f <- file.path(pathOut, paste0("yield.", a, ".tif"))

#predict using the interaction of npk dataframe
for(i in 1:nrow(np)) {
  print(i)
  if (file.exists(f[i])) next
  NP <- data.frame(n_rate2 = np$N[i], p_rate2 = np$P[i])
  terra::predict(
    above,
    model,
    const = NP,
    filename = f[i],
    overwrite = TRUE,
    wopt = list(filetype = "GTiff", names = a[i], verbose = T),
    na.rm = T
  )
}
