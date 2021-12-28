# ------------------------------------------------------------------------------
# Data pre-processing for prediction

dataPreparation_Predict <-
  function(rast_path, csv, csv_path, aoi, aoi_path) {

  require(rgdal)
  require(raster)
  require(caret)
  require(dplyr)
  require(sf)
  require(spatial.tools)
  require(sp)
  
# ------------------------------------------------------------------------------ 
# list raster files
  dir.create("output", showWarnings = F)
  rfiles <- list.files(path = rast_path, pattern = ".tif$", all.files = T)
  rasters <- lapply(rfiles, raster)

# ------------------------------------------------------------------------------ 
# sync and stack rasters
  for(i in 1:length(rasters)){
    if("dem" %in% names(rasters[[i]])){ # setting dem as a reference
      dem <- rasters[[i]] 
    }
    }
  
  for(i in 1:length(rasters)){
    rasters[[i]] <-
      spatial_sync_raster(rasters[[i]], dem, method = 'bilinear')
  }
  grid <- stack(rasters)
  
# ------------------------------------------------------------------------------ 
# Mask rasters
  aoi <- readOGR(dsn = aoi_path, layer = aoi)
  new_crs <- proj4string(grid)
  aoi <- spTransform(aoi, new_crs)
  grid2 <- mask(grid, aoi)
  grid2 <- crop(grid, aoi)#removes NA values
  
  grid3 <- as(grid2, "SpatialGridDataFrame") #used later for prediction
  grid4 <- as.data.frame(grid3)
  grid4$landform <- as.factor(grid4$landform)

# ------------------------------------------------------------------------------
# rename columns of spatial grid data frame
  tmax_data_sgdf <- dplyr::select(grid4, contains("tmax"))
  for(i in 1:ncol(tmax_data_sgdf)){
    colnames(tmax_data_sgdf)[i] <- paste0("tmax", i)
  }
  
  tmin_data_sgdf <- dplyr::select(grid4, contains("tmin"))
  for(i in 1:ncol(tmin_data_sgdf)){
    colnames(tmin_data_sgdf)[i] <- paste0("tmin", i)
  }
  
  prcp_data_sgdf <- dplyr::select(grid4, contains("prcp"))
  for(i in 1:ncol(prcp_data_sgdf)){
    colnames(prcp_data_sgdf)[i] <- paste0("prcp", i)
  }
  
  srad_data_sgdf <- dplyr::select(grid4, contains("srad"))
  for(i in 1:ncol(srad_data_sgdf)){
    colnames(srad_data_sgdf)[i] <- paste0("srad", i)
  }
  
# ------------------------------------------------------------------------------ 
# Combine climatic data of new variable names
  clim_data_sgdf <-
    data.frame(tmax_data_sgdf,
               tmin_data_sgdf,
               prcp_data_sgdf,
               srad_data_sgdf)

# ------------------------------------------------------------------------------ 
# other than climatic data - sgdf
  other_data_sgdf <- grid4 %>% 
    dplyr::select(-(contains(c("prcp","tmin", "srad", "tmax")))) 
  
# ------------------------------------------------------------------------------
# Creating binary variables from factor variables - from landform
  bin_var <- model.matrix(as.formula(paste0("~","landform","+0")), other_data_sgdf)
  bin_var2 <- as.data.frame(bin_var)
  bin_var3 <- lapply(bin_var2, FUN = factor)
  bin_var_final <- as.data.frame(bin_var3)
  
  sgdf <- data.frame(other_data_sgdf, clim_data_sgdf, bin_var_final)# spatial grid data frame
  sgdf <- dplyr::select(sgdf, -"landform") # remove original landform
  
# ------------------------------------------------------------------------------
# excluding zero & near zero variance variables
  nzv <- nearZeroVar(sgdf, saveMetrics = TRUE) #df of predictor info returned
  sgdf <- sgdf[,!nzv[,4]]#nzv stored in the 4th column of nzv - remove

# ------------------------------------------------------------------------------
# Rearrange the columns
  sgdf_s1s2 <- dplyr::select(sgdf, "s1", "s2")
  sgdf_other <- dplyr::select(sgdf, -c("s1", "s2"))
  sgdf_final <- data.frame(sgdf_s1s2, sgdf_other)
  sgdf_final <- unique(na.omit(sgdf_final))
  save(sgdf_final, file = paste0("./output/","sgdf_final_output", ".rda"))

# ------------------------------------------------------------------------------
# Convert back to Spatial Grid Data frame
  gridded(sgdf_final) <- ~s1+s2
  proj4string(sgdf_final) <- new_crs
  
# ------------------------------------------------------------------------------
# read csv data and convert to spatial
  points <- read.csv(file = paste(csv_path,csv, sep = "/"), header = T)
  coordinates(points) <- ~x+y
  proj4string(points) <- new_crs

# ------------------------------------------------------------------------------
# extracting values by points
  grid_val <- over(points, sgdf_final) #overlay
  grid_val2 <- as.data.frame(grid_val)
  
# ------------------------------------------------------------------------------
# binding response variables and predictor data
  yield <- dplyr::select(as.data.frame(points), "yield")
  yield <- as.numeric(yield$yield)
  # save(yield, file = paste0("./output/","yield", ".rda"))
  
  data_pred_npks <- cbind(points@coords, yield, grid_val2)
  data_pred_npks <- unique(na.omit(data_pred_npks))
  covs_yield <- data_pred_npks[, 3:ncol(data_pred_npks)]#covariates & yield
  save(covs_yield, file = paste0("./output/","regression_matrix", ".rda"))
}

