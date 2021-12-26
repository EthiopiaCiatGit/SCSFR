# ------------------------------------------------------------------------------
# Reads rda output to predict optimal yield
predictFertilizer <- function(rda_file_path, nuitrent){
  require(rgdal)
  require(raster)
  require(caret)
  require(ranger)
  require(dplyr)
  require(e1071)
  require(randomForest)
  require(sp)
  require(hydroGOF)
  require(Metrics)
  require(ggplot2)
  
  dir.create("predict_output", showWarnings = F)
  
# ------------------------------------------------------------------------------
# loading the rda outputs and adding to the global env
  rda_files <- list.files(path = rda_file_path, pattern = ".rda")
  lapply(rda_files,load,.GlobalEnv)
  
  new_crs <- "+proj=longlat +datum=WGS84 +no_defs"
  
# ------------------------------------------------------------------------------
# Create training and testing data set
  mtry <- as.integer((ncol(covs_yield))/3)
  mtry <- seq(mtry-10, mtry+10, by = 3)

# ------------------------------------------------------------------------------ 
# model calibration and validation
  set.seed(7)
  
  
  rf_fitControl <- trainControl(method = "repeatedcv",
                             number = 5,
                             repeats = 3)
  
  rf_tuneGrid <- expand.grid(.mtry = mtry,
                             .splitrule =  "maxstat",
                             .min.node.size = c(20, 30))
  
  dir.create("tmp", showWarnings = F)
  if(nuitrient == "N"){
    trail_leng <- 200
    by <- 10
  }
  else if(nuitrient == "P"){
    trail_leng <- 60
    by <- 5
  }
  else if(nuitrient == "K"){
    trail_leng <- 60
    by <- 5
  }
  else (nuitrient == "S"){
    trail_leng <- 50
    by <- 5
  }
  
  vals <- c(seq(0, trail_leng, by = by))
  d <- data.frame(matrix(ncol = length(vals), nrow = nrow(covs_yield)))
  d_test <- data.frame(matrix(ncol = length(vals), nrow = nrow(covs_yield)))
  for(i in 1:length(vals)){
    names(d)[i] <- paste0(nuitrient,vals[i])
    d[,i] <- vals[i]
  }
  
  for(i in 1:ncol(d)){
    set.seed(7)
    rate <- d[,i]
    data_pred <- cbind(covs_yield,rate)
    inTrain <- createDataPartition(y =  data_pred$yield, p = 0.75, list = FALSE)
    training <- data_pred[inTrain,]
    testing <- data_pred[-inTrain,]
    
      mod_fit <- train(
        yield ~ .,
        data = training,
        method = "ranger",
        trControl = rf_fitControl,
        importance = 'impurity',
        tuneGrid = rf_tuneGrid,
        preProcess = c('scale', 'center')
      )

# ------------------------------------------------------------------------------  
# Plot and save variable importance
  var_imp <- varImp(mod_fit)
  ggplot(var_imp)
  
  ggsave(filename =
           paste0(
             "./predict_output/",
             paste("variable_importance", colnames(d)[i], sep = "_"),
             ".png"
           ))
  
# ------------------------------------------------------------------------------
#checking the stats and saving
  stat <- ggof(sim = pred_test, obs = testing$yield, ftype = "o", gofs = c("ME", "MAE", "RMSE", "NRMSE", "PBIAS", "R2"))
  ggplot(stat)
  ggsave(filename = paste0("./predict_output/", paste("stats", colnames(d)[i], sep = "_"), ".png"))

# ------------------------------------------------------------------------------
# predict using spatial grid data frame
  data <- data.frame(matrix(c(d[1,i]), ncol = 1, nrow = nrow(sgdf_final)))
  colnames(data) <- paste(nuitrent, "rate", sep = "_")
  data <- cbind(sgdf_final,data)
  pred_grid <- predict(mod_fit, data)

# ------------------------------------------------------------------------------
# converting to spatial grid data frame Exporting the prediction as raster
  pred_grid_sp <-
    data.frame(x = data$s1, y = data$s2, pred_grid)
  gridded(pred_grid_sp) <- ~x+y
  proj4string(pred_grid_sp) <- new_crs
 
  writeGDAL(
    pred_grid_sp,
    fname = paste0("./tmp/", paste(nuitrent,"predicted", i, sep = "_"),".tif"),
    drivername = "Gtiff",
    type = "Float32",
    options = "COMPRESS=DEFLATE"
  )
  }

# ------------------------------------------------------------------------------
# generating optimal yield
  setwd("tmp")
  yield <- list.files(path = ".", pattern = ".tif$", all.files = T)
  yield_ras <- lapply(yield, raster)
  yield_stack <- stack(yield_ras)
  optimal_yield <- max(yield_stack)
  writeRaster(
    optimal_yield,
    filename = paste(nuitrent, "optimal_yield", sep = "_"),
    format = "GTiff",
    overwrite = TRUE
  )
}

  