#-------------------------------------------------------------------------------   
# Function model training, testing and prediction
# ------------------------------------------------------------------------------

#' @param covs_csv_path - path for the extracted covariates and trial data
#' @param covs_csv - covariates data
#' @param stack_path - path for the raster of different scenarios
#' @param stack - raster stacks for different scenario
#' @param dest_path - path for the output

model_caliPredict <- function(covs_csv_path, covs_csv, stack_path, stack, 
                              dest_path){
  
  
  packages_required <- c("terra","dplyr", "sf", "caret", "HydroGOF","ggplot2", 
                         "Metrics", "randomforest", "rgdal", "ranger")
  
  # check and install packages that are not yet installed
  installed_packages <- packages_required %in% rownames(installed.packages())
  if(any(installed_packages == FALSE)){
    install.packages(packages_required[!installed_packages])}
  
  # load required packages
  invisible(lapply(packages_required, library, character.only = TRUE))
  
  # read Covariates with trial data and filter the data
  setwd(covs_csv_path)
  covs <- read.csv(covs_csv, header = T, sep = ",") %>%
          select(-c(id, long, lat, year, crop)) %>% na.omit() %>%
          unique() %>% filter(gr_yield > 0 & gr_yield <= 15000)
  
  #read the raster stack (for instance the stack normal)
  setwd(stack_path)
  normal <- rast(stack)
  names(stack) <-
    names(stack_normal) 
  c(
    "bdod",
    "cec",
    "cfvo",
    "clay",
    "dem",
    "landform",
    "nitrogen",
    "ocd",
    "ocs",
    "phh2o",
    "sand",
    "silt",
    "slope",
    "soc",
    "tpi",
    "tri", 
    "prcp1",
    "prcp2",
    "prcp3",       
    "srad1",
    "srad2",        
    "srad3",
    "tmax1",        
    "tmax2",       
    "tmax3",       
    "tmin1",      
    "tmin2",       
    "tmin3"
  )
  
  #training
  set.seed(123)
  mtry <- as.integer((ncol(covs))/3) #this will be optimized
  mtry <- seq(mtry - 8, mtry + 8, by = 2)
  
  rf_fitControl <- trainControl(method = "repeatedcv",
                                number = 10,
                                repeats = 5)
  
  rf_tuneGrid <- expand.grid(.mtry = mtry,
                             .splitrule =  "maxstat",
                             .min.node.size = c(20, 30))
  
  inTrain <- createDataPartition(y =  covs$gr_yield, p = 0.70, list = FALSE)
  training <- cov_train[inTrain,]
  testing <- cov_train[-inTrain,]
  
  message(noquote("Training the model..."))
  mod_fit <- train(
    gr_yield ~ .,
    data = training,
    method = "ranger",
    trControl = rf_fitControl,
    importance = 'impurity',
    tuneGrid = rf_tuneGrid,
    preProcess = c('scale', 'center'))
  # ------------------------------------------------------------------------------  
  # Plot and save variable importance and testing
  var_imp <- varImp(mod_fit)
  ggplot(var_imp)
  
  #Model validation to be included??
  
  # ------------------------------------------------------------------------------  
  # Multiple iteration
  N <- c(seq(0, 75, 15), seq(100, 200, 25)) #this will be optimized
  P <- K <- N[1:7]
  npk <- expand.grid(N = N, P = P, K = K)
  
  #save important parameters
  setwd(dest_path)
  save(cov_train, file = paste0("regression_matrix", ".RData"))
  save(mod_fit, file = paste0("model_with_soil", ".RData"), Overwrite = T)
  ggsave(filename = "variable_importance.png", width = 20, height = 10)
  
  path <- "yield_normal"
  dir.create(path, FALSE, TRUE)
  a <- apply(npk, 1, function(i) paste(i, collapse="."))
  f <- file.path(path, paste0("yield.", a, ".tif"))
  #message(noquote("Predicting yield"))
  #pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
  
  #loop for climate forcast scenario
  for (i in 1:nrow(npk)) {
    if (file.exists(f[i])) next
    NPK <- data.frame(n = npk$N[i], p = npk$P[i], k = npk$K[i])
    predict(
      stack_normal, # use below normal, average and above average
      mod_fit,
      const = NPK,
      filename = f[i],
      overwrite = TRUE,
      wopt = list(datatype = "INT2S", names = a[i])
    )
  }
}