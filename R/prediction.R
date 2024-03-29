# ------------------------------------------------------------------------------
# souring data preparation script

  #rm(list = ls())
  setwd(path_to_input_script)
  source("input.R")

# ------------------------------------------------------------------------------
# extracting values by points
  message(noquote("Extracting point vaues..."))
  cov_soil_land <- stack_normal[[c(1:17)]]#stack of soil and landform
  grid_val <- extract(cov_soil_land, points)
  cov_train <- cbind(grid_val, points_train) # covariates, yield & nps
  cov_train <- cov_train[complete.cases(cov_train),]
  #cov_train$landform <- as.factor(cov_train$landform)
  #cov_train$soil_type <- as.factor(cov_train$soil_type)
  yield_nps <- dplyr::select(cov_train, c("n","p","k", "yield"))
  cov_train <- dplyr::select(cov_train, -c("n","p","k", "yield"))
  cov_train <- cbind(cov_train, yield_nps)
  cov_train <- unique(na.omit(cov_train[, 1:ncol(cov_train)])) #removing NAs and duplicates
  cov_train <- cov_train[(which(cov_train$yield <= 15000)), ]
  
  setwd(path_to_output)
  save(cov_train, file = paste0("regression_matrix", ".RData"))
  
  # training
  set.seed(123)
  mtry <- as.integer((ncol(cov_train))/3) #this will be optimized
  mtry <- seq(mtry - 8, mtry + 8, by = 2)
  
  rf_fitControl <- trainControl(method = "repeatedcv",
                                number = 10,
                                repeats = 5)
  
  rf_tuneGrid <- expand.grid(.mtry = mtry,
                             .splitrule =  "maxstat",
                             .min.node.size = c(20, 30))
  
  inTrain <- createDataPartition(y =  cov_train$yield, p = 0.70, list = FALSE)
  training <- cov_train[inTrain,]
  testing <- cov_train[-inTrain,]
  
  message(noquote("Training the model..."))
  mod_fit <- train(
    yield ~ .,
    data = training,
    method = "ranger",
    trControl = rf_fitControl,
    importance = 'impurity',
    tuneGrid = rf_tuneGrid,
    preProcess = c('scale', 'center'))
  save(mod_fit, file = paste0("model_with_soil", ".RData"), Overwrite = T)
# ------------------------------------------------------------------------------  
# Plot and save variable importance and testing
  var_imp <- varImp(mod_fit)
  ggplot(var_imp)
  
  ggsave(filename = "variable_importance.png", width = 20, height = 10)
  
# ------------------------------------------------------------------------------  
# Multiple iteration
  N <- c(seq(0, 75, 15), seq(100, 200, 25)) #this will be optimized
  P <- K <- N[1:7]
  # P <- K <- N[1:3]
  npk <- expand.grid(N = N, P = P, K = K)
  
#loop for climate forcast scenario
  setwd(path_to_output)
  path <- "yield_normal"
  dir.create(path, FALSE, TRUE)
  a <- apply(npk, 1, function(i) paste(i, collapse="."))
  f <- file.path(path, paste0("yield.", a, ".tif"))
  setwd(path_to_output)
  message(noquote("Predicting yield"))
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
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

  
