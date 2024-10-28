library(tidyverse) 
library(caret) 
library(doParallel) 
library(hydroGOF)
library(randomForest)
library(h2o)

rm(list = ls())

h2o.init()
setwd("~/01My_Docs/01CIAT/02Fertilizer/prediction_2024/maize/data/final/")
data <- read.csv("analysis_ready_maize_soiltype.csv", header = T, sep = ",")
dim(data)
colnames(data)
train_data <- data[, -c(1:8)] %>% na.omit() %>% unique()
dim(train_data)
glimpse(train_data)
colnames(train_data)

#change factor variables
train_data$soil_type <- as.factor(train_data$soil_type)

set.seed(2021)
rfe_maize <- readRDS("~/01My_Docs/01CIAT/02Fertilizer/prediction_2024/maize/data/workspace/rfe_maize.rds")

#training
list_of_var <- rfe_maize$optVariables[1:25]
list_of_var <- append(c("grain_yield_kgpha"), list_of_var)
plot(rfe_maize)
# select most important variables only

#hyper parameter tuning 


data <- data[,list_of_var]
data$soil_type <- as.factor(data$soil_type)
colnames(data)
glimpse(data)

ML_inputData.h2o <- as.h2o(data)
ML_inputData_split <- h2o.splitFrame(data = ML_inputData.h2o, ratios = 0.7, seed = 1234)
training_data <- ML_inputData_split[[1]]
test_data <- ML_inputData_split[[2]]

response <- "grain_yield_kgpha"
#predictors <- ML_inputData %>% dplyr::select(-c(TLID, P_base_supply, K_base_supply, lon, lat)) %>% names()
predictors <- data |> select(-grain_yield_kgpha) |> names()


hyper_params <- list(
  ntrees = seq(20, 200, 10),
  max_depth = seq(4, 12, 2),
  mtries = c(2, 3, 4, 5, 6)
)

# Train and tune the random forest model
grid <- h2o.grid(
  algorithm = "randomForest",
  x = predictors,
  y = response,
  grid_id = "rf_grid",
  hyper_params = hyper_params,
  training_frame = training_data,
  validation_frame = test_data,
  seed = 444
)

# Get the best model from the grid search
best_model <- h2o.getModel(grid@model_ids[[1]])

# View the hyperparameters of the best model
print(best_model@parameters) ## mtry = 6, ntrees = 80, max_Depth = 12


ML_randomForest <- h2o.randomForest(x = predictors,
                                    y = response,
                                    ntrees =  best_model@parameters$ntrees,
                                    max_depth = best_model@parameters$max_depth,
                                    mtries =   best_model@parameters$mtries,
                                    training_frame = training_data,
                                    validation_frame = test_data,
                                    keep_cross_validation_predictions = TRUE,
                                    nfolds = 5,
                                    seed = 444)


rmse_r2_randomforest <- data.frame( rmse =h2o.rmse(ML_randomForest, train=TRUE, valid=TRUE),
                                    R_sq = c(h2o.r2(ML_randomForest, train=TRUE, valid=TRUE)))
rmse_r2_randomforest


rf_valid <- test_data
rf_valid$predResponse <- h2o.predict(object = ML_randomForest, newdata = test_data)
rf_valid <- as.data.frame(rf_valid)
rf_valid$Response <- rf_valid[,which(names(rf_valid)==response)]

ggplot(rf_valid, aes(Response, predResponse)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, col = "blue")+
  xlab("Measured yield") + ylab("predicted yield")+
  ggtitle("Random forest") +
  #xlim(0,max(ML_inputData$Yield)) + ylim(0,max(ML_inputData$Yield))+ 
  theme_minimal()


#the variable importance plot
par(mar = c(1, 1, 1, 1))#Expand the plot layout pane
h2o.varimp_plot(ML_randomForest)


# shap values = the direction of the relationship between our features and target
# e.g., high vlaues of total rainfall has positive contribution
h2o.shap_summary_plot(ML_randomForest, test_data, num_of_features = 25)
