
library(tidyverse)
setwd("C:\\Users\\ATilaye\\Documents\\01My_Docs\\01CIAT\\02Fertilizer\\prediction_2024\\maize\\data\\input")
list.files()

rm(list = ls())
#read maize data
maize <- read.csv("maize_dst_team_shared_v2.csv", header = T, sep = ",")
dim(maize)

#remove SG-2000 data
maize_less_sg2000 <- maize |> dplyr::filter(source != "SG 2000")
dim(maize_less_sg2000)
view(maize_less_sg2000)

#remove yield more than 14ton/ha
maize_less_14ton <- maize_less_sg2000 |> filter(grain_yield_kgpha <= 14000)
dim(maize_less_14ton)

#write the final filtered csv
setwd("../final/")
write.csv(maize_less_14ton, "final_filtered_maize.csv", col.names = T, row.names = F)
summary(maize_less_14ton)
