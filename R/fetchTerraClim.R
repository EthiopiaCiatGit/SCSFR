# ------------------------------------------------------------------------------
# Fetches climate data from TerraClim
# ------------------------------------------------------------------------------ 

fetchTerraClimate <-
  function(param,
           aoi,csv, 
           aoi_path,
           csv_path) {
    
    library(AOI)
    library(rgdal)
    library(climateR)
    library(sf)
    library(raster)
    library(dplyr)
    
    file_name <- as.character(param)
    aoi <- st_read(paste(dsn = aoi_path, layer = aoi, sep = "/"))
    csv <- read.csv(file = paste(csv_path, csv, sep = "/"), header = TRUE, sep = ",")
    # }
    
data <- 
  dplyr::select (csv,
                 "Region",
                 "District",
                 "Site",
                 "Year",
                 "Planting_date",
                 "Harvesting_date")
data$Year <- as.character(data$Year) 
data$Year <- paste(data$Year, "01", "01", sep = "-")
data$Year <- as.Date(data$Year)
data2 <-
  data %>% mutate(
    plating_date2 = Year + Planting_date - 10,
    harvesting_date2 = Year + Harvesting_date + 10
  )
data3 <- 
  dplyr::select (data2,
                 "Year",
                 "plating_date2",
                 "harvesting_date2")
data3 <- unique.data.frame(data3)
start_date <- data3$plating_date2
end_date <- data3$harvesting_date2
output <- list()

for(i in 1:length(start_date)){
  pout = getTerraClim(AOI = aoi,
                      param = param,
                      startDate = start_date[i],
                      endDate = end_date[i])
  pout <- pout[[1]]
  output  <- append(output, pout)
  for(j in 1:length(output)){
    pout <- output[[j]]
    rname <- c(names(pout))
    for (k in 1:length(rname)) {
      writeRaster(
        pout,
        filename = paste(file_name, start_date[i], end_date[i], sep = "_"),
        bylayer = T,
        format = "GTiff",
        overwrite = TRUE
      )
    }
  }
}
  }


