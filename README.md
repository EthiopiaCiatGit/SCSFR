**1. About the Repo**

    This repository consists of R Scripts that can be used to downloaded different covariates that will be used for predicting fertilizer requirements and crop yields for an area of interest.
      
      **Scripts**
      
        **1.1 FetchSoilGRIDS.R**
        
            - used to access and download ISRIC Soil grids of 250m resolution

             _ arguements _
              
              a) voi - the soil properties to be download
               
              b) depth - six standard depth intervals which soil predictions are made
              
              c) quantile - prediction quantiles (5% quantile, median of the distribution, mean of the distribution and 95% quantile)
              
              
       ** 1.2 converSoilGRIDS.R**
        
            - used to convert download ISRIC Soil grids to their conventional units using a conversion factor

             _ arguements _
              
              a) sou_path - the source path where the downloaded soil grids are located
               
              b) dest_path - the destination path for writing the output after conversion
              
          **1.3 fetchTerraClim.R**
        
            -  used to download meteorological TerraClimate dataset.

              _arguements_ 
              
              a) param - the climatic parameter to be downloaded
              
              b) aoi - an area of interest mostly a shapefile
              
              c) csv - a csv that contains the date intervals of year, planting date and harvesting date
              
              d) aoi_path - the path of the area of interest
              
              e) csv_path - the path of the csv file
              
          **1.3 fetchElevetion_slope_aspect_tpi.R**
        
            -  used to download elevation data and calculate slope, aspect, topographic position index, topographic ruggedness index and landform.

              _arguements_ 
              
              a) aoi  - an area of interest mostly a shapefile
              
              b) aoi_path - the path of the area of interest
              
**2. Steps to be taken**


**4. Requirements**


**6. Examples**


**8. Useful Links**
