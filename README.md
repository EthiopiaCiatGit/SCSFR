## 1. About the Repo

    This repository consists of R Scripts that can be used to downloaded different covariates 
    that will be used for predicting fertilizer requirements and crop yields for an area of interest.
      
## 2. Scripts
      
#### 2.1 FetchSoilGRIDS.R
        
            used to access and download ISRIC Soil grids of 250m resolution usning a Web Coverage Service (WCS).
**arguments** 
              
              a) voi - the soil properties to be download
               
              b) depth - six standard depth intervals which soil predictions are made
              
              c) quantile - prediction quantiles (5% quantile, median of the distribution, mean of the distribution and 95% quantile)
              
              
#### 2.2 converSoilGRIDs.R
        
            - used to convert download ISRIC Soil grids to their conventional units using a conversion factor

**arguments**
              
              a) sou_path - the source path where the downloaded soil grids are located
               
              b) dest_path - the destination path for writing the output after conversion
              
##### 2.3 fetchTerraClim.R
        
            -  used to download meteorological TerraClimate dataset.

**arguments**
              
              a) param - the climatic parameter to be downloaded
              
              b) aoi - an area of interest mostly a shapefile
              
              c) csv - a csv that contains the date intervals of year, planting date and harvesting date
              
              d) aoi_path - the path of the area of interest
              
              e) csv_path - the path of the csv file
              
#### 2.4 fetchElevetion_slope_aspect_tpi.R
        
        used to download elevation data and calculate slope, aspect, topographic position index, 
        topographic ruggedness index and landform.

**arguments** 
              
              a) aoi - an area of interest mostly a shapefile
              
              b) aoi_path - the path of the area of interest
              
## 3. Steps to be taken

    fetchElevetion_slope_aspect_tpi.R and fetchTerraClim.R can be executed arbitrarily. 
    However FetchSoilGRIDS.R & converSoilGRIDS.R should be executed in their order before
    the prevoius scripts. 

## 4. Requirements

    - All the required packages should be installed before using the scripts. 
        In addition GDAL (Geospatial Data Abstraction Library) should be installed; 
        wcs_version ="VERSION=2.0.1" works for gdal >=2.3.
    - Install AIO and ClimateR using below command 
            #remotes::install_github("mikejohnson51/AOI") 
            #remotes::install_github("mikejohnson51/climateR")

## 5. Examples
### i. fetchSoilGRID.R
    
        fetchSoilGrid("clay", "0-5cm", "mean")
### ii. converSoilGRIDs.R
    
        a) converSoilGRIDs("D:/soil", "D:/final_output")
        b) converSoilGRIDs(sou_path = "D:/soil" , dest_path = "D:/final_output")
        c) converSoilGRIDs(sou_path = "D:/soil")
        
### iii. fetchTerraClim.R
        
        a) fetchTerraClimate(param = "srad", aoi ="eth.shp", aoi_path = "D:/test", csv ="test_data.csv", csv_path = "D:/test_data")
        
#### iv. fetchElevetion_slope_aspect_tpi.R
        
        a) fetchElevation(aoi = "tigray_prj", "D:/test")
        b) fetchElevation(aoi = "bg_prj")

## 6. Useful Links

1. [https://www.isric.org/explore/soilgrids/faq-soilgrids](url)
2. [https://www.climatologylab.org/terraclimate.html](url)
3. [https://github.com/mikejohnson51/climateR](url)
