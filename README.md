# SCSFR (Site and Context Specific Fertilizer Recommendation)

The fertilizer recommendation focuses on predicting fertilizer nutrient and crop yield. The script can be used to automatically access grided covariates such as ISRIC soil data, TerraClim climate data, SRTM topographic data from amazon web services and apply random forest machine learning algorithm to predict the required fertilizer nutrient for optimal crop yield.
      
## 1. Scripts
      
#### 1.1 fetchSoilGRID.R
        
#### - Used to access and download ISRIC Soil grids of 250m resolution usning a Web Coverage Service (WCS).
### **arguments** 
              
    a) voi - the soil properties to be downloaded

    b) depth - six standard depth intervals which soil predictions are made

    c) quantile - prediction quantiles (5% quantile, median of the distribution, mean of the distribution and 95% quantile)

    d) aoi - an area of interest mostly a shapefile

    e) aoi_path - the path of the area of interest

              
#### 1.2 convertSoilGRID.R
        
#### - used to convert download ISRIC Soil grids to their conventional units using a conversion factor

**arguments**
              
    a) sou_path - the source path where the downloaded soil grids are located

    b) dest_path - the destination path for writing the output after conversion
              
##### 1.3 fetchTerraClim.R
        
#### -  Used to download meteorological TerraClimate dataset.

**arguments**
              
    a) param - the climatic parameter to be downloaded

    b) aoi - an area of interest mostly a shapefile

    c) csv - a csv that contains the date intervals of year, planting date and harvesting date

    d) aoi_path - the path of the area of interest

    e) csv_path - the path of the csv file
              
#### 1.4 fetchElevation.R
        
#### - Used to download elevation data and calculate slope, aspect, topographic position index, topographic ruggedness index and landform.

**arguments** 
              
    a) aoi - an area of interest mostly a shapefile

    b) aoi_path - the path of the area of interest
              
## 2. Steps to be taken

#### - fetchElevetion and fetchTerraClim.R can be executed arbitrarily. However FetchSoilGRIDS.R & converSoilGRIDS.R 
should be executed in their order before the prevoius scripts. 

## 3. Requirements

#### - All the required packages should be installed before using the scripts. In addition GDAL (Geospatial Data Abstraction Library)should be installed for dwonloading soil covariates. wcs_version ="VERSION=2.0.1" works for gdal >=2.3.
#### - Install AIO and ClimateR for downloading climate covariates using below commands
    remotes::install_github("mikejohnson51/AOI") 
    remotes::install_github("mikejohnson51/climateR")

## 4. Examples
#### i. fetchSoilGRID.R
    fetchSoilGrid("clay", "0-5cm", "mean")
#### ii. convertSoilGRID.R
    a) converSoilGRIDs("D:/Test", "D:/Test2")
    b) converSoilGRIDs(sou_path = "D:/Test" , dest_path = "D:/Test2")
    c) converSoilGRIDs(sou_path = "D:/Test")
#### iii. fetchTerraClim.R
        
    fetchTerraClimate(param = "srad", aoi ="eth.shp", aoi_path = "D:/Test", csv ="test_data.csv", csv_path = "D:/Test_csv")
        
#### iv. fetchElevetion.R

    a) fetchElevation(aoi = "eth", "D:/Test")
    b) fetchElevation(aoi = "eth")

## 5. Useful Links

1. [www.isric.org/explore/soilgrids/faq-soilgrids](https://www.isric.org/explore/soilgrids/faq-soilgrids)
2. [www.climatologylab.org/terraclimate.html](https://www.climatologylab.org/terraclimate.html)
3. [github.com/mikejohnson51/climateR](https://github.com/mikejohnson51/climateR)
