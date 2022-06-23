########################################################################################################################
##                                                                                                                    ##
##  Download and process soil data                                                                                    ##
##  Date: Jun, 2022                                                                                                   ##
##                                                                                                                    ##
##  Using large soybean historical data to study genotype by environment variation and identify mega-environments     ##
##  with the integration of genetic and non-genetic factors                                                           ##
##  Krause et al. (2022)                                                                                              ##
##                                                                                                                    ##
##  Authors:    MD Krause        <krause.d.matheus@gmail.com>                                                         ##
##              WD Beavis        <wdbeavis@iastate.edu>                                                               ##
##                                                                                                                    ##
########################################################################################################################

########################################################################################################################
##                                                  DISCLAIMER                                                        ##
##      This code was forked from https://github.com/zecojls/downloadSoilGridsV2¨ (Dr. José Lucas Safanelli) in 2020  ##
##  We basically inplemented a grid search approach to make sure we would obtain information for all locations, added ##
##  some tricks, and implemented parallel computing                                                                   ##  
##                                                                                                                    ##
##  The original one can be founded at:                                                                               ##   
##  https://github.com/zecojls/downloadSoilGridsV2/blob/master/script_soilGrids_download.R                            ## 
########################################################################################################################

# Setup
options(stringsAsFactors = FALSE)
options(pillar.sigfig=3)

# loading packages
library(foreach)
library(doParallel)
library(dplyr)

## Directories
dir.proj <- getwd() # must create a folder cal
dir.create("raster")
# saving temporary files, that will be automatically removed
dir.export <- paste0("./raster/") 

# phenotypic data to obtain latitude, longitude, and locations names
library(SoyURT)
data(pheno)
geo <- pheno %>% group_by(location) %>% summarise(latitude, longitude)
geo <- distinct(geo)
colnames(geo)<- c('Location', 'lat', 'long')

################################################################################
# If you are interested in downloading soil data for another locations, just 
# replace the step below for your own location name, latitude, and
# longitude, and voilà! :)
#
# For example:
# geo <- data.frame(Location = c("Pato Branco", "Lyon", "Accra"),
#                   lat = c(-26.2295, 45.767592, 5.614818),
#                   long = c(-52.6716, 4.676925, -0.205874))
#
################################################################################


# Obtaining the soil data:
cps <- detectCores() - 1
cl <- parallel::makeCluster(cps)
registerDoParallel(cl)

infoSoil <- foreach(ENV=1:nrow(geo), .errorhandling='pass', .combine = 'rbind',
                    .packages = c('curl','XML','dplyr', 'EnvRtype', 'raster', 'stringr'))%dopar% {

  min.long <- min(geo$long[ENV])
  min.lat <- min(geo$lat[ENV])
  max.long <- max(geo$long[ENV])
  max.lat <- max(geo$lat[ENV])
  
  seq.long <- seq(min.long-0.01, max.long+0.01, 0.01)
  seq.lat <- seq(min.lat-0.01, max.lat+0.01, 0.01)
  
  combination.min <- expand.grid(seq.long[-length(seq.long)], seq.lat[-length(seq.lat)])
  combination.max <- expand.grid(seq.long[-1], seq.lat[-1])
  
  full.combination <- tibble(min.long = combination.min[,1],
                             max.long = combination.max[,1],
                             min.lat = combination.min[,2],
                             max.lat = combination.max[,2])
  
  full.combination <- full.combination %>%
    mutate(min.long = min.long - 0.01,
           max.long = max.long + 0.01,
           min.lat = min.lat - 0.01,
           max.lat = max.lat + 0.01)
  
  full.combination <- as.data.frame(full.combination)
  
  bbox.coordinates <- full.combination %>%
    mutate(left.coord = paste0(ifelse(min.long < 0, "W", "E"), round(abs(min.long), 0)),
           top.coord = paste0(ifelse(max.lat < 0, "S", "N"), round(abs(max.lat), 0)))
  
  # Download links
  
  # WRB
  #"https://maps.isric.org/mapserv?map=/map/wrb.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=MostProbable&FORMAT=image/tiff&SUBSET=long(-54.2280,-52.2280)&SUBSET=lat(-22.0906,-20.0906)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326"
  # pH
  #'https://maps.isric.org/mapserv?map=/map/phh2o.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=phh2o_0-5cm_mean&FORMAT=image/tiff&SUBSET=long(-51.8169,-49.8169)&SUBSET=lat(-20.9119,-18.9119)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326'
  # SOC
  #"https://maps.isric.org/mapserv?map=/map/soc.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=soc_0-5cm_mean&FORMAT=image/tiff&SUBSET=long(-52.0848,-50.0848)&SUBSET=lat(-17.2684,-15.2684)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326"
  # N
  #"https://maps.isric.org/mapserv?map=/map/nitrogen.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=nitrogen_0-5cm_mean&FORMAT=image/tiff&SUBSET=long(-49.0307,-47.0307)&SUBSET=lat(-20.4832,-18.4832)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326"
  # CTC
  #"https://maps.isric.org/mapserv?map=/map/cec.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=cec_0-5cm_mean&FORMAT=image/tiff&SUBSET=long(-49.2986,-47.2986)&SUBSET=lat(-23.7516,-21.7516)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326"
  # Silt
  #"https://maps.isric.org/mapserv?map=/map/silt.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=silt_0-5cm_mean&FORMAT=image/tiff&SUBSET=long(-51.1739,-49.1739)&SUBSET=lat(-20.1082,-18.1082)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326"
  # Clay
  #"https://maps.isric.org/mapserv?map=/map/clay.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=clay_0-5cm_mean&FORMAT=image/tiff&SUBSET=long(-52.8950,-50.8950)&SUBSET=lat(-19.4116,-17.4116)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326"
  # Sand
  #"https://maps.isric.org/mapserv?map=/map/sand.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=sand_0-5cm_mean&FORMAT=image/tiff&SUBSET=long(-48.3342,-46.3342)&SUBSET=lat(-19.1437,-17.1437)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326"
  # BD
  #'https://maps.isric.org/mapserv?map=/map/bdod.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=bdod_0-5cm_mean&FORMAT=image/tiff&SUBSET=long(-50.9661,-48.9661)&SUBSET=lat(-18.0721,-16.0721)&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326'
  
  # Automatic download
  
  attributes <- c("wrb.map", "phh2o.map", "soc.map", "nitrogen.map",
                  "cec.map", "silt.map", "clay.map", "sand.map", "bdod.map")
  
  #layers <- c("0-5cm_mean", "5-15cm_mean", "15-30cm_mean", "30-60cm_mean")
  layers <- c("5-15cm_mean")
  
  for(a in 1:length(attributes)) {
    
    attribute <- attributes[a]
    
    attribute.prefix <- gsub(".map", "", attribute)
    
    if(attribute == "wrb.map") {
      
      layer <- "MostProbable"
      
        for(t in 1:nrow(bbox.coordinates)) {
        
        min.long = bbox.coordinates[t,"min.long"]
        max.long = bbox.coordinates[t,"max.long"]
        min.lat = bbox.coordinates[t,"min.lat"]
        max.lat = bbox.coordinates[t,"max.lat"]
        left.coord <- bbox.coordinates[t,"left.coord"]
        top.coord <- bbox.coordinates[t,"top.coord"]
        
        wcs <- paste0("https://maps.isric.org/mapserv?map=/map/", attribute, "&",
                      "SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=", layer, "&",
                      "FORMAT=image/tiff&",
                      "SUBSET=long(", min.long, ",", max.long, ")&",
                      "SUBSET=lat(", min.lat, ",", max.lat, ")&",
                      "SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326")
        
        destination.file <- paste0(dir.export, "/SoilGrids_",
                                   paste(attribute.prefix, layer,
                                         left.coord, top.coord, "env", ENV,sep = "_"),
                                   ".tif")
        
        if(file.exists(destination.file)) {
          
          next
          
        } else {
          
          cat("Downloading: ", destination.file, "\n")
          download.file(wcs, destfile = destination.file, mode = 'wb')
          
        }
        
      }
      
    } else {
      
      for(l in 1:length(layers)) {
        
        layer <- layers[l]
        
          for(t in 1:nrow(bbox.coordinates)) {
          
          min.long = bbox.coordinates[t, "min.long"]
          max.long = bbox.coordinates[t, "max.long"]
          min.lat = bbox.coordinates[t, "min.lat"]
          max.lat = bbox.coordinates[t, "max.lat"]
          left.coord <- bbox.coordinates[t, "left.coord"]
          top.coord <- bbox.coordinates[t, "top.coord"]
          
          wcs <- paste0("https://maps.isric.org/mapserv?map=/map/", attribute, "&",
                        "SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=", attribute.prefix, "_", layer, "&",
                        "FORMAT=image/tiff&",
                        "SUBSET=long(", min.long, ",", max.long, ")&",
                        "SUBSET=lat(", min.lat, ",", max.lat, ")&",
                        "SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326")
          
          destination.file <- paste0(dir.export, "/SoilGrids_",
                                     paste(attribute.prefix, layer,
                                           left.coord, top.coord, "env", ENV,sep = "_"),
                                     ".tif")
          
          if(file.exists(destination.file)) {
            
            next
            
          } else {
            
            cat("Downloading: ", destination.file, "\n")
            download.file(wcs, destfile = destination.file, mode = 'wb')
            
          }
        }
      }
    }
  }
  

  # after downloading the .tif files, let's process them!
  dir = dir.export
  soil_grid = list.files(path = dir,pattern = paste0("env_", ENV, ".tif"))
  soil_name = gsub(soil_grid,pattern = paste0("_env_", ENV, ".tif"),replacement = '')
  
  env.data = data.frame(env = geo$Location[ENV],LAT = geo$lat[ENV], LON = geo$long[ENV])
  
  soil_data = c()
  for(i in 1:length(soil_grid)){
    soil_data = rbind(soil_data,data.frame(Feature = soil_name[i],
                                           extract_GIS(covraster = raster(paste0(dir,'/',soil_grid[i])),
                                                       name.out = 'Soil_Grid',env.data = env.data)))
  }

  ### grid search by walking through lat and long by 0.00001
  if(any(soil_data$Soil_Grid == 0 | is.na(soil_data$Soil_Grid))){
    env.data = data.frame(env = geo$Location[ENV],
                          LAT = min(bbox.coordinates$min.lat),
                          LON = min(bbox.coordinates$min.long))
  add = 0
  while(any(soil_data$Soil_Grid == 0 | is.na(soil_data$Soil_Grid))){
    env.data$LAT <- env.data$LAT + add
    env.data$LON <- env.data$LON + add
    soil_data = c()
    for(i in 1:length(soil_grid)){
      soil_data = rbind(soil_data,data.frame(Feature = soil_name[i],
                                             extract_GIS(covraster = raster(paste0(dir,'/',soil_grid[i])),
                                                         name.out = 'Soil_Grid', env.data = env.data)))
    }
    add = add + 0.00001
        if(env.data$LAT > max(bbox.coordinates$max.lat) | env.data$LON > max(bbox.coordinates$max.long)){
      cat('Please replace the coordinates for', env.data$env)
      break}
  }
}
  
  # some formatting
  soil_data$Feature[str_detect(soil_data$Feature, "bdod_5-15cm_mean")] <- "bdod_5-15cm_mean"
  soil_data$Feature[str_detect(soil_data$Feature, "cec_5-15cm_mean")] <- "cec_5-15cm_mean"
  soil_data$Feature[str_detect(soil_data$Feature, "clay_5-15cm_mean")] <- "clay_5-15cm_mean"
  soil_data$Feature[str_detect(soil_data$Feature, "nitrogen_5-15cm_mean")] <- "nitrogen_5-15cm_mean"
  soil_data$Feature[str_detect(soil_data$Feature, "phh2o_5-15cm_mean")] <- "phh2o_5-15cm_mean"
  soil_data$Feature[str_detect(soil_data$Feature, "sand_5-15cm_mean")] <- "sand_5-15cm_mean"
  soil_data$Feature[str_detect(soil_data$Feature, "silt_5-15cm_mean")] <- "silt_5-15cm_mean"
  soil_data$Feature[str_detect(soil_data$Feature, "soc_5-15cm_mean")] <- "soc_5-15cm_mean"
  soil_data$Feature[str_detect(soil_data$Feature, "wrb_MostProbabl")] <- "wrb_MostProbabl"
  
  unlink(x = paste0(dir.export,"/",soil_grid))
  
  soil_data <- distinct(soil_data) 
  #write.csv(soil_data, file = paste0("./raw_data_all/env",ENV,".csv"))
  return(soil_data)
}
stopCluster(cl)

# the data:

## If there is missing data for a given location, please slightly modify its 
## geographic coordinates and run the code again

head(infoSoil)
