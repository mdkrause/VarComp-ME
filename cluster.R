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

# necessary packages
library(dplyr)
library(factoextra)
library(ggsci)
library(pcaMethods)
library(cluster)
library(SoyURT)

# loading pheno and soil data
data(pheno)
data(soil)

# formatting from tidy to messy
soil <- soil %>% select(location, Feature, Soil_Grid) %>% 
  tidyr::pivot_wider(names_from = Feature, values_from = Soil_Grid) %>%
  as.data.frame()

# adding altitude to soil data
alt <- pheno %>% group_by(location) %>% summarise(alt = unique(altitude))
soil <- left_join(soil, alt)

# non-linear PCA
PCA <- prep(soil[,-1], scale = 'uv', center = TRUE)
PCA <- pca(PCA, method="nipals", nPcs=7) 

# selecting the first 5 PC's because they explained > 90% of variation
PCA
PCA <- scores(PCA)[, c(1:5)] 
rownames(PCA) <- soil[,1]

## Silhouette
fviz_nbclust(PCA, kmeans, method = "silhouette", k.max = 10)

## Elbow
fviz_nbclust(PCA, kmeans, method = 'wss', k.max = 10)

# Finally, clustering
set.seed(1)
cluster <- kmeans(PCA, centers = 2, iter.max = 1000, nstart = 100)

fviz_cluster(cluster, data = PCA, axes = c(1,2), 
             show.clust.cent = FALSE, shape = FALSE,
             ellipse.type = "convex") 

# getting indexes
cluster$cluster
