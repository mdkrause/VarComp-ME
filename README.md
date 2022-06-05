# Using large soybean historical data to study genotype by environment variation and identify mega-environments with the integration of genetic and non-genetic factors
Krause, M.D., Dias, K.O.G., Singh, A.K., Beavis, W.D. Field Crop Research. 
https://doi.org/10.1007/s00122-022-04041-y 

The following files contain:

CERIS_func: Script for the Critical Environmental Regressor through Informed Search (CERIS).

CERIS_search_URT: Script to implement CERIS.

MET_models: Script for the METs models M3_1 to M3_20.

cluster: Script for the non-linear PCA and k-means clustering. This script implements the analysis performed for the soil data, which can be replicated for both weather and loadings from FA models.

soil_data: Script for downloading soil data from SoilGrids.

jackknife_distributions: Script for computing modified jackknife resampling technique and estimating density functions.

How to download weather data? We do have R script for it, but we do recommend users to drink from the well at https://github.com/allogamous/EnvRtype.

All used data sets (pheno + soil + weather) are available in the SoyURT package, which can be obtained on CRAN.
