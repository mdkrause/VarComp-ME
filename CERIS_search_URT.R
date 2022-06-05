########################################################################################################################
##                                                                                                                    ##
##  CERIS search                                                                                                      ##
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

# Loading CERIS function
source("CERIS_func.R") # does not change

# Loading packages
library(SoyURT)
library(dplyr)
library(data.table)
library(stringr)
library(asreml)

# data preparation 
data(pheno)
pheno$location <- as.factor(pheno$location)
pheno$year <- as.factor(pheno$year)
pheno$G <- as.factor(pheno$G)
pheno$average_planting_date <- as.Date(pheno$average_planting_date,format= "%m/%d/%y")
pheno$weights <- 1/pheno$SE^2 # weights as 1/SE^2

# obtaining genotype by location deviations within years
res <- c()
pheno <- pheno %>% group_by(year, location)
for(i in levels(pheno$year)){
model <- asreml(fixed = eBLUE ~ G + location, 
                family = asr_gaussian(dispersion = 1),
                weights = weights, 
                trace = FALSE,
                data = pheno %>% filter(year == i))
res <- rbind(res,model$residuals)
print(i)
}
pheno$GL <- as.vector(res^2)

# weather data
data(weather)
pheno$env <- as.factor(paste0(pheno$location, pheno$year)) # index for pheno (location x year)
weather$env <- as.factor(paste0(weather$location, # index for weather
                                format(as.Date(weather$YYYYMMDD,format= "%Y-%m-%d"), "%Y")))

# getting latitude/longitude info from pheno
lat <- pheno %>% group_by(env) %>% summarise(lat = unique(latitude))
long <- pheno %>% group_by(env) %>% summarise(long = unique(longitude))

# getting environment ID
env.i = levels(pheno$env) 

# computing planting / aproximate harvesting date and creating the metadata file "env_meta_info_0"
plant.date <- pheno %>% group_by(env) %>% summarise(pd = unique(average_planting_date)) %>% na.omit()
plant.date <- plant.date %>% group_by(env) %>% summarise(pd = mean.Date(pd,format=c("%Y-%m-%d"))) 
days_to_mature <- pheno %>% group_by(env) %>% summarise(dtm = unique(average_maturity_date)) %>% na.omit()
days_to_mature <- days_to_mature %>% group_by(env) %>% summarise(dtm = mean(dtm)) 
days_to_mature$dtm = round(days_to_mature$dtm,0)
plant.date <- left_join(plant.date, days_to_mature)
plant.date <- left_join(plant.date, lat)
env_meta_info_0 <- left_join(plant.date, long)
harv.date = as.Date(plant.date$pd) + (days_to_mature$dtm)
harv.date = as.character(harv.date)
plant.date = as.character(plant.date)

# extracting observed environments from downloaded weather data
weather <- weather %>% filter(env %in% unique(env_meta_info_0$env)) 

# extracting only necessary columns from phenotypic file
dat <- data.frame(pheno[,c('env','G','GL')])

# removnig the environment Booone - IA in 2016. For some reason, we were not able to
# download weather data for this one
dat <- dat %>% filter(env != 'boone_ia2016')

# The CERIS code requires character and more formatting
dat$G <- as.character(dat$G)
dat$env <- as.character(dat$env)
env_meta_info_0$pd <- env_meta_info_0$pd
env_meta_info_0$env <- as.character(env_meta_info_0$env)
weather$env <- as.character(weather$env)
weather$YYYYMMDD <- as.integer(format(as.Date(weather$YYYYMMDD), "%Y%m%d"))

# computing the mean GL value that will represent an environment
env_mean_trait <- dat %>% group_by(env) %>% summarise(meanY = mean(GL)) %>% as.data.frame()

# number of days after planting for CERIS to search, limited by the smallest one
# when computed across all environments
searching_daps <- min(days_to_mature$dtm)

# Parameters names
Paras <- colnames(weather)[-c(1:6, 26)]

## CERIS ##
# The Exhaustive_search() function must be loaded first! #

# parallel computing 
library(foreach)
library(doParallel)
cps <- detectCores() - 1
cl <- parallel::makeCluster(cps)
registerDoParallel(cl)

search.results <- Exhaustive_search(env_mean_trait, 
                                    weather,
                                    searching_daps, 
                                    dat$GL, 
                                    searching_daps,
                                    searching_daps, 
                                    LOO = 0,
                                    Paras, 
                                    window = 7)
stopCluster(cl)
head(search.results)
