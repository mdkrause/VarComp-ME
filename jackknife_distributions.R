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

# loading packages
library(SoyURT)
library(asreml)
library(foreach)
library(doParallel)
library(stringr)
library(dplyr)
library(ForestFit)

# phenotypic data
str(pheno)
pheno$location <- as.factor(pheno$location)
pheno$year <- as.factor(pheno$year)
pheno$G <- as.factor(pheno$G)
pheno$weights <- 1/(pheno$SE)^2

# example for the first period, from 1989 to 1995
dat <- filter(pheno, between(year, 1989, 1995)) %>% droplevels()
dat$env <- paste0(dat$year, dat$location)
Env <- unique(dat$env) # 181

##### jakknife in parallel #####
cps <- detectCores() - 1
cl <- parallel::makeCluster(cps)
registerDoParallel(cl)
        
results <- foreach(l=1:length(Env), .packages = c('asreml','dplyr', 'stringr'))  %dopar% {
  
  model <- asreml(fixed = Value ~ 1,
                  random = ~ location + G + year +
                    location:year + location:G + 
                    year:G + location:germplasmId:year,
                  family = asr_gaussian(dispersion = 1), 
                  weights = weights,
                  workspace = '1gb',
                  trace = FALSE, 
                  data = dat %>% filter(!env %in% Env[l]) %>% droplevels())
  
  res <- summary(model)$varcomp
  return(res)
}
stopCluster(cl)
           
##### probability distributions #####
## assuming there is a .rds files saved for each group of years

# extracting info from variances
permutation1 <- readRDS('years_1989_1995.rds') 
permutation2 <- readRDS('years_1996_2003.rds')
permutation3 <- readRDS('years_2004_2011.rds')
permutation4 <- readRDS('years_2012_2019.rds')
    
vG1 <- c() ; vG2 <- c(); vG3 <- c(); vG4 <- c()
vGL1 <- c() ; vGL2 <- c(); vGL3 <- c(); vGL4 <- c()
vGY1 <- c() ; vGY2 <- c(); vGY3 <- c(); vGY4 <- c()
vGLY1 <- c();vGLY2 <- c();vGLY3 <- c();vGLY4 <- c()

for(i in 1:length(permutation1)){
  
  vG1[i] = permutation1[[i]]['GermplasmId','component']
  vGL1[i] = permutation1[[i]]['Location:GermplasmId','component']
  vGY1[i] = permutation1[[i]]['year:GermplasmId','component']
  vGLY1[i] = permutation1[[i]]['Location:GermplasmId:year','component']
}

for(i in 1:length(permutation2)){
  
  vG2[i] = permutation2[[i]]['GermplasmId','component']
  vGL2[i] = permutation2[[i]]['Location:GermplasmId','component']
  vGY2[i] = permutation2[[i]]['year:GermplasmId','component']
  vGLY2[i] = permutation2[[i]]['Location:GermplasmId:year','component']
}

for(i in 1:length(permutation3)){
  vG3[i] = permutation3[[i]]['GermplasmId','component']
  vGL3[i] = permutation3[[i]]['Location:GermplasmId','component']
  vGY3[i] = permutation3[[i]]['year:GermplasmId','component']
  vGLY3[i] = permutation3[[i]]['Location:GermplasmId:year','component']
  
}

for(i in 1:length(permutation4)){
  vG4[i] = permutation4[[i]]['GermplasmId','component']
  vGL4[i] = permutation4[[i]]['Location:GermplasmId','component']
  vGY4[i] = permutation4[[i]]['year:GermplasmId','component']
  vGLY4[i] = permutation4[[i]]['Location:GermplasmId:year','component']
  
}

# getting together genotypic variance (can be any other)
values <- c(vG1, vG2, vG3, vG4)

# estimating probability distributions with ForestFit 

burr2 <- fitmixture(values, 'burr', 2)
f2 <- fitmixture(values, 'f', 2)
gamma2 <- fitmixture(values, 'gamma', 2)
llogis2 <- fitmixture(values, 'log-logistic', 2)
lnormal2 <- fitmixture(values, 'log-normal', 2)

burr3 <- fitmixture(values, 'burr', 3)
f3 <- fitmixture(values, 'f', 3)
gamma3 <- fitmixture(values, 'gamma', 3)
llogis3 <- fitmixture(values, 'log-logistic', 3)
lnormal3 <- fitmixture(values, 'log-normal', 3)

burr4 <- fitmixture(values, 'burr', 4)
f4 <- fitmixture(values, 'f', 4)
gamma4 <- fitmixture(values, 'gamma', 4)
llogis4 <- fitmixture(values, 'log-logistic', 4)
lnormal4 <- fitmixture(values, 'log-normal', 4)

burr5 <- fitmixture(values, 'burr', 5)
f5 <- fitmixture(values, 'f', 5)
gamma5 <- fitmixture(values, 'gamma', 5)
llogis5 <- fitmixture(values, 'log-logistic', 5)
lnormal5 <- fitmixture(values, 'log-normal', 5)

burr6 <- fitmixture(values, 'burr', 6)
f6 <- fitmixture(values, 'f', 6)
gamma6 <- fitmixture(values, 'gamma', 6)
llogis6 <- fitmixture(values, 'log-logistic', 6)
lnormal6 <- fitmixture(values, 'log-normal', 6)

# saving statistics 
stats <- rbind(gamma2$measures, llogis2$measures, lnormal2$measures, burr2$measures, f2$measures,
               gamma3$measures, llogis3$measures, lnormal3$measures, burr3$measures, f3$measures,
               gamma4$measures, llogis4$measures, lnormal4$measures, burr4$measures, f4$measures,
               gamma5$measures, llogis5$measures, lnormal5$measures, burr5$measures, f5$measures,
               gamma6$measures, llogis6$measures, lnormal6$measures, burr6$measures, f6$measures)

rownames(stats) <- c('gamma2', 'llogis2', 'lnormal2', 'burr2', 'F2',
                     'gamma3', 'llogis3', 'lnormal3', 'burr3', 'F3',
                     'gamma4', 'llogis4', 'lnormal4', 'burr4', 'F4',
                     'gamma5', 'llogis5', 'lnormal5', 'burr5', 'F5',
                     'gamma6', 'llogis6', 'lnormal6', 'burr6', 'F6')

capture.output(stats, file = "statsG_LOO.txt")
