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

# Loading packages
library(asreml)
asreml.options(pworkspace="3gb", workspace="3gb", 
               maxit=200, uspd=F, 
               rotate.fa = FALSE) # why FALSE? Because we will rotate it later on! :) 
library(dplyr) 
library(SoyURT)

# phenotypic data
str(pheno)
pheno$location <- as.factor(pheno$location)
pheno$year <- as.factor(pheno$year)
pheno$G <- as.factor(pheno$G)
pheno$weights <- 1/(pheno$SE)^2
L <- length(unique(pheno$location))
Y <- length(unique(pheno$year))

# Models

M3.1 <- asreml(fixed = Value ~ location,
               random = ~ year + location:year + G + location:G + 
                          year:G + location:year:G, 
               family = asr_gaussian(dispersion = 1),
               weights = weights,
               data = pheno)

M3.2 <- asreml(fixed = Value ~ location,
               random = ~ year + location:year + G + diag(location):G + 
                          year:G + location:year:G, 
              family = asr_gaussian(dispersion = 1),
              weights = weights,
              data = pheno)

M3.3 <- asreml(fixed = Value ~ location,
               random = ~ year + location:year + G + location:G + 
                          diag(year):G + location:year:G, 
               family = asr_gaussian(dispersion = 1),
               weights = weights,
               data = pheno)

M3.4 <- asreml(fixed = Value ~ location,
               random = ~ year + location:year + G + diag(location):G + 
                          diag(year):G + location:year:G, 
               family = asr_gaussian(dispersion = 1),
               weights = weights,
               data = pheno)

M3.5 <- asreml(fixed = Value ~ location,
               random = ~ year + location:year + G + fa(location,1):G + 
                          diag(year):G + location:year:G, 
                family = asr_gaussian(dispersion = 1),
                weights = weights,
                data = pheno)

## obtaining variance explained by FA for GxL from model M3.5
res <- summary(M3.5, all = TRUE)$varcomp
Load <- matrix(summary(M3.5, all = TRUE)$varcomp['component'][c((L+3+Y+1):(nrow(res)-2)),],L,1)
Psi <- diag(summary(M3.5, all = TRUE)$varcomp['component'][(3+Y+1):(3+Y+L),])
VCOV <- (Load %*% t(Load)) + Psi
Diag.Matrix <- diag(sqrt(diag(VCOV)))
Corr <- solve(Diag.Matrix) %*% VCOV %*% solve(Diag.Matrix) # Correlation matrix
Var <-(100*sum(diag(Load %*% t(Load))))/ sum(diag(Load %*% t(Load) + Psi)) # Variance explained

#### Models M3.6, M3.7, and M3.8 are coded similarly ####

M3.9 <- asreml(fixed = Value ~ location,
               random = ~  year + location:year + G + fa(location,1):G + 
                           fa(year,1):G + location:year:G, 
               family = asr_gaussian(dispersion = 1),
               weights = weights,
               data = pheno)

## obtaining variance explained by FA for both GxL and GxY

res <- summary(M3.9, all = TRUE)$varcomp

## GxY
Load <- matrix(summary(M3.95, all = TRUE)$varcomp['component'][c((3+Y+1):(3+Y+Y)),],Y,1)
Psi <- diag(summary(M3.9, all = TRUE)$varcomp['component'][(4):(3+Y),])
VCOV <- (Load %*% t(Load)) + Psi
Diag.Matrix <- diag(sqrt(diag(VCOV)))
Corr <- solve(Diag.Matrix) %*% VCOV %*% solve(Diag.Matrix)
Var <-(100*sum(diag(Load %*% t(Load))))/ sum(diag(Load %*% t(Load) + Psi))

## GxL
Load <- matrix(summary(M3.9, all = TRUE)$varcomp['component'][c((3+2*Y+L+1):(3+2*Y+L+L)),],L,1)
Psi <- diag(summary(M3.9, all = TRUE)$varcomp['component'][c((3+2*Y+1):(3+2*Y+L)),])
VCOV <- (Load %*% t(Load)) + Psi
Diag.Matrix <- diag(sqrt(diag(VCOV)))
Corr <- solve(Diag.Matrix) %*% VCOV %*% solve(Diag.Matrix)
Var <-(100*sum(diag(Load %*% t(Load))))/ sum(diag(Load %*% t(Load) + Psi))

#### Models M3.10, M3.11, ..., and M3.20 are coded similarly ####