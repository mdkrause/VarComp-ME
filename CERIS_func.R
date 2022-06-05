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

########################################################################################################################
##                                                  DISCLAIMER                                                        ##
##                 We basically removed all the plotting options and added parallel computing                         ##   
## The original one can be founded at: https://github.com/jmyu/CERIS-JGRA-Workshop/blob/main/Sub_functions_Workshop.r ##
########################################################################################################################


Exhaustive_search <- function(env_mean_trait, env_paras, searching_daps, FTdaps, dap_x, dap_y, LOO, Paras, window) {
        
        nParas <- length(Paras);
       
                dap_win <- searching_daps * searching_daps  / 2;
                
                pop_cors_matrix <- matrix(ncol = 3 + (2 * nParas), nrow = dap_win * 1);
                
                colnames(pop_cors_matrix) <- c('Day_x', 'Day_y', 'window', paste('R_', Paras, sep = ''), 
                                               paste('nR_', Paras, sep = ''));
                
                results <-foreach(d1 = 1:(dap_y - window), .combine = "rbind")%:%
                        
                        foreach(d2= (d1 + window):dap_y, .combine = "rbind")  %dopar% {
                        
                                days <- c(d1:d2); 
                                env_facts_matrix <- matrix(nrow = nrow(env_mean_trait), ncol = nParas);
                                
                                
                                for (e_i in 1:nrow(env_mean_trait)) {
                                        e <- env_mean_trait$env_code[e_i];
                                        env_para <- subset(env_paras, env_paras$env_code == e);
                                        env_mean <- colMeans(env_para[days, Paras], na.rm = T);
                                        env_facts_matrix[e_i,] <- env_mean;
                                }
                                
                                #n <- n + 1;
                                ### leave one environment out and get the median correlation
                                Ymean_envPara <- cbind(env_facts_matrix, env_mean_trait$meanY);
                                rs <- c();
                                if (LOO == 0) {
                                        for (k in 1:nParas) {
                                                rs[k] <- round(cor(Ymean_envPara[,nParas + 1], Ymean_envPara[,k]), digits = 4)
        
                                        }
                                } else {
                                        loo_rs_matrix <- matrix(nrow = nrow(Ymean_envPara)+ 0, ncol = nParas);
                                        for (k in 1:nParas) { 
                                                for (e_x in c(1:nrow(Ymean_envPara))) {
                                                        t_matrix <- Ymean_envPara[-e_x,];
                                                        loo_rs_matrix[e_x, k] <- round(cor(t_matrix[,nParas + 1], t_matrix[,k]), digits = 4)
                                                }
                                        }
                                        rs <- apply(loo_rs_matrix, 2, median);
                                }
                                pop_cors_matrix <- c(d1, d2, d2 - d1, rs, 0 - rs);
                                return(pop_cors_matrix)
                        }
        
                colnames(results) <- c('Day_x', 'Day_y', 'window', paste('R_', Paras, sep = ''), paste('nR_', Paras, sep = ''));
               
               return(results)                
}