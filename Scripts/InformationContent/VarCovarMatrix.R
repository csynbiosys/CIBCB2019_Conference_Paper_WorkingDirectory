
# ---------------------- COVARIANCE MATRIX CALCULATION ---------------------- #

# Script to compute the determinant of the covariance matrix for a set of stanfit objects (experiment) and saving the results as a csv file. 

# String vector with the stanfit objects generated in this work
experiment <- c("Calibration_4","Calibration_5","Calibration_6",
                 "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_8", "DynStim_9",
                 "DynStim_11", "DynStim_14", "ALL_Long_Model3.stan")

varcovarDet <- function(experiment){

  srdvc <<- matrix(data=NaN, nrow=length(experiment), ncol=1)
  
  for(y in 1:length(experiment)){
    
    # Extraction of samples from stanfit object for each parameter
    x <- as.array(readRDS(paste("fit_", experiment[y], ".rds", sep="")))
    
    k_IPTG <- c(x[,1,15], x[,2,15], x[,3,15], x[,4,15])
    k_aTc <- c(x[,1,16], x[,2,16], x[,3,16], x[,4,16])
    k_L_pm0 <- c(x[,1,17], x[,2,17], x[,3,17], x[,4,17])
    k_L_pm <- c(x[,1,18], x[,2,18], x[,3,18], x[,4,18])
    theta_T <- c(x[,1,19], x[,2,19], x[,3,19], x[,4,19])
    theta_aTc <- c(x[,1,20], x[,2,20], x[,3,20], x[,4,20])
    n_aTc <- c(x[,1,21], x[,2,21], x[,3,21], x[,4,21])
    n_T <- c(x[,1,22], x[,2,22], x[,3,22], x[,4,22])
    k_T_pm0 <- c(x[,1,23], x[,2,23], x[,3,23], x[,4,23])
    k_T_pm <- c(x[,1,24], x[,2,24], x[,3,24], x[,4,24])
    theta_L <- c(x[,1,25], x[,2,25], x[,3,25], x[,4,25])
    theta_IPTG <- c(x[,1,26], x[,2,26], x[,3,26], x[,4,26])
    n_IPTG <- c(x[,1,27], x[,2,27], x[,3,27], x[,4,27])
    n_L <- c(x[,1,28], x[,2,28], x[,3,28], x[,4,28])
    
    # Add results of each parameter for one experiment into one same matrix
    M <- cbind(k_IPTG, k_aTc, k_L_pm0, k_L_pm, theta_T, theta_aTc, n_aTc, n_T, k_T_pm0, k_T_pm, theta_L, theta_IPTG, n_IPTG, n_L)
    # Computation of the variance-covariance matrix
    varcovar <- cov(M)
    # Computation of the determinant for the covariance matrix
    res <- (det(varcovar))
    
    srdvc[y,1] = res
    # Write csv file with the results
    write.table(srdvc, file = "DetCovarM.csv",row.names=experiment, na="", sep=",")

  }
  
}







































































