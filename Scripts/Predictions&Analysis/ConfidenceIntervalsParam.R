
# --------------------------- PARAMETER CONFIDENCE INTERVAL PLOTS ---------------------------#

# Function to plot the 95% confidence intervals plots for each parameter from a stanfit object, where a plot for each 
# parameter is generated comparing the posteriors for each inference. The input exper takes a string vector
# with the names of the stanfit objects to consider and exper2 the name a string vector with the names for the experiments considered

# Package to load for the plots
library(bayesplot)

experiment <- c("fit_Calibration_4","fit_Calibration_5","fit_Calibration_6",
                 "fit_DynStim_1", "fit_DynStim_2", "fit_DynStim_3", "fit_DynStim_8","fit_DynStim_9", 
                 "fit_DynStim_11", "fit_DynStim_14", "fit_ALL_Long_Model3.stan")
experiment2 <- c("EC1","EC2","EC3",
                 "ED1", "ED2", "ED3", "ED4","ED5", 
                 "ED6", "ED7", "Multi")

ConfidInter <- function (exper, exper2){

  # Parameter names
  param <- c("k_IPTG", "k_aTc", "k_L_pm0", "k_L_pm", "theta_T", "theta_aTc", "n_aTc", "n_T", 
             "k_T_pm0", "k_T_pm", "theta_L", "theta_IPTG", "n_IPTG", "n_L")
  
  ################# k_IPTG #################
  
  k_IPTG_inter <- c()
  
  for (ex in 1:length(exper)){
  
    # Extract stanfit results
    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))
    # Consider only the ones for a specific parameter
    interm <- c(x[,1,15], x[,2,15], x[,3,15], x[,4,15])
    # Join the results into one same matrix
    k_IPTG_inter <- cbind(k_IPTG_inter, interm)
  }
  # Add experimental identifier names to each result
  colnames(k_IPTG_inter) = exper2
  
  ################# k_aTc #################

  k_aTc_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,16], x[,2,16], x[,3,16], x[,4,16])

    k_aTc_inter <- cbind(k_aTc_inter, interm)
  }

  colnames(k_aTc_inter) = exper2


  ################# k_L_pm0 #################

  k_L_pm0_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,17], x[,2,17], x[,3,17], x[,4,17])

    k_L_pm0_inter <- cbind(k_L_pm0_inter, interm)
  }

  colnames(k_L_pm0_inter) = exper2


  ################# k_L_pm #################

  k_L_pm_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,18], x[,2,18], x[,3,18], x[,4,18])

    k_L_pm_inter <- cbind(k_L_pm_inter, interm)
  }

  colnames(k_L_pm_inter) = exper2


  ################# theta_T #################

  theta_T_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,19], x[,2,19], x[,3,19], x[,4,19])

    theta_T_inter <- cbind(theta_T_inter, interm)
  }

  colnames(theta_T_inter) = exper2


  ################# theta_aTc #################

  theta_aTc_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,20], x[,2,20], x[,3,20], x[,4,20])

    theta_aTc_inter <- cbind(theta_aTc_inter, interm)
  }

  colnames(theta_aTc_inter) = exper2


  ################# n_aTc #################

  n_aTc_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,21], x[,2,21], x[,3,21], x[,4,21])

    n_aTc_inter <- cbind(n_aTc_inter, interm)
  }

  colnames(n_aTc_inter) = exper2


  ################# n_T #################

  n_T_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,22], x[,2,22], x[,3,22], x[,4,22])

    n_T_inter <- cbind(n_T_inter, interm)
  }

  colnames(n_T_inter) = exper2


  ################# k_T_pm0 #################

  k_T_pm0_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,23], x[,2,23], x[,3,23], x[,4,23])

    k_T_pm0_inter <- cbind(k_T_pm0_inter, interm)
  }

  colnames(k_T_pm0_inter) = exper2


  ################# k_T_pm #################

  k_T_pm_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,24], x[,2,24], x[,3,24], x[,4,24])

    k_T_pm_inter <- cbind(k_T_pm_inter, interm)
  }

  colnames(k_T_pm_inter) = exper2


  ################# theta_L #################

  theta_L_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,25], x[,2,25], x[,3,25], x[,4,25])

    theta_L_inter <- cbind(theta_L_inter, interm)
  }

  colnames(theta_L_inter) = exper2


  ################# theta_IPTG #################

  theta_IPTG_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,26], x[,2,26], x[,3,26], x[,4,26])

    theta_IPTG_inter <- cbind(theta_IPTG_inter, interm)
  }

  colnames(theta_IPTG_inter) = exper2


  ################# n_IPTG #################

  n_IPTG_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,27], x[,2,27], x[,3,27], x[,4,27])

    n_IPTG_inter <- cbind(n_IPTG_inter, interm)
  }

  colnames(n_IPTG_inter) = exper2


  ################# n_L #################

  n_L_inter <- c()

  for (ex in 1:length(exper)){

    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))

    interm <- c(x[,1,28], x[,2,28], x[,3,28], x[,4,28])

    n_L_inter <- cbind(n_L_inter, interm)
  }

  colnames(n_L_inter) = exper2

  
  
  # Add all results into one same list
  fin <- list()
  names2 <- list(k_IPTG_inter, k_aTc_inter, k_L_pm0_inter, k_L_pm_inter, theta_T_inter, theta_aTc_inter, n_aTc_inter, n_T_inter,
                 k_T_pm0_inter, k_T_pm_inter, theta_L_inter, theta_IPTG_inter, n_IPTG_inter, n_L_inter)
  
  for(y in 1:14){
    fin[param[y]] = names2[y]
  }
  
  # Generate a plot for each parameter with x in a log scale for better visualisation
  for(indx in 1:14){
    
    imn <- paste("ConfInter_", param[indx], ".png", sep = "")
    png(imn, width = 9, height = 5, units = 'in', res = 300)
    # change mcmc_intervals for mcmc_dens to generate plots with the densities
    plot(mcmc_intervals(fin[indx])+ggtitle(param[indx])+ scale_x_continuous(trans = "log"))
    dev.off()
    
  }
  
}












