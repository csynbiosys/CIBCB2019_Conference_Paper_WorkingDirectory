
# --------------------- GAUSSIAN MIXTURE APROXIMATION VISUAL CHECK FOR POSTERIORS --------------------- #

# Script to plot the marginal kernel densities of the stan sample results and samples from the Gaussian Mixture approximation. 
# The input value for indx is either nothing or else the string "Reparam" to load the results trom the script ReparamGaussMixt.R.
# exper takes a character string indicating the name of the experiment for which the stanfit object result will be analysed. 
# Apro is a character string to add aditional information to the plot and the saved file. The function requires
# to have run the relativeEntropyFunction or ReparamGaussMixt for the results to be analysed. 


DensAproxComp <- function(exper, apro, indxw = ""){
  
  # Load Gaussian Mixtures results from MClust
  if (!require("mclust",character.only = TRUE)){
    install.packages(x,dep=TRUE)
    if(!require("mclust",character.only = TRUE)) stop("Package not found")
  }
  library(mclust)
  
  mvpdfPosterior <- readRDS(paste(indxw, "mvpdfPostRefined_", exper, ".rds", sep = ""))
  
  
  # Extraction of parameters (14D)
  
  # Means
  mvMean <- mvpdfPosterior$parameters$mean
  # Covariances matrices
  mvCovar <- mvpdfPosterior$parameters$variance
  # Proportionality vectors
  mvPro <- mvpdfPosterior$parameters$pro
  # Optimum number of components
  gm <- mvpdfPosterior$G
  
  # Function of Density Estimation for Gaussian Mixtures
  
  GK <- function (mu, s, c){
    
    d <- seq(min(mu)-s[which.min(mu)]*4, max(mu)+s[which.max(mu)]*4, length=10000)
    
    kd <- 0
    
    
    for (x in 1:length(mu)) {
      
      k <- c[x]*((1/(s[x]*sqrt(2*pi)))*exp(-((d-mu[x])^2)/(2*s[x]^2)))
      
      kd <- kd+k
      
    }
    
    resu <- list(de = kd, le = d)
    
    return(resu)
    
  }
  
  
  
  # Extraction and check of the marginal posteriors
  
  ###################################### K_IPTG
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[1]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[1,1,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_k_IPTG_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(k_IPTG), main = "Kernel Density Estimation for k_IPTG", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### K_aTc
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[2]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[2,2,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_k_aTc_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(k_aTc), main = "Kernel Density Estimation for k_aTc", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### k_L_pm0
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[3]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[3,3,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_k_L_pm0_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(k_L_pm0), main = "Kernel Density Estimation for k_L_pm0", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### k_L_pm
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[4]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[4,4,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_k_L_pm_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(k_L_pm), main = "Kernel Density Estimation for k_L_pm", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### theta_T
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[5]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[5,5,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_theta_T_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(theta_T), main = "Kernel Density Estimation for theta_T", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### theta_aTc
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[6]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[6,6,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_theta_aTc_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(theta_aTc), main = "Kernel Density Estimation for theta_aTc", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### n_aTc
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[7]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[7,7,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_n_aTc_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(n_aTc), main = "Kernel Density Estimation for n_aTc", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### n_T
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[8]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[8,8,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_n_T_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(n_T), main = "Kernel Density Estimation for n_T", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### k_T_pm0
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[9]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[9,9,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_k_T_pm0_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(k_T_pm0), main = "Kernel Density Estimation for k_T_pm0", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### k_T_pm
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[10]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[10,10,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_k_T_pm_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(k_T_pm), main = "Kernel Density Estimation for k_T_pm", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### theta_L
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[11]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[11,11,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_theta_L_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(theta_L), main = "Kernel Density Estimation for theta_L", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### theta_IPTG
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[12]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[12,12,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_theta_IPTG_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(theta_IPTG), main = "Kernel Density Estimation for theta_IPTG", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### n_IPTG
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[13]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[13,13,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_n_IPTG_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(n_IPTG), main = "Kernel Density Estimation for n_IPTG", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  ###################################### n_L
  mu1 <- c()
  
  for(x in 1:gm){
    m <- mvMean[,x][[14]]
    mu1 <- c(mu1, m)
  }
  
  va1 <- c()
  
  for(x in 1:gm){
    m <- mvCovar$sigma[14,14,x]
    va1 <- c(va1, m)
  }
  sd1 <- sqrt(va1)
  
  tot <- c()
  for(x in 1:gm){
    s <- rnorm(10000*mvPro[x], mu1[x], sd1[x])
    tot <- c(tot, s)
  }
  
  png(paste("GMM_n_L_", exper, "_", apro, ".png", sep = ""), width = 6, height = 4, units = 'in', res = 300)
  plot(density(n_L), main = "Kernel Density Estimation for n_L", col = "red")
  lines(density(tot), col = "green")
  dev.off()
  
  
  
}

