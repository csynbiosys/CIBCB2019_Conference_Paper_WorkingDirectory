

# --------------------------------- RELATIVE ENTROPY FUNCTION --------------------------------- #

# Calculation of the entropy of prior and posterior results for an inference result using the Taylor series approximation 
# for Gaussian Mixtures proposed in "On Entropy Approximation for Gaussian Mixture Random Vectors". For all entropy calculations, 
# x represents a random vector, MU the vectors of means for the Gaussian mixtures, E the covariance matrices for the Gaussian
# mixtures and w the vector of weights for the Gaussian Mixtures. As inputs the function takes a character string with the 
# name of the experimental profile to be analysed. 

RelativeEntropy <- function(experiment){
  
  #### ENTROPY OF PRIORS ####
  
  
  # Parameters of the joint prior multivariate Gaussian
  
  MU_pri <- matrix(c(0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0),
                   nrow = 14,ncol = 1)
  
  sd_pri <- c(1,1,1,1,1,
              1,1,1,1,1,1,1,1,1)
  
  vr_pri <- (sd_pri)^2
  
  E_pri <- diag(14)*vr_pri
  
  
  # Entropy Calculation Function for the priors (analytical form for a multivariate Gaussian)
  
  HxPri <- function(E){
    
    dim <- length(E[,1])
    
    H <- (dim/2)*log(2*pi*exp(1))+0.5*log(det(E))
    
    return(H)
    
  }
  
  # Entropy for the multivariate prior
  H_prior <<- HxPri(E_pri)
  
  
  #### ENTROPY OF POSTERIORS ####
  
  
  # Extraction of all posterior samples from MCMC (STAN) results
  
  fitName <- paste("fit_", experiment, ".rds", sep="") 
  x <- as.array(readRDS(fitName))
  
  k_IPTG <- c(x[,1,1], x[,2,1], x[,3,1], x[,4,1])
  k_aTc <- c(x[,1,2], x[,2,2], x[,3,2], x[,4,2])
  k_L_pm0 <- c(x[,1,3], x[,2,3], x[,3,3], x[,4,3])
  k_L_pm <- c(x[,1,4], x[,2,4], x[,3,4], x[,4,4])
  theta_T <- c(x[,1,5], x[,2,5], x[,3,5], x[,4,5])
  theta_aTc <- c(x[,1,6], x[,2,6], x[,3,6], x[,4,6])
  n_aTc <- c(x[,1,7], x[,2,7], x[,3,7], x[,4,7])
  n_T <- c(x[,1,8], x[,2,8], x[,3,8], x[,4,8])
  k_T_pm0 <- c(x[,1,9], x[,2,9], x[,3,9], x[,4,9])
  k_T_pm <- c(x[,1,10], x[,2,10], x[,3,10], x[,4,10])
  theta_L <- c(x[,1,11], x[,2,11], x[,3,11], x[,4,11])
  theta_IPTG <- c(x[,1,12], x[,2,12], x[,3,12], x[,4,12])
  n_IPTG <- c(x[,1,13], x[,2,13], x[,3,13], x[,4,13])
  n_L <- c(x[,1,14], x[,2,14], x[,3,14], x[,4,14])
  
  # Introduce MCMC results into a dataframe to work with 
  
  s <- matrix(data = 0, nrow = length(k_IPTG), ncol = 14)
  
  s[,1] <- k_IPTG
  s[,2] <- k_aTc
  s[,3] <- k_L_pm0
  s[,4] <- k_L_pm
  s[,5] <- theta_T
  s[,6] <- theta_aTc
  s[,7] <-n_aTc
  s[,8] <-n_T
  s[,9] <-k_T_pm0
  s[,10] <-k_T_pm
  s[,11] <-theta_L
  s[,12] <-theta_IPTG
  s[,13] <-n_IPTG
  s[,14] <-n_L
  
  y <- data.frame(s)
  
  # Estimate joint posterior as Gaussian Mixtures with MClust
  if (!require("mclust",character.only = TRUE)){
    install.packages(x,dep=TRUE)
    if(!require("mclust",character.only = TRUE)) stop("Package not found")
  }
  library(mclust)
  
  # Perform clustering in a range between 1 to 30 components and use BIC to asses the optimum value
  mvpdfPost <- densityMclust(y, G=1:30)
  
  # Save MClust results
  mcres <- paste("mvpdfPost_", experiment, ".rds", sep = "")
  saveRDS(mvpdfPost, mcres)
  
  
  # Extraction of parameters from MClust results
  
  # Split of Gaussians to 4 times the optimum value to reduce variances
  gm <- mvpdfPost$G
  mvpdfRefined <- densityMclust(y, G=gm*4)
  # Save MClust results
  mcres2 <- paste("mvpdfPostRefined_", experiment, ".rds", sep = "")
  saveRDS(mvpdfRefined, mcres2)
  
  # Means of each component
  mvMean <- mvpdfRefined$parameters$mean
  # Covariances matrices of each component
  mvCovar <- mvpdfRefined$parameters$variance$sigma
  # Proportionality vectors of each component
  mvPro <- mvpdfRefined$parameters$pro
  # Optimum number of components according to the BIC
  gm <- mvpdfRefined$G
  
  
  # Calculation of upper bound for the posterior Entropy
  
  H_Upper <- function(w, E){
    
    comp <- length(w)
    dims <- length(E[,,1][,1])
    Hu <- c()
    for (i in 1:comp){
      
      hu <- w[i]*(-log(w[i])+0.5*log(det(E[,,i])*(2*pi*exp(1))^dims)) 
      Hu <- c(Hu, hu)
      
    }
    return(sum(Hu))
  }
  
  # Upper bound for the posterior entropy
  H_U_Post <- H_Upper(mvPro, mvCovar)
  
  
  # Calculation of lower bound for the posterior Entropy
  
  # PDF of a multivariate Gaussian distribution 
  mvGauss <- function(x, MU, E){
    
    gx <- (1/sqrt(det(E)*(2*pi)^2))*exp((-1/2)*t(x-MU)%*%solve(E)%*%(x-MU))
    
    return(gx)
    
  }
  
  H_Lower <- function(w, E, MU){
    
    comp <- length(w)
    Hl <- c()
    
    for(i in 1:comp){
      
      inter <- c()
      
      for(j in 1:comp){
        
        insu <- w[j]*mvGauss(MU[,i], MU[,j], E[,,i]+E[,,j])
        inter <- c(inter, insu)
        
      }
      Hl <- c(Hl, w[i]*log(sum(inter)))
      
    }
    
    return(-sum(Hl))
  }
  
  # Lower bound for the posterior entropy
  H_L_Post <- H_Lower(mvPro, mvCovar, mvMean)
  
  
  # Function for a multivariate PDF of a Gaussian Mixture
  GaussMix <- function(x, MU, E, w){
    
    comp <- length(w)
    FX <- c()
    
    for(i in 1:comp){
      
      fx <- w[i]*(1/sqrt(det(E[,,i])*(2*pi)^2))*exp((-1/2)*t(x-MU[,i])%*%solve(E[,,i])%*%(x-MU[,i]))
      FX <- c(FX, fx)
      
    }
    
    return(sum(FX))
    
  }
  
  # Zero-order Taylor Series Expansion
  ZOTSE <- function(MU, E, w){
    
    comp <- length(w)
    ZO <- c()
    
    for(i in 1:comp){
      
      zo <- w[i]*log(GaussMix(MU[,i], MU, E, w))
      ZO <- c(ZO, zo)
      
    }
    
    return(-sum(ZO))
    
  }
  
  # Zero-order Taylor series result for the approximation of the posterior entropy
  ZO_Expan <- ZOTSE(mvMean, mvCovar, mvPro)
  
  
  # Second-order Taylor Series Expansion
  if (!require("pracma",character.only = TRUE)){
    install.packages(x,dep=TRUE)
    if(!require("pracma",character.only = TRUE)) stop("Package not found")
  }
  library(pracma)
  
  # Function to calculate the most computationaly expensive part of the Taylor Expansion component
  FMix <- function(x, MU, E, w){
    comp <- length(w)
    dims <- length(E[,,1][,1])
    
    GaussMix2 <- function(x){
      
      comp <- length(w)
      FX <- c()
      MU <- mvMean
      E <- mvCovar
      w <- mvPro
      for(i in 1:comp){
        fx <- w[i]*(1/sqrt(det(E[,,i])*(2*pi)^2))*exp((-1/2)*t(x-MU[,i])%*%solve(E[,,i])%*%(x-MU[,i]))
        FX <- c(FX, fx)
      }
      return(sum(FX))
    }
    
    F_x <- list()
    
    for(j in 1:comp){
      
      inter <- (1/GaussMix(x, MU, E, w))*(x-MU[,j])%*%t(grad(GaussMix2, x))+(x-MU[,j])%*%t(solve(E[,,j])%*%(x-MU[,j]))-diag(dims)
      
      fx <- w[j]*solve(E[,,j])%*%(inter)*as.vector(mvGauss(x , MU[,j], E[,,j]))
      
      F_x[[j]] <- fx
      
    }
    
    return(Reduce("+", F_x)/GaussMix(x, MU, E, w))
  }
  
  # Calculation of the Second-Order Taylor Series expanssion term
  SOTSE <- function(MU, E, w){
    
    comp <- length(w)
    dims <- length(E[,,1][,1])
    
    H2 <- c()
    
    for(i in 1:comp){
      print(paste(" Iteration ", i," of ", comp, sep=""))
      inter2 <- (w[i]/2)*sum(FMix(MU[,i],MU,E,w)*E[,,i])
      H2 <- c(H2, inter2)
      
    }
    
    return(sum(H2))
    
  } 
  
  # Second-order Taylor series result term for the approximation of the posterior entropy
  SO_Expan <- SOTSE(mvMean, mvCovar, mvPro)
  
  
  # Aproximation of the Entropy for the posterior
  H_posterior <- ZO_Expan-SO_Expan
  
  
  # Approximation of the relative entropy between prior and posterior
  Relative_H <- H_prior-H_posterior
  
  
  # Save results as a CSV file
  tit <- paste("EntropyResults_", experiment, ".csv", sep = "")
  
  write.table(t(c(H_prior, H_U_Post, H_L_Post, H_posterior, Relative_H)), file = tit, row.names = "Entropy", 
              col.names = c("H_Prior", "H_Upper", "H_Lower", "H_Posterior", "Relative_Entropy"), sep=",")
  
  # Save value of the entropy for the posterior as an rds file
  mcres <- paste("PosteriorEntropy", experiment, ".rds", sep = "")
  saveRDS(H_posterior, mcres)
  
}


