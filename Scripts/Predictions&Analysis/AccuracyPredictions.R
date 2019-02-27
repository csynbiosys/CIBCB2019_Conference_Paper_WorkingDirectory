
# ----------------------------- ACCURACY OF PREDICTIONS ----------------------------- #

if (!require("fpc",character.only = TRUE)){
  install.packages(x,dep=TRUE)
  if(!require("fpc",character.only = TRUE)) stop("Package not found")
}
library(fpc)

expert <- c("Calibration_4","Calibration_5","Calibration_6",
             "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_8","DynStim_9", 
             "DynStim_11", "DynStim_14")
expert2 <- c("Calibration_4","Calibration_5","Calibration_6",
             "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_8","DynStim_9", 
             "DynStim_11", "DynStim_14", "ALL_Long_Model3.stan")

# Function to compute the normalised root mean squared error (equation 7) distributions (nRMSE for each parameter vector sampled from the
# stanfit object) and save the results in csv files for post processing. As inputs it takes a string vector with the name of the
# experimental squemes to be used (experiment)  and the name a string vector with the stanfit objects to be used (experiment2)
# without the characters "fit_".

AccurPred1 <- function(experiment, experiment2){
  
  # Iteration over experimental schemes
  for(i in 1:length(experiment)){
    
    
    fileObservables <- paste(experiment[i], "_Observables.csv", sep="")
    # Extract observables data and stored into global variables
    observables <- read.csv(file=fileObservables, header=TRUE, sep=",")
    samplingT <- round(observables[,1])
    GFPmean <- observables[,2]
    GFPstd <- observables[,3]
    RFPmean <- observables[,5]
    RFPstd <- observables[,6]
    
    # Iteration over stanfit objects
    for(j in 1:length(experiment2)){
      
      # Load prediction results
      print(paste("InferenceResults_", experiment2[j], "_Simulation_", experiment[i], "_RFP.csv", sep = ""))
      C4C5r <- read.csv(paste("InferenceResults_", experiment2[j], "_Simulation_", experiment[i], "_RFP.csv", sep = ""))
      C4C5g <- read.csv(paste("InferenceResults_", experiment2[j], "_Simulation_", experiment[i], "_GFP.csv", sep = ""))
      
      RMSE1 <- c()
      RMSE2 <- c()
      
      # Compute the normalised root mean square error distributions for each parameter vector sample
      for(w in 1:length(C4C5r[1,])){
        
        inel <- ((C4C5r[,w]-RFPmean)^2)/RFPstd^2
        inel2 <- ((C4C5g[,w]-GFPmean)^2)/GFPstd^2
        
        RMSE1 <- c(RMSE1, sqrt(sum(inel)/length(C4C5r[,w])))
        RMSE2 <- c(RMSE2, sqrt(sum(inel2)/length(C4C5g[,w])))
        
      }
      
      # Save results in csv files specifying the stanfit object used and the experimental scheme used
      resnam1 <- paste("RMSE_Parameters_", experiment2[j], "_SimulationVar_", experiment[i], "_RFP.csv", sep = "")
      resnam2 <- paste("RMSE_Parameters_", experiment2[j], "_SimulationVar_", experiment[i], "_GFP.csv", sep = "")
      
      write.table(RMSE1, file = resnam1, sep=",")
      write.table(RMSE2, file = resnam2, sep=",")
      
    }
    
  }
  
}

# Execute previous function
AccurPred1(expert, expert2)


# Function to compute the Bhattacharyya distance (equation 8) between the nRMSE distribution computed with the stanfit object predicting the
# same experimental scheme as has been inferred with the nRMSE distribution of the same experimental profile predicted by the other
# posteriors from the different stanfit objects. As input it takes a vector string with the name of the different experimental profiles. 

AccurPred2 <- function(experiment){
  
  APmrLog <- matrix(data = 0, length(experiment), length(experiment))
  APmgLog <- matrix(data = 0, length(experiment), length(experiment))
  APmr <- matrix(data = 0, length(experiment), length(experiment))
  APmg <- matrix(data = 0, length(experiment), length(experiment))
  
  # Iteration over the nRMSE distributions where prediction has been done to the experimental profile that the stanfit object has been 
  # infered from.
  for(i in 1:length(experiment)){
    
    normTer1 <- read.csv(paste("RMSE_Parameters_", experiment[i], "_SimulationVar_", experiment[i], "_RFP.csv", sep = ""))
    normTer2 <- read.csv(paste("RMSE_Parameters_", experiment[i], "_SimulationVar_", experiment[i], "_GFP.csv", sep = ""))
    
    # Iteration over the nRMSE distributions for different experimental profile predictions j for each posterior parameter i used
    for(j in 1:length(experiment)){
      
      # Calculation of the Bhattacharyya distance (equation 8) between elements from first loop and each element from second loop
      C4C5r <- read.csv(paste("RMSE_Parameters_", experiment[i], "_SimulationVar_", experiment[j], "_RFP.csv", sep = ""))
      klr <- bhattacharyya.dist(mean(C4C5r[,1]), mean(normTer1[,1]), diag(1)*std(C4C5r[,1])^2, diag(1)*std(normTer1[,1])^2)
      
      C4C5g <- read.csv(paste("RMSE_Parameters_", experiment[i], "_SimulationVar_", experiment[j], "_GFP.csv", sep = ""))
      klg <- bhattacharyya.dist(mean(C4C5g[,1]), mean(normTer2[,1]), diag(1)*std(C4C5g[,1])^2, diag(1)*std(normTer2[,1])^2)
      
      # Introduce logarithm of results (for better visualisation) where the diagonal elements stay as zero
      if(klr == 0){
        APmrLog[i,j] <- (klr)
      } else {
        APmrLog[i,j] <- log(klr)
      }
      
      if(klg == 0){
        APmgLog[i,j] <- (klg)
      } else {
        APmgLog[i,j] <- log(klg)
      }
      
      # Introduce results into a matrix
      APmr[i,j] <- klr
      APmg[i,j] <- klg
      
    }
    
  }
  
  # Save all the resultant matrices as CSV files
  resnam1 <- paste("AccuracyMatrix_RFP.csv", sep = "")
  resnam2 <- paste("AccuracyMatrix_GFP.csv", sep = "")
  resnam3 <- paste("AccuracyMatrix_RFP_Log.csv", sep = "")
  resnam4 <- paste("AccuracyMatrix_GFP_Log.csv", sep = "")
  
  write.table(APmr, file = resnam1, sep=",", row.names = experiment, col.names = experiment)
  write.table(APmg, file = resnam2, sep=",", row.names = experiment, col.names = experiment)
  write.table(APmrLog, file = resnam3, sep=",", row.names = experiment, col.names = experiment)
  write.table(APmgLog, file = resnam4, sep=",", row.names = experiment, col.names = experiment)
  
}

# Execute previous function
AccurPred2(expert)


