
# ----------------------------- SENSITIVITY FILE PROCESSING ----------------------------- #

expert <- c("Calibration_4","Calibration_5","Calibration_6",
             "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_8","DynStim_9", 
             "DynStim_11", "DynStim_14")

expert2 <- c("ALL_Long_Model3.stan", "ALL_Long_Model3.stan", "ALL_Long_Model3.stan", "ALL_Long_Model3.stan", 
             "ALL_Long_Model3.stan", "ALL_Long_Model3.stan", "ALL_Long_Model3.stan", "ALL_Long_Model3.stan", 
             "ALL_Long_Model3.stan", "ALL_Long_Model3.stan")


# Function to extrac the correct matrices from the parameter draws and predictions of the system and save them in .mat files
# to be used by the Olivia Eriksson Sensitivity index calculation published on the paper "Uncertainty quantification, 
# propagation and characterization by Bayesian analysis combined with global sensitivity analysis applied to dynamical 
# intracellular pathway model". As inputs it takes a string vector with the name of the experimental schemes used (experiment), and 
# the stan fit object without the carachters "fit_" (fitRes). Both vectors need to have the same length where in this case the combination
# expert and expert will give the singe inference result files and the combination expert and  expert2 the multi-experimental
# inference results. 

SenseDatConv <- function(experiment, fitRes){
  
  if (!require("R.matlab",character.only = TRUE)){
    install.packages(x,dep=TRUE)
    if(!require("R.matlab",character.only = TRUE)) stop("Package not found")
  }
  library(R.matlab)
  
  # Extraction of a matrix with all the samples for all the parameters for a specified stanfit object and safe the results in
  # a .mat file
  for(ind in 1:length(fitRes)){
    
    x <- as.array(readRDS(paste("fit_", fitRes[ind], ".rds", sep="")))
    draws <- matrix(data = NaN, nrow = 8000, ncol = 14)
    
    draws[,1] <- c(x[,1,15], x[,2,15], x[,3,15], x[,4,15])
    draws[,2] <- c(x[,1,16], x[,2,16], x[,3,16], x[,4,16])
    draws[,3] <- c(x[,1,17], x[,2,17], x[,3,17], x[,4,17])
    draws[,4] <- c(x[,1,18], x[,2,18], x[,3,18], x[,4,18])
    draws[,5] <- c(x[,1,19], x[,2,19], x[,3,19], x[,4,19])
    draws[,6] <- c(x[,1,20], x[,2,20], x[,3,20], x[,4,20])
    draws[,7] <- c(x[,1,21], x[,2,21], x[,3,21], x[,4,21])
    draws[,8] <- c(x[,1,22], x[,2,22], x[,3,22], x[,4,22])
    draws[,9] <- c(x[,1,23], x[,2,23], x[,3,23], x[,4,23])
    draws[,10] <- c(x[,1,24], x[,2,24], x[,3,24], x[,4,24])
    draws[,11] <- c(x[,1,25], x[,2,25], x[,3,25], x[,4,25])
    draws[,12] <- c(x[,1,26], x[,2,26], x[,3,26], x[,4,26])
    draws[,13] <- c(x[,1,27], x[,2,27], x[,3,27], x[,4,27])
    draws[,14] <- c(x[,1,28], x[,2,28], x[,3,28], x[,4,28])
    
    fn <- paste("draws_", fitRes[ind], ".mat", sep = "")
    writeMat(con=fn, x=draws)
    
  }
  
  # Extract predictions of a stanfit object result to the experimental scheme it has been inferred and save the transpose matrix
  # into an .mat file
  for(ind2 in 1:length(experiment)){
    fn1 <- paste("InferenceResults_", fitRes[ind2], "_Simulation_", experiment[ind2], "_RFP.csv", sep = "")
    fn2 <- paste("InferenceResults_", fitRes[ind2], "_Simulation_", experiment[ind2], "_GFP.csv", sep = "")
    outp <- t(read.csv(fn1)[])
    outp2 <- t(read.csv(fn2)[])
    
    if(experiment[ind2]==fitRes[ind2]){
      fm1 <- paste("outpR_", experiment[ind2], "_RFP.mat", sep = "")
      fm2 <- paste("outpG_", experiment[ind2], "_GFP.mat", sep = "")
    } else {
      fm1 <- paste("outpR_", experiment[ind2],"_", fitRes[ind2], "_RFP.mat", sep = "")
      fm2 <- paste("outpG_", experiment[ind2],"_", fitRes[ind2], "_GFP.mat", sep = "")
    }
    
    writeMat(con=fm1, x=outp)
    writeMat(con=fm2, x=outp)
    
  }
  
}
