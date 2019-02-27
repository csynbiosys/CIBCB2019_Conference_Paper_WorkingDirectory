
# ------------------------- PARAMETER SAMPLES EXTRACTION ------------------------- #

# Function to extract all the samples for a specific parameter from a stanfit object to generate a csv file with the results
# all toghether to be used for the boxplot. The input exper takes a string vector with the names of the stanfit objects 
# to consider and exper2 the name a string vector with the names for the experiments considered. parNam is a string caracter
# containing the parameter name to be extracted and parIn the index position for the parameter samples in the stanfit
# object. Order of parameters is as specified in the STAN model script, where index from 0 to 14 are the values in the
# unit normal mapping and 15 to 28 the samples after the reparameterisation. 


experiment <- c("fit_Calibration_4","fit_Calibration_5","fit_Calibration_6",
                "fit_DynStim_1", "fit_DynStim_2", "fit_DynStim_3", "fit_DynStim_8","fit_DynStim_9", 
                "fit_DynStim_11", "fit_DynStim_14", "fit_ALL_Long_Model3.stan")
experiment2 <- c("EC1","EC2","EC3",
                 "ED1", "ED2", "ED3", "ED4","ED5", 
                 "ED6", "ED7", "Multi")

ConfidInter2 <- function (exper, exper2, parNam, parIn){
  
  k_inter <- c()
  matX <- matrix(data=0, nrow = 512, ncol = 11)
  matY <- matrix(data=0, nrow = 512, ncol = 11)
  nameEx <- c()
  sampEx <- c()
  allRES <- matrix(data = 0, nrow = 8000*length(exper), ncol = 2)
  
  for (ex in 1:length(exper)){
    
    x <- as.array(readRDS(paste(exper[ex], ".rds", sep="")))
    
    interm <- c(x[,1,parIn], x[,2,parIn], x[,3,parIn], x[,4,parIn])
    
    k_inter <- cbind(k_inter, interm)
    sampEx <- c(sampEx, interm)
    nameEx <- c(nameEx, rep(exper2[ex], 8000))
    
  }
  
  allRES[,1] <- nameEx
  allRES[,2] <- sampEx
  
  # File with all results in one same column
  n1 <- paste("Samples_Toghether_", parNam, ".csv", sep = "")
  write.table(allRES, file = n1, sep = ",")
  
  colnames(k_inter) = exper2
  # File with all results in different columns
  n2 <- paste(parNam, "_Samples.csv", sep = "")
  write.table(k_inter, file = n2, sep=",")
  
}

# Run comand for the script. Comment if necessary. 
ConfidInter2(experiment, experiment2)










