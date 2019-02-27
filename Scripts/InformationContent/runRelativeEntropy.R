
# ------------------- RUN RELATIVE ENTROPY COMPUTATION ------------------- #

# Script to run the function RelativeEntropy from the script RelativeEntropyFunction.R for a specified set of stan inference results

# String vector with the stanfit object names (without the element "fit_" at the beginning) to compute the relative ntropy to

files4<-c("Calibration_4","Calibration_5","Calibration_6",
          "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_6", "DynStim_8",
          "DynStim_11", "DynStim_14", "ALL_Long_Model3.stan")

# For loop to iterate over all the stanfit object results selected, computing the relative entropy for it and saving the results into csv and rds files
for(na in files4){
  
  RelativeEntropy(na)
  
}



