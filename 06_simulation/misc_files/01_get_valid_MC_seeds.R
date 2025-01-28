
start = Sys.time()

#Load libraries:
library(readxl)
library(writexl)
library(dplyr)


get_valid_seeds = function(N, mode, n_MC, corstr_pres, corstr_sev)
{
  ages = c(9,13,17,23)
  
  inSLURM = !is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))
  
  path_prefix = "W:/"
  
  # Check if running in a SLURM environment
  if (inSLURM) {
    # If in SLURM environment
    path_prefix = "/blue/"
  } 
  
  setwd(paste0(path_prefix, "somnath.datta/shoumisarkar/Fluorosis/"))
  
  source("Codes/functions.R")
  
  #############################################
  ############ Set age & mc_seed ##############
  #############################################
  
  age = 9
  mc_seed = 1
  
  
  
  setwd(paste0(path_prefix, "somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode, "/",
               corstr_pres, ",", corstr_sev, "/age", age))
  
  file = paste0("std_coef_pres_MC_", mc_seed, ".Rdata")
  
  if(mode=="severity"){
    file = paste0("std_coef_sev_MC_", mc_seed, ".Rdata")
  }
  
  obj = load(file)
  assign("std_coef_obj", get(obj))
  
  nparam = length(as.vector(std_coef_obj[,1]))
  
  #as.vector(std_coef_obj[,1])
  
  mc_std_coefs = matrix(, nrow = nparam, ncol = n_MC)
  rownames(mc_std_coefs) = rownames(std_coef_obj)
  
  mc_std_coefs_list = list(mc_std_coefs, mc_std_coefs, mc_std_coefs, mc_std_coefs)
  names(mc_std_coefs_list) = paste0("age", ages)
  
  
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_MC, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  
  for(mc_seed in 1:n_MC)
  {
    
    
    for(age in ages)
    {
      setwd(paste0(path_prefix, "somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode, "/",
                   corstr_pres, ",", corstr_sev, "/age", age))
      
      file = paste0("std_coef_pres_MC_", mc_seed, ".Rdata")
      
      if(mode=="severity"){
        file = paste0("std_coef_sev_MC_", mc_seed, ".Rdata")
      }
      
      if(!file.exists(file)) #if file does not exist, skip it
      {
        next
      }
      
      obj = load(file)
      assign("std_coef_obj", get(obj))
      
      age_ind = which(ages %in% age)
      
      mc_std_coefs_list[[age_ind]][,mc_seed] = as.vector(std_coef_obj[,1])
      
      
      setTxtProgressBar(pb, mc_seed) #progress bar
      
    }
  }
  
  # Counts how many sets of valid MC samples generated estimates per age group.
  non_na_columns_counts <- lapply(mc_std_coefs_list, function(df) {
    sum(apply(df, 2, function(x) any(!is.na(x))))
  })
  
  non_na_columns_counts
  
  
  #Do we have at least 100 valid samples?
  flags = unlist(lapply(non_na_columns_counts, function(x){x>100}))
  valid = !(any(!flags)) #if TRUE, we have at least 100 samples.
  
  ################################################################################
  
  if(valid==FALSE)
  {
    print("Not enough MC samples...")
    
    print(non_na_columns_counts)
  }else{
    
    # Collect the seeds for the valid MC samples per age group.
    non_na_MC_seeds <- lapply(mc_std_coefs_list, function(df) {
      which(apply(df, 2, function(x) any(!is.na(x))))
    })
    
    
    setwd(paste0(path_prefix, "somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation"))
    
    filename = paste0("validMCseeds_", corstr_pres, "_", corstr_sev, "_N", N, "_", mode, ".RData")
    save(non_na_MC_seeds, file = filename)
    
  }
}


get_valid_seeds(N = 30, mode = "presence", n_MC = 700, corstr_pres = "exchangeable", corstr_sev = "independence")
get_valid_seeds(N = 30, mode = "severity", n_MC = 700, corstr_pres = "exchangeable", corstr_sev = "independence")

get_valid_seeds(N = 50, mode = "presence", n_MC = 400, corstr_pres = "exchangeable", corstr_sev = "independence")
get_valid_seeds(N = 50, mode = "severity", n_MC = 400, corstr_pres = "exchangeable", corstr_sev = "independence")

#add more function calls as more corstrs are covered...



# row_means_list <- lapply(mc_std_coefs_list, function(df) {
#   rowMeans(df, na.rm = TRUE)  # na.rm = TRUE to remove NA values from the calculation
# })
# 
# 
# 
# row_sds_list <- lapply(mc_std_coefs_list, function(df) {
#   apply(df, 1, sd, na.rm = TRUE)  # Applying sd function across rows (MARGIN = 1)
# })
# 
# row_means_list; row_sds_list
# 
# 
# non_na_columns_counts