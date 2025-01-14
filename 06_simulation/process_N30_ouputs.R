
start = Sys.time()

#Load libraries:
library(readxl)
library(writexl)
library(dplyr)

N = 30
mode = "severity" #"severity" #for mode="combined", edits needed: need to specify the pres and sev files in the file reading step of the loop.

#mc_seed = 1 #change this as needed
#age = 9 #change this as needed

corstr_pres = "exchangeable"
corstr_sev = "independence"

#exch_rho_pres = 0.6 #0.3 #0.6
#exch_rho_sev = 0.6 #0.8 #0.6

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

mc_std_coefs = matrix(, nrow = nparam, ncol = 600)
rownames(mc_std_coefs) = rownames(std_coef_obj)

mc_std_coefs_list = list(mc_std_coefs, mc_std_coefs, mc_std_coefs, mc_std_coefs)
names(mc_std_coefs_list) = paste0("age", ages)


for(mc_seed in 1:600)
{
  print(mc_seed)
  
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
  
    
  }
}

# Counts how many sets of valid MC samples generated estimates per age group.
non_na_columns_counts <- lapply(mc_std_coefs_list, function(df) {
  sum(apply(df, 2, function(x) any(!is.na(x))))
})

non_na_columns_counts



row_means_list <- lapply(mc_std_coefs_list, function(df) {
  rowMeans(df, na.rm = TRUE)  # na.rm = TRUE to remove NA values from the calculation
})



row_sds_list <- lapply(mc_std_coefs_list, function(df) {
  apply(df, 1, sd, na.rm = TRUE)  # Applying sd function across rows (MARGIN = 1)
})

row_means_list; row_sds_list


non_na_columns_counts
