start = Sys.time()

#############################################
####### Set up age, corstr_sev, kappa ######
#############################################

#Read the age and corstr_sev
args = commandArgs(trailingOnly=TRUE) #read in age (9,13,17,23) as needed for time-specific results
age = as.numeric(args[1])
corstr_sev = args[2]
#kappa = 1

######################
### Load packages ####
######################

library(writexl)
library(readxl)

#######################
###### Set path #######
#######################

# Check if running in a SLURM environment
if (!is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))) {
  # If in SLURM environment
  setwd("/blue/somnath.datta/shoumisarkar/Fluorosis/")
} else {
  # If not in SLURM environment
  setwd("W:/somnath.datta/shoumisarkar/Fluorosis/")
}

source(file = "Codes/functions.R")

filename = paste0("Results/02_preprocess_data/preprocessed_IFS_data_age", age, ".xlsx")
IFS_dat <- read_xlsx(filename)

#######################################
####### Updation of beta_p ############
#######################################

fit = fit_GEE_severity(dat=IFS_dat, corstr_sev = corstr_sev)

coef_sev = fit$coefs
filename = paste0("Results/04_severity_modelling/", corstr_sev, "/whole_data_based/coefs_sev_age", age, ".RData")
save(coef_sev, file = filename) 

iter=fit$iter
filename = paste0("Results/04_severity_modelling/", corstr_sev, "/whole_data_based/iter_sev_age", age, ".RData")
save(iter, file = filename) 

rho_sev = fit$rho_sev
filename = paste0("Results/04_severity_modelling/", corstr_sev, "/whole_data_based/rho_sev_age", age, ".RData")
save(rho_sev, file = filename) 

#######################################
###### Get jackknifed estimates #######
#######################################

all_IDs = unique(IFS_dat$SUBJECT_ID)

coefs_JK_sev = matrix( , nrow=(0+0+2+(ncol(IFS_dat)-2)),
                        ncol= length(all_IDs) )

rho_sev_JK = c()
iter_JK = c()

#Jackknifing to calculate SDs...

for(id in all_IDs)
{
  JK_subdat = IFS_dat %>% filter(IFS_dat$SUBJECT_ID != id)  #leave cluster i out
  
  if(length(unique(JK_subdat$FRI_Score))<4)
  {
    next
  }
  
  temp_fit = fit_GEE_severity(dat = JK_subdat, corstr_sev = corstr_sev)
  temp_coefs = temp_fit$coefs
  ind = which(all_IDs == id)
  coefs_JK_sev[,ind] = temp_coefs[,2]
  
  if(corstr_sev!="jackknifed")
  {
    rho_sev_JK = c(rho_sev_JK, temp_fit$rho_sev)
  }
  
  iter_JK = c(iter_JK, temp_fit$iter) 
  
  print(paste0("Jackknifing... ", ind, "/", length(all_IDs)))
  
  rownames(coefs_JK_sev) = temp_coefs[,1]
  
  filename = paste0("Results/04_severity_modelling/", corstr_sev, "/whole_data_based/coefs_JK_sev_age", age, ".RData")
  save(coefs_JK_sev, file = filename)
  
}

#For ar1/exchangeable, save the jackknifed rho's...
if(corstr_sev!="jackknifed")
{ 
  filename = paste0("Results/04_severity_modelling/", corstr_sev, "/whole_data_based", 
                    "/rho_sev_JK_age", age,".RData")
  save(rho_sev_JK, file = filename)
  
}

##################################################
###### Get SD of the jackknifed estimators #######
##################################################

JK_SD_sev = data.frame(apply(coefs_JK_sev, MARGIN=1, FUN= compute_JK_SE))
colnames(JK_SD_sev) = "SD"

#Write the file:
filename = paste0("Results/04_severity_modelling/", corstr_sev, "/whole_data_based/JK_SD_sev_age", age, ".RData")
save(JK_SD_sev, file = filename)

##################################################
#########  Standardize the estimators:   #########
##################################################

std_coef_sev = coef_sev[,2]/JK_SD_sev
colnames(std_coef_sev) = "Standardized Estimate"


#Write the file:
filename = paste0("Results/04_severity_modelling/", corstr_sev, "/whole_data_based/std_coefs_sev_age", age, ".RData")
save(std_coef_sev, file = filename)


end = Sys.time()
elapsed = end - start; elapsed

