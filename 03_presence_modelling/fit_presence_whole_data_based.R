start = Sys.time()

#############################################
####### Set up age, corstr_pres, kappa ######
#############################################

#Read the age and corstr_pres
args = commandArgs(trailingOnly=TRUE) #read in age (9,13,17,23) as needed for time-specific results
age = as.numeric(args[1])
corstr_pres = args[2]
#kappa = 1

######################
### Load packages ####
######################

library(writexl)
library(readxl)
#library(geepack)
#library(MASS)
#library(dplyr)

#######################
###### Set path #######
#######################

# Check if running in a SLURM environment
if (!is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))) {
  # If in SLURM environment
  setwd("/path/to/Fluorosis/")
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

fit = fit_GEE_presence(dat=IFS_dat, corstr_pres = corstr_pres, maxIter = 100)

coef_pres = fit$coefs
filename = paste0("Results/03_presence_modelling/", corstr_pres, "/whole_data_based/coefs_pres_age", age, ".RData")
save(coef_pres, file = filename) 

iter=fit$iter
filename = paste0("Results/03_presence_modelling/", corstr_pres, "/whole_data_based/iter_pres_age", age, ".RData")
save(iter, file = filename) 

rho_pres = fit$rho_pres
filename = paste0("Results/03_presence_modelling/", corstr_pres, "/whole_data_based/rho_pres_age", age, ".RData")
save(rho_pres, file = filename) 

#######################################
############ Jackknifing ##############
#######################################

all_IDs = unique(IFS_dat$SUBJECT_ID)

coefs_JK_pres = matrix( , nrow=(0+1+0+(ncol(IFS_dat)-2)),
                            ncol= length(all_IDs) )

rho_pres_JK = c()
iter_JK = c()

#Jackknifing to calculate SDs...

for(id in all_IDs)
{
  JK_subdat = IFS_dat %>% filter(IFS_dat$SUBJECT_ID != id)  #leave cluster i out
  
  if(length(unique(JK_subdat$FRI_Score))<4)
  {
    next
  }
  
  temp_fit = fit_GEE_presence(dat = JK_subdat, corstr_pres = corstr_pres, maxIter = 100)
  temp_coefs = temp_fit$coefs
  ind = which(all_IDs == id)
  coefs_JK_pres[,ind] = temp_coefs[,2]
  
  if(corstr_pres!="jackknifed")
  {
    rho_pres_JK = c(rho_pres_JK, temp_fit$rho_pres)
  }
    
  iter_JK = c(iter_JK, temp_fit$iter) 
  
  print(paste0("Jackknifing... ", ind, "/", length(all_IDs)))
  
  rownames(coefs_JK_pres) = temp_coefs[,1]
  
  filename = paste0("Results/03_presence_modelling/", corstr_pres, "/whole_data_based/coefs_JK_pres_age", age, ".RData")
  save(coefs_JK_pres, file = filename)
  
}

#For ar1/exchangeable, save the jackknifed rho's...
if(corstr_pres!="jackknifed")
{ 
  filename = paste0("Results/03_presence_modelling/", corstr_pres, "/whole_data_based", 
                    "/rho_pres_JK_age", age,".RData")
  save(rho_pres_JK, file = filename)
  
}

##################################################
############# Get the jackknifed SD ##############
##################################################

JK_SD_pres = data.frame(apply(coefs_JK_pres, MARGIN=1, FUN= compute_JK_SE))
colnames(JK_SD_pres) = "SD"

#Write the file:
filename = paste0("Results/03_presence_modelling/", corstr_pres, "/whole_data_based/JK_SD_pres_age", age, ".RData")
save(JK_SD_pres, file = filename)

##################################################
#########  Standardize the estimators:   #########
##################################################

std_coef_pres = coef_pres[,2]/JK_SD_pres
colnames(std_coef_pres) = "Standardized Estimate"


#Write the file:
filename = paste0("Results/03_presence_modelling/", corstr_pres, "/whole_data_based/std_coefs_pres_age", age, ".RData")
save(std_coef_pres, file = filename)


end = Sys.time()
elapsed = end - start; elapsed

