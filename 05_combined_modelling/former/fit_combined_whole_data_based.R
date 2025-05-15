start = Sys.time()

########################################
####### Set up age, corstr, kappa ######
########################################


#Read the age and corstr
args = commandArgs(trailingOnly=TRUE) #read in age (9,13,17,23) as needed for time-specific results
age = as.numeric(args[1])
corstr = args[2]

######################
### Load packages ####
######################

library(writexl)
library(readxl)
library(geepack)
library(MASS)
library(dplyr)

#######################
###### Set path #######
#######################

inSLURM = !is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))

# Check if running in a SLURM environment
if (inSLURM) {
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

#Extract gamma and coefs from separate modelling:

path_pres = paste0(ifelse(inSLURM, "/blue/", "W:/"), 
                   "somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/", 
                   corstr, "/whole_data_based")

path_sev = paste0(ifelse(inSLURM, "/blue/", "W:/"), 
                  "somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/", 
                  corstr, "/whole_data_based")


obj_pres = load(paste0(path_pres,"/coeffs_pres_age", age, ".RData"))
assign("coef_pres", get(obj_pres))

obj_sev = load(paste0(path_sev,"/coeffs_sev_age", age, ".RData"))
assign("coef_sev", get(obj_sev))

fit = fit_GEE(dat=IFS_dat, corstr_pres=corstr, corstr_sev=corstr, 
              init_pres_coef = coef_pres$Estimates, init_sev_coef = coef_sev$Estimates)

coef = fit$coefs

coef_pres = coef
filename = paste0("Results/05a_combined_presence_modelling/", corstr, "/whole_data_based/coeffs_pres_age", age, ".RData")
save(coef_pres, file = filename) 

coef_sev = coef
coef_sev[-c(1:4),2] = coef_sev[1,2] * coef_sev[-c(1:4),2] 
filename = paste0("Results/05b_combined_severity_modelling/", corstr, "/whole_data_based/coeffs_sev_age", age, ".RData")
save(coef_sev, file = filename) 


rho_pres = fit$rho_pres
filename = paste0("Results/05a_combined_presence_modelling/", corstr, "/whole_data_based/rho_pres_age", age, ".RData")
save(rho_pres, file = filename)

rho_sev = fit$rho_sev
filename = paste0("Results/05b_combined_severity_modelling/", corstr, "/whole_data_based/rho_sev_age", age, ".RData")
save(rho_sev, file = filename)

#####################################################


#######################################
###### Get jackknifed estimates #######
#######################################

obj_pres_JK = load(paste0(path_pres,"/coeffs_JK_pres_age", age, ".RData"))
assign("pres_JK", get(obj_pres_JK))

obj_sev_JK = load(paste0(path_sev,"/coeffs_JK_sev_age", age, ".RData"))
assign("sev_JK", get(obj_sev_JK))

coefs_JK_pres = matrix( , nrow=(1+1+2+(ncol(IFS_dat)-2)),
                        ncol= length(unique(IFS_dat$SUBJECT_ID)) )

coefs_JK_sev = matrix( , nrow=(1+1+2+(ncol(IFS_dat)-2)),
                       ncol= length(unique(IFS_dat$SUBJECT_ID)) )

rho_pres_JK = c()
rho_sev_JK = c()


for(id in unique(IFS_dat$SUBJECT_ID))
{
  print(id)
  
  JK_subdat = IFS_dat %>% filter(IFS_dat$SUBJECT_ID != id)  #leave cluster i out
  
  if(length(unique(JK_subdat$FRI_Score))<4)
  {
    next
  }
  
  JK_ind = which(unique(IFS_dat$SUBJECT_ID) %in% id)
  
  pres_coef_JK = pres_JK[, JK_ind]
  sev_coef_JK = sev_JK[, JK_ind]
  
  temp_fit = fit_GEE(dat = JK_subdat, corstr_pres=corstr, corstr_sev=corstr, 
                     init_pres_coef = pres_coef_JK, init_sev_coef = sev_coef_JK) 
  
  temp_coefs = temp_fit$coefs
  
  ind = which(unique(IFS_dat$SUBJECT_ID) == id)
  
  coefs_JK_pres[,ind] = temp_coefs[,2] 
  
  coefs_JK_sev[,ind] = temp_coefs[ ,2] 
  coefs_JK_sev[-c(1:4),ind] = coefs_JK_sev[1,ind] * coefs_JK_sev[-c(1:4),ind]
  
  rho_pres_JK = c(rho_pres_JK, temp_fit$rho_pres)
  rho_sev_JK = c(rho_sev_JK, temp_fit$rho_sev)
  
  print(paste0(ind, "/", length(unique(IFS_dat$SUBJECT_ID))))
  
  rownames(coefs_JK_pres) = rownames(coefs_JK_sev) = as.vector(temp_coefs[,1])
  
  filename = paste0("Results/05a_combined_presence_modelling/", corstr, "/whole_data_based/coeffs_JK_pres_age", age, ".RData")
  save(coefs_JK_pres, file = filename)
  
  filename = paste0("Results/05b_combined_severity_modelling/", corstr, "/whole_data_based/coeffs_JK_sev_age", age, ".RData")
  save(coefs_JK_sev, file = filename)
}

##################################################
###### Get SD of the jackknifed estimators #######
##################################################

JK_SD_pres = data.frame(apply(coefs_JK_pres, MARGIN=1, FUN= compute_JK_SE))
colnames(JK_SD_pres) = "SD"

JK_SD_sev = data.frame(apply(coefs_JK_sev, MARGIN=1, FUN= compute_JK_SE))
colnames(JK_SD_sev) = "SD"

#Write the file:

filename = paste0("Results/05a_combined_presence_modelling/", corstr, "/whole_data_based/JK_SD_pres_age", age, ".RData")
save(JK_SD_pres, file = filename)

filename = paste0("Results/05b_combined_severity_modelling/", corstr, "/whole_data_based/JK_SD_sev_age", age, ".RData")
save(JK_SD_sev, file = filename)


##################################################
#########  Standardize the estimators:   #########
##################################################

std_coef_pres = coef_pres[,2]/JK_SD_pres

std_coef_sev = coef_sev[,2]/JK_SD_sev

#Write the file:

#filename = paste0("Results/05a_combined_presence_modelling/", corstr, "/whole_data_based/std_coeffs_age", age, ".xlsx")
#write_xlsx(std_coef, path = filename)

filename = paste0("Results/05a_combined_presence_modelling/", corstr, "/whole_data_based/std_coeffs_pres_age", age, ".RData")
save(std_coef_pres, file = filename)

filename = paste0("Results/05b_combined_severity_modelling/", corstr, "/whole_data_based/std_coeffs_sev_age", age, ".RData")
save(std_coef_sev, file = filename)


filename = paste0("Results/05a_combined_presence_modelling/", corstr, "/whole_data_based", 
                  "/rho_pres_JK_age", age,".RData")
save(rho_pres_JK, file = filename)

filename = paste0("Results/05b_combined_severity_modelling/", corstr, "/whole_data_based", 
                  "/rho_sev_JK_age", age,".RData")
save(rho_sev_JK, file = filename)


end = Sys.time()
elapsed = end - start; elapsed

