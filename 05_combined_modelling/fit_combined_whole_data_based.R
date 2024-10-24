#############################################
####### Set up age, corstr_sev, kappa ######
#############################################

#Read the age and corstr_sev
args = commandArgs(trailingOnly=TRUE) #read in age (9,13,17,23) as needed for time-specific results
age = as.numeric(args[1])
corstr_pres = args[2]
corstr_sev = args[3]
#kappa = 1

######################
### Load packages ####
######################

library(writexl)
library(readxl)

start = Sys.time()

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

###########################################################
####### Read in estimates from separate models ############
###########################################################

path_pres = paste0(ifelse(inSLURM, "/blue/", "W:/"), 
                   "somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/", 
                   corstr_pres, "/whole_data_based")

path_sev = paste0(ifelse(inSLURM, "/blue/", "W:/"), 
                  "somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/", 
                  corstr_sev, "/whole_data_based")


obj_pres = load(paste0(path_pres,"/coefs_pres_age", age, ".RData"))
assign("coef_pres", get(obj_pres))

obj_sev = load(paste0(path_sev,"/coefs_sev_age", age, ".RData"))
assign("coef_sev", get(obj_sev))


######################################################################
####### Read in jackknifed estimates from separate models ############
######################################################################

obj_pres_jk = load(paste0(path_pres,"/coefs_JK_pres_age", age, ".RData"))
assign("coef_pres_jk", get(obj_pres_jk))

obj_sev_jk = load(paste0(path_sev,"/coefs_JK_sev_age", age, ".RData"))
assign("coef_sev_jk", get(obj_sev_jk))


#gamma = mean(coef_sev$Estimates[-c(1,2)]/coef_pres$Estimates[-c(1)])

fit = fit_GEE_combined(dat=IFS_dat, corstr_pres = corstr_pres, corstr_sev = corstr_sev,
                       init_pres_coef = coef_pres$Estimates, init_sev_coef = coef_sev$Estimates
                       )

#Save presence results...
combined_coef_pres = fit$pres_coefs
filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev, 
                  "/whole_data_based/coefs_pres_age", age, ".RData")
save(combined_coef_pres, file = filename) 

iter=fit$iter
filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev, 
                  "/whole_data_based/iter_pres_age", age, ".RData")
save(iter, file = filename) 

#Save severity results...
combined_coef_sev = fit$sev_coefs
filename = paste0("Results/05b_combined_severity_modelling/", corstr_pres, ",", corstr_sev, 
                  "/whole_data_based/coefs_sev_age", age, ".RData")
save(combined_coef_sev, file = filename) 


## as this is already saved in the 05a_combined_presence_... folder, 
## we do not need to save a duplicate for the common `iter` in the 05b_combined_severity_... folder.
#iter=fit$iter
#filename = paste0("Results/05b_combined_severity_modelling/", corstr_sev, "/whole_data_based/iter_sev_age", age, ".RData")
#save(iter, file = filename) 


#Save the estimates of rho, if applicable
#Presence
if(corstr_pres=="ar1" || corstr_pres=="exchangeable")
{
  rho_pres = fit$rho_pres
  filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev,
                    "/whole_data_based/rho_pres_age", age, ".RData")
  save(rho_pres, file = filename) 
}

#Severity
if(corstr_sev=="ar1" || corstr_sev=="exchangeable")
{
  rho_sev = fit$rho_sev
  filename = paste0("Results/05b_combined_severity_modelling/", corstr_pres, ",", corstr_sev,
                    "/whole_data_based/rho_sev_age", age, ".RData")
  save(rho_sev, file = filename) 
}


#######################################
###### Get jackknifed estimates #######
#######################################

all_IDs = unique(IFS_dat$SUBJECT_ID)

combined_coefs_JK_pres = matrix( , nrow=(1+1+2+(ncol(IFS_dat)-2)),
                       ncol= length(all_IDs) )

combined_coefs_JK_sev  = matrix( , nrow=(1+1+2+(ncol(IFS_dat)-2)),
                        ncol= length(all_IDs) )

rho_pres_JK = c()
rho_sev_JK = c()

iter_JK = c()

#Jackknifing to calculate SDs...

for(id in all_IDs)
{
  JK_subdat = IFS_dat %>% filter(IFS_dat$SUBJECT_ID != id)  #leave cluster i out
  
  ind = which(all_IDs == id)
  
  if(length(unique(JK_subdat$FRI_Score))<4 || any(is.na(c(unlist(coef_pres_jk[,ind]), unlist(coef_sev_jk[,ind])))) ) 
  {
    next
  }
  
  temp_fit = fit_GEE_combined(dat = JK_subdat, corstr_pres = corstr_pres, corstr_sev = corstr_sev,
                     init_pres_coef = unlist(coef_pres_jk[,ind]), 
                     init_sev_coef = unlist(coef_sev_jk[,ind])
                     )
  
  temp_pres_coefs = temp_fit$pres_coefs
  combined_coefs_JK_pres[,ind] = temp_pres_coefs[,2]
  rownames(combined_coefs_JK_pres) = temp_pres_coefs[,1]
  
  temp_sev_coefs = temp_fit$sev_coefs
  combined_coefs_JK_sev[,ind] = temp_sev_coefs[,2]
  rownames(combined_coefs_JK_sev) = temp_sev_coefs[,1]
  
   if(corstr_pres=="ar1" || corstr_pres=="exchangeable")
   {
     rho_pres_JK = c(rho_pres_JK, temp_fit$rho_pres)
   }
   
   if(corstr_sev=="ar1" || corstr_sev=="exchangeable")
   {
     rho_sev_JK = c(rho_sev_JK, temp_fit$rho_sev)
   }
  
  iter_JK = c(iter_JK, temp_fit$iter) 
  
  print(paste0("Jackknifing... ", ind, "/", length(all_IDs)))
  
  filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev,
                    "/whole_data_based/coefs_JK_pres_age", age, ".RData")
  save(combined_coefs_JK_pres, file = filename)
  
  filename = paste0("Results/05b_combined_severity_modelling/", corstr_pres, ",", corstr_sev,
                    "/whole_data_based/coefs_JK_sev_age", age, ".RData")
  save(combined_coefs_JK_sev, file = filename)
  
}

#Save rho, if applicable...

if(corstr_pres=="ar1" || corstr_pres=="exchangeable")
{
  filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev, 
                    "/whole_data_based",
                    "/rho_pres_JK_age", age,".RData")
  save(rho_pres_JK, file = filename)
}

if(corstr_sev=="ar1" || corstr_sev=="exchangeable")
{
  filename = paste0("Results/05b_combined_severity_modelling/", corstr_pres, ",", corstr_sev, 
                    "/whole_data_based",
                    "/rho_sev_JK_age", age,".RData")
  save(rho_sev_JK, file = filename)

}

##################################################
###### Get SD of the jackknifed estimators #######
##################################################

JK_SD_pres = data.frame(apply(combined_coefs_JK_pres, MARGIN=1, FUN= compute_JK_SE))
colnames(JK_SD_pres) = "SD"

JK_SD_sev = data.frame(apply(combined_coefs_JK_sev, MARGIN=1, FUN= compute_JK_SE))
colnames(JK_SD_sev) = "SD"

#Save the file:
filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev, 
                  "/whole_data_based/JK_SD_pres_age", age, ".RData")
save(JK_SD_pres, file = filename)

#Save the file:
filename = paste0("Results/05b_combined_severity_modelling/", corstr_pres, ",", corstr_sev,
                  "/whole_data_based/JK_SD_sev_age", age, ".RData")
save(JK_SD_sev, file = filename)

##################################################
#########  Standardize the estimators:   #########
##################################################

std_coef_pres = combined_coef_pres[,2]/JK_SD_pres
colnames(std_coef_pres) = "Standardized Estimate"

std_coef_sev = combined_coef_sev[,2]/JK_SD_sev
colnames(std_coef_sev) = "Standardized Estimate"


#Save the file:
filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev, 
                  "/whole_data_based/std_coefs_pres_age", age, ".RData")
save(std_coef_pres, file = filename)

filename = paste0("Results/05b_combined_severity_modelling/", corstr_pres, ",", corstr_sev, 
                  "/whole_data_based/std_coefs_sev_age", age, ".RData")
save(std_coef_sev, file = filename)


end = Sys.time()
elapsed = end - start; elapsed

