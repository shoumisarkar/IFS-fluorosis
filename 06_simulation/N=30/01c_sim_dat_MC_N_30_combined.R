
#Read the MC_seed #
args = commandArgs(trailingOnly=TRUE)
mc_seed = as.numeric(args[1])
age = as.numeric(args[2])


start = Sys.time()

#Load libraries:
library(readxl)
library(dplyr)
library(MASS)

N = 30
#mc_seed = 1 #change this as needed
#age = 9 #change this as needed

corstr_pres = "jackknifed" #"independence" #"jackknifed" # "exchangeable" #"independence"
corstr_sev = "jackknifed" #"jackknifed" # "exchangeable"  #"independence" #"exchangeable"

exch_rho_pres = 0.3
exch_rho_sev = 0.8

print(paste0(corstr_pres, ",", corstr_sev, ",", exch_rho_pres, "," , exch_rho_sev))

setwd("/blue/somnath.datta/shoumisarkar/Fluorosis")
source("Codes/functions.R")


print(paste0("Running with MC seed ", mc_seed))


#Load the data for the chosen age group:
dat <- read_excel(paste0("Results/02_preprocess_data/preprocessed_IFS_data_age", age, ".xlsx"))

#Keep selected columns:
dat <- dat %>% dplyr::select(SUBJECT_ID, FRI_Score, dental_age, Avg_homeppm, 
                             #SugarAddedBeverageOzPerDay, BrushingFrequencyPerDay, 
                             Tooth8, Tooth9, Tooth10, ZoneM, ZoneI, ZoneO)


#Filter to keep only SUBJECT_IDs with full information (16 rows)
complete_dat <- dat %>%
  group_by(SUBJECT_ID) %>%
  filter(n() == 16) %>%
  ungroup()

#Get the superset of SUBJECT_IDs corresponding to complete information.
ids = unique(complete_dat$SUBJECT_ID)

#############################################################################################################################
## Create a "MC dataset" of size N by sampling with replacement from these IDs to get "N" IDs. Rename the IDs to 1,...,N  ###
#############################################################################################################################

#Sample N IDs with replacement from the complete SUBJECT_IDs
set.seed(mc_seed)  # For reproducibility
sampled_ids <- sample(ids, N, replace = TRUE)

#Initialize a dataframe with required columns:
MC_dataset = complete_dat[-c(1:nrow(complete_dat)),]

for(i in 1:N)
{
  subdat = complete_dat %>% filter(SUBJECT_ID==sampled_ids[i])
  subdat$SUBJECT_ID = i
  
  MC_dataset = rbind(MC_dataset, subdat) 
}
#This takes care of the covariates; but we still need to generate the Y and update the FRI_Score column in MC_dataset.

#############################################
## Set "true" parameters in `true_params`: ##
#############################################

sim_varnames = c("gamma", "alpha", "alpha_1|2", "alpha_2|3", colnames(MC_dataset)[-c(1:2)]) #variable names

obj = load(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/05a_combined_presence_modelling/exchangeable,exchangeable/whole_data_based/coefs_pres_age9.RData"))
assign("coef_obj", get(obj))

true_params = coef_obj[coef_obj$Variables %in% sim_varnames,2] #these are the true values we need.
names(true_params) = sim_varnames


temp = update_MC_dataset(MC_dataset, mc_seed)

if(all(c(1,2) %in% unique(temp$FRI_Score)))
{
  MC_dataset = temp
  
  #Load files...initial estimates of presence, severity:
  
  path_pres = paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N",
                     N, "/", "presence" , "/", corstr_pres, ",", corstr_sev, "/age", age, "/",
                     "coef_pres_MC_", mc_seed, ".Rdata")
  
  pres_obj = load(path_pres) ; assign("sep_pres_coef", get(pres_obj))
  sep_pres_coef = sep_pres_coef$Estimates
  
  
  path_sev = paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N",
                    N, "/", "severity" , "/", corstr_pres, ",", corstr_sev, "/age", age, "/",
                    "coef_sev_MC_", mc_seed, ".Rdata")
  
  sev_obj = load(path_sev) ; assign("sep_sev_coef", get(sev_obj))
  sep_sev_coef = sep_sev_coef$Estimates
  
  
  #Fit the GEE:
  
  MC_output = fit_GEE_combined(dat=MC_dataset, corstr_pres = corstr_pres,
                               corstr_sev = corstr_sev, sep_pres_coef = sep_pres_coef,
                               sep_sev_coef = sep_sev_coef, kappa=0.25)
  
  coef_pres = MC_output$pres_coefs
  coef_sev = MC_output$sev_coefs
  
  gamma = MC_output$pres_coefs[1,2]
  
  rho_pres = 1; rho_sev = 1
  JK_corr_mat_pres = JK_corr_mat_sev = NA
  
  if(corstr_pres=="ar1" || corstr_pres=="exchangeable")
  {
    rho_pres = MC_output$rho_pres
    
  }else if(corstr_pres=="jackknifed")
  {
    JK_corr_mat_pres = MC_output$JK_corr_mat_pres
  }
  
  
  if(corstr_sev=="ar1" || corstr_sev=="exchangeable")
  {
    rho_sev = MC_output$rho_sev
    
  }else if(corstr_sev=="jackknifed")
  {
    JK_corr_mat_sev = MC_output$JK_corr_mat_sev
  }
  
  
  ##########################################################
  ## Get jackknifed estimates (for SD) and rhos (for CIs) ##
  ##########################################################
  
  coef_pres_JK = matrix( , nrow=length(true_params), ncol=N)
  rownames(coef_pres_JK) = names(true_params)
  
  coef_sev_JK = coef_pres_JK
  
  rho_pres_JK = rho_sev_JK = c()
  
  count = 1
  
  for(ID in unique(MC_dataset$SUBJECT_ID))
  {
    MC_dataset_JK = MC_dataset %>% filter(SUBJECT_ID != ID)
    
    #save only valid outputs
    if(all(c(1,2) %in% unique(MC_dataset_JK$FRI_Score)))
    {
      MC_output_JK = 
        fit_GEE_combined(dat=MC_dataset_JK, corstr_pres = corstr_pres,
                         sep_pres_coef = sep_pres_coef,
                         sep_sev_coef = sep_sev_coef,
                         kappa=0.25,
                         reuse = T, 
                         rho_pres = rho_pres, rho_sev = rho_sev,
                         JK_corr_mat_pres = JK_corr_mat_pres, 
                         JK_corr_mat_sev = JK_corr_mat_sev)
      
      #store JK outputs:
      
      coef_pres_JK[ -c(1,3,4), count] = MC_output_JK$pres_coefs$Estimates[-c(1,3,4)]
      coef_sev_JK[ -c(1,2), count] = MC_output_JK$sev_coefs$Estimates[-c(1,2)]
    }
    
    count = count+1
    
    print(ID)
  }
  
  #Get SD:
  JK_SD_pres = data.frame(apply(coef_pres_JK, MARGIN=1, FUN=compute_JK_SE)) 
  colnames(JK_SD_pres) = "SD"
  JK_SD_pres = JK_SD_pres[-c(1,3,4), , drop=F]
  
  JK_SD_sev = data.frame(apply(coef_sev_JK, MARGIN=1, FUN=compute_JK_SE)) 
  colnames(JK_SD_sev) = "SD"
  JK_SD_sev = JK_SD_sev[-c(1,2), , drop=F]
  
  #Get standardized estimates:
  std_coef_pres = coef_pres[,2]/JK_SD_pres
  std_coef_sev = coef_sev[,2]/JK_SD_sev
  
  mode = "combined"
  
  dir.create(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N))
  
  dir.create(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode))
  
  dir.create(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode, "/",
                    corstr_pres, ",", corstr_sev))
  
  dir.create(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode, "/",
                    corstr_pres, ",", corstr_sev, "/age", age))
  
  setwd(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode, "/",
               corstr_pres, ",", corstr_sev, "/age", age))
  
  save(gamma, file = paste0("gamma_MC_", mc_seed, ".Rdata") )
  
  save(coef_pres, file = paste0("coef_pres_MC_", mc_seed, ".Rdata") )
  save(coef_sev, file = paste0("coef_sev_MC_", mc_seed, ".Rdata") )
  
  if(corstr_pres=="ar1" || corstr_pres=="exchangeable"){
    save(rho_pres, file = paste0("rho_pres_MC_", mc_seed, ".Rdata") )
    save(rho_sev, file = paste0("rho_sev_MC_", mc_seed, ".Rdata") )
  }
  
  save(JK_SD_pres, file = paste0("pres_JK_SD_MC_", mc_seed, ".Rdata") )
  save(JK_SD_sev, file = paste0("sev_JK_SD_MC_", mc_seed, ".Rdata") )
  
  save(std_coef_pres, file = paste0("std_coef_pres_MC_", mc_seed, ".Rdata") )
  save(std_coef_sev, file = paste0("std_coef_sev_MC_", mc_seed, ".Rdata") )
  
}else{
  print("MC dataset lacks both 1/2 FRI Scores.")
}

end = Sys.time(); elapsed = end - start; elapsed