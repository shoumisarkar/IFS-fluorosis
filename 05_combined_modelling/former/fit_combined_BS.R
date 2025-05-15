
#####################
###### Read b #######
#####################

#Read the bootstrap sample #
args = commandArgs(trailingOnly=TRUE)

b = as.numeric(args[1])
age = as.numeric(args[2])
corstr = args[3]

# Next: need to cluster-based bootstrap this process B=1000 times.
start = Sys.time()

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
if(inSLURM) {
  # If in SLURM environment
  setwd("/blue/somnath.datta/shoumisarkar/Fluorosis/")
} else {
  # If not in SLURM environment
  setwd("W:/somnath.datta/shoumisarkar/Fluorosis/")
}

source("Codes/functions.R")

file_pres = paste0(ifelse(inSLURM, "/blue/", "W:/"), 
                   "somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/", 
                   corstr, "/bootstrapping/age", age, "/BS_std_estimators_pres_age", 
                   age, "_b", b, ".RData") 
  
file_sev = paste0(ifelse(inSLURM, "/blue/", "W:/"),
                  "somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/", 
                  corstr, "/bootstrapping/age", age, "/BS_std_estimators_sev_age",
                  age, "_b", b, ".RData")

if(file.exists(file_pres) & file.exists(file_sev))
{
  filename = paste0("Results/02_preprocess_data/preprocessed_IFS_data_age", age, ".xlsx")
  IFS_dat <- read_xlsx(filename)
  
  ##### Bootstrap:
  
  #BS_std_estimators = matrix( , nrow = (1+1+2+(ncol(IFS_dat)-2)),
  #                                  ncol = B )
  
  ids = unique(IFS_dat$SUBJECT_ID)
  
  set.seed(b)
  
  selected_ids = sample(x = ids, size = length(ids), replace = T)
  
  IFS_dat_BS = IFS_dat %>% filter(SUBJECT_ID == selected_ids[1])
  IFS_dat_BS = IFS_dat_BS[-c(1:nrow(IFS_dat_BS)),] #make it an empty dataframe
  
  new_sid = 1
  
  for(t in 1:length(selected_ids))
  {
    temp = IFS_dat %>% filter(SUBJECT_ID == selected_ids[t])
    temp$SUBJECT_ID = new_sid #relabel the subject ID, as we treat repeated clusters as new clusters
    IFS_dat_BS = rbind(IFS_dat_BS, temp)
    
    new_sid = new_sid + 1
  }
  
  obj_pres = load(file_pres)
  assign("coefs_pres", get(obj_pres))
  
  obj_sev = load(file_sev)
  assign("coefs_sev", get(obj_sev))

  #gamma_b = mean(coefs_sev[-c(1:2)]/coefs_pres[-c(1)])
  
  fit = fit_GEE(dat=IFS_dat_BS, corstr_pres=corstr, corstr_sev=corstr,
                init_pres_coef = coefs_pres, init_sev_coef = coefs_sev); coef = fit$coefs
  
  coef_pres = coef
  coef_sev = coef; coef_sev[-c(1:4),2] = coef_sev[1,2] * coef_sev[-c(1:4),2]
  
  rho_pres = fit$rho_pres
  rho_sev = fit$rho_sev
  
  filename = paste0("Results/05a_combined_presence_modelling/", corstr, "/bootstrapping/age", age, 
                    "/BS_rho_pres_age", age, "_b", b,".RData")
  save(rho_pres, file = filename)
  
  
  filename = paste0("Results/05b_combined_severity_modelling/", corstr, "/bootstrapping/age", age, 
                    "/BS_rho_sev_age", age, "_b", b,".RData")
  save(rho_sev, file = filename)
  
  #######################################
  ###### Get jackknifed estimates #######
  #######################################
  
  obj1 = load(paste0(ifelse(inSLURM, "/blue/", "W:/"), 
                "somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/", 
                corstr, "/bootstrapping/age", age, "/BS_coefs_JK_pres_age", 
                age, "_b", b, ".RData"))
  assign("BS_coefs_JK_pres", get(obj1))
  
  obj2 = load(paste0(ifelse(inSLURM, "/blue/", "W:/"), 
                     "somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/", 
                     corstr, "/bootstrapping/age", age, "/BS_coefs_JK_sev_age", 
                     age, "_b", b, ".RData"))
  assign("BS_coefs_JK_sev", get(obj2))
  
  coefs_JK_pres = matrix( , nrow=(1+1+2+(ncol(IFS_dat_BS)-2)),
                          ncol= length(unique(IFS_dat_BS$SUBJECT_ID)) )
  
  coefs_JK_sev = matrix( , nrow=(1+1+2+(ncol(IFS_dat_BS)-2)),
                         ncol= length(unique(IFS_dat_BS$SUBJECT_ID)) )
  
  rho_pres_JK = c()
  rho_sev_JK = c()
  
  for(id in unique(IFS_dat_BS$SUBJECT_ID))
  {
    print(paste0("JK leave ",id, " out"))
    
    JK_subdat = IFS_dat_BS %>% filter(IFS_dat_BS$SUBJECT_ID != id)  #leave cluster i out
    
    if(length(unique(JK_subdat$FRI_Score))<4)
    {
      next
    }
    
    mod_jk = fit_GEE(JK_subdat, corstr_pres=corstr, corstr_sev=corstr, 
                     init_pres_coef = as.vector(BS_coefs_JK_pres[,b]),
                     init_sev_coef = as.vector(BS_coefs_JK_sev[,b]))
    
    ind = which(unique(IFS_dat_BS$SUBJECT_ID) == id)
    #coefs_JK[,ind] = (mod_jk$coefs)[,2]
    
    coefs_JK_pres[,ind] = (mod_jk$coefs)[,2]
    
    coefs_JK_sev[,ind] = (mod_jk$coefs)[,2]
    coefs_JK_sev[-c(1:4),ind] = coefs_JK_sev[1,ind] * coefs_JK_sev[-c(1:4),ind] 
    
    rho_pres_JK = c(rho_pres_JK, mod_jk$rho_pres)
    rho_sev_JK = c(rho_sev_JK, mod_jk$rho_sev)
  }
  
  rownames(coefs_JK_pres) =
    rownames(coefs_JK_sev) =
    rownames(mod_jk$coefs)
  
  ##################################################
  ###### Get SD of the jackknifed estimators #######
  ##################################################
  
  JK_SD_pres = data.frame(apply(coefs_JK_pres, MARGIN=1, FUN= compute_JK_SE))
  JK_SD_sev = data.frame(apply(coefs_JK_sev, MARGIN=1, FUN= compute_JK_SE))
  #colnames(JK_SD) = "SD"
  
  ##################################################
  #########  Standardize the estimators:   #########
  ##################################################
  
  std_coef_pres = data.frame(Variables = rownames(mod_jk$coefs),
                             StdCoef = coef_pres[,2]/JK_SD_pres)
  
  std_coef_sev = data.frame(Variables = rownames(mod_jk$coefs),
                            StdCoef = coef_sev[,2]/JK_SD_sev)
  
  #BS_std_estimators[,b] =  
  BS_std_estimators_pres = na.omit(as.numeric(unlist(std_coef_pres)))
  BS_std_estimators_sev = na.omit(as.numeric(unlist(std_coef_sev)))
  
  filename = paste0("Results/05a_combined_presence_modelling/", corstr, "/bootstrapping/age", age, 
                    "/BS_std_estimators_pres_age", age, "_b", b,".RData")
  save(BS_std_estimators_pres, file = filename)
  
  
  filename = paste0("Results/05b_combined_severity_modelling/", corstr, "/bootstrapping/age", age, 
                    "/BS_std_estimators_sev_age", age, "_b", b,".RData")
  save(BS_std_estimators_sev, file = filename)
  
  
  filename = paste0("Results/05a_combined_presence_modelling/", corstr, "/bootstrapping/age", age, 
                    "/BS_rho_pres_JK_age", age, "_b", b,".RData")
  save(rho_pres_JK, file = filename)
  
  
  filename = paste0("Results/05b_combined_severity_modelling/", corstr, "/bootstrapping/age", age, 
                    "/BS_rho_sev_JK_age", age, "_b", b,".RData")
  save(rho_sev_JK, file = filename)
  
  print(b)
  
  end = Sys.time()
  
  elapsed = end - start; elapsed
}else{

  if(!file.exists(file_pres))
    {
      print("Presence file missing")
    }
  if(!file.exists(file_sev))
    {
      print("Severity file missing")
    }
}

