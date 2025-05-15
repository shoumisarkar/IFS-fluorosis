#####################
###### Read b #######
#####################

#Read the bootstrap sample #
args = commandArgs(trailingOnly=TRUE)

b = as.numeric(args[1])
age = as.numeric(args[2])
corstr_sev = args[3]

# Next: need to cluster-based bootstrap this process B=1000 times.
start = Sys.time()

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

filename = paste0("Results/04_severity_modelling/", corstr_sev, "/bootstrapping/age", age, 
                  "/std_coefs_sev_age", age, "_b", b,".RData")

if(!file.exists(filename))
{
  source("Codes/functions.R")
  
  filename = paste0("Results/02_preprocess_data/preprocessed_IFS_data_age", age, ".xlsx")
  IFS_dat <- read_xlsx(filename)
  
  ##### Bootstrap:
  
  all_IDs = unique(IFS_dat$SUBJECT_ID)
  
  set.seed(b)
  
  selected_ids = sample(x = all_IDs, size = length(all_IDs), replace = T)
  
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
  
  #####################################
  ###### Get severity estimates #######
  #####################################
  
  if(length(unique(IFS_dat_BS$FRI_Score))==4)
  {
    
    fit = fit_GEE_severity(IFS_dat_BS, corstr_sev=corstr_sev, maxIter = 100)
    
    coef_sev = fit$coefs
    
    filename = paste0("Results/04_severity_modelling/", corstr_sev, "/bootstrapping/age", age, 
                      "/coefs_sev_age", age, "_b", b,".RData")
    save(coef_sev, file = filename)
    
    iter = fit$iter
    filename = paste0("Results/04_severity_modelling/", corstr_sev, "/bootstrapping/age", age, 
                      "/iter_sev_age", age, "_b", b,".RData")
    save(iter, file = filename)
    
    
    rho_sev = fit$rho_sev
    filename = paste0("Results/04_severity_modelling/", corstr_sev, "/bootstrapping/age", age, 
                      "/rho_sev_age", age, "_b", b,".RData")
    save(rho_sev, file = filename)
    
    #######################################
    ###### Get jackknifed estimates #######
    #######################################
    
    all_IDs = unique(IFS_dat_BS$SUBJECT_ID)
    
    coefs_JK_sev = matrix( , nrow=(0+0+2+(ncol(IFS_dat_BS)-2)),
                           ncol= length(all_IDs) )
    
    rho_sev_JK = c()
    iter_JK = c()
    
    for(id in all_IDs)
    {
      print(paste0("JK leave ",id, " out"))
      
      JK_subdat = IFS_dat_BS %>% filter(IFS_dat_BS$SUBJECT_ID != id)  #leave cluster i out
      
      if(length(unique(JK_subdat$FRI_Score))<4)
      {
        next
      }
      
      mod_jk = fit_GEE_severity(JK_subdat, corstr_sev=corstr_sev, maxIter = 100)
      temp_coefs = mod_jk$coefs
      ind = which(all_IDs == id)
      coefs_JK_sev[,ind] = temp_coefs[,2]
      
      if(corstr_sev!="jackknifed")
      {
        rho_sev_JK = c(rho_sev_JK, mod_jk$rho_sev)
      }
      
      iter_JK = c(iter_JK, mod_jk$iter)
      
      print(paste0("Jackknifing... ", ind, "/", length(all_IDs)))
      
      rownames(coefs_JK_sev) = temp_coefs[,1]
      
      filename = paste0("Results/04_severity_modelling/", corstr_sev, 
                        "/bootstrapping/age", age, "/coefs_JK_sev_age", age, "_b", 
                        b,".RData")
      save(coefs_JK_sev, file = filename)
    }
    
    #For ar1/exchangeable, save the jackknifed rho's...
    if(corstr_sev!="jackknifed")
    { 
      filename = paste0("Results/04_severity_modelling/", corstr_sev, 
                        "/bootstrapping/age", age, "/rho_sev_JK_age", age, "_b",
                        b,".RData")
      save(rho_sev_JK, file = filename)
      
    }
    
    ##################################################
    ###### Get SD of the jackknifed estimators #######
    ##################################################
    
    JK_SD_sev = data.frame(apply(coefs_JK_sev, MARGIN=1, FUN= compute_JK_SE))
    colnames(JK_SD_sev) = "SD"
    
    filename = paste0("Results/04_severity_modelling/", corstr_sev, "/bootstrapping/age", age, 
                      "/JK_SD_sev_age", age, "_b", b,".RData")
    save(JK_SD_sev, file = filename)
    
    ##################################################
    #########  Standardize the estimators:   #########
    ##################################################
    
    std_coef_sev = coef_sev[,2]/JK_SD_sev
    colnames(std_coef_sev) = "Standardized Estimate"
    
    filename = paste0("Results/04_severity_modelling/", corstr_sev, "/bootstrapping/age", age, 
                      "/std_coefs_sev_age", age, "_b", b,".RData")
    save(std_coef_sev, file = filename)
    
    print(b)
    
  }else{
    print("length(unique(IFS_dat_BS$FRI_Score))<4")  
  }
  
  end = Sys.time()

}else{
  print("Output for this case has already exists.")
}

elapsed = end - start; elapsed