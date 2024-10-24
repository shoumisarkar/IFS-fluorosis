#####################
###### Read b #######
#####################

#Read the bootstrap sample #
args = commandArgs(trailingOnly=TRUE)

b = as.numeric(args[1])
age = as.numeric(args[2])
corstr_pres = args[3]

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

filename = paste0("Results/03_presence_modelling/", corstr_pres, "/bootstrapping/age", age, 
                  "/std_coefs_pres_age", age, "_b", b,".RData")


if(!file.exists(filename)) #run only if the referenced output file does not exist yet 
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
  ###### Get presence estimates #######
  #####################################
  
  if(length(unique(IFS_dat_BS$FRI_Score))==4)
  {
    
    fit = fit_GEE_presence(IFS_dat_BS, corstr_pres=corstr_pres)
    
    coef_pres = fit$coefs
    
    filename = paste0("Results/03_presence_modelling/", corstr_pres, "/bootstrapping/age", age, 
                      "/coefs_pres_age", age, "_b", b,".RData")
    save(coef_pres, file = filename)
    
    iter = fit$iter
    filename = paste0("Results/03_presence_modelling/", corstr_pres, "/bootstrapping/age", age, 
                      "/iter_pres_age", age, "_b", b,".RData")
    save(iter, file = filename)
    
    
    rho_pres = fit$rho_pres
    filename = paste0("Results/03_presence_modelling/", corstr_pres, "/bootstrapping/age", age, 
                      "/rho_pres_age", age, "_b", b,".RData")
    save(rho_pres, file = filename)
    
    #######################################
    ###### Get jackknifed estimates #######
    #######################################
    
    all_IDs = unique(IFS_dat_BS$SUBJECT_ID)
    
    coefs_JK_pres = matrix( , nrow=(0+1+0+(ncol(IFS_dat_BS)-2)),
                            ncol= length(all_IDs) )
    
    rho_pres_JK = c()
    iter_JK = c()
    
    for(id in all_IDs)
    {
      print(paste0("JK leave ",id, " out"))
      
      JK_subdat = IFS_dat_BS %>% filter(IFS_dat_BS$SUBJECT_ID != id)  #leave cluster i out
      
      if(length(unique(JK_subdat$FRI_Score))<4)
      {
        next
      }
      
      mod_jk = fit_GEE_presence(JK_subdat, corstr_pres=corstr_pres)
      temp_coefs = mod_jk$coefs
      ind = which(all_IDs == id)
      coefs_JK_pres[,ind] = temp_coefs[,2]
      
      if(corstr_pres!="jackknifed")
      {
        rho_pres_JK = c(rho_pres_JK, mod_jk$rho_pres)
      }
      
      iter_JK = c(iter_JK, mod_jk$iter)
      
      print(paste0("Jackknifing... ", ind, "/", length(all_IDs)))
      
      rownames(coefs_JK_pres) = temp_coefs[,1]
      
      filename = paste0("Results/03_presence_modelling/", corstr_pres, 
                        "/bootstrapping/age", age, "/coefs_JK_pres_age", age, "_b", 
                        b,".RData")
      save(coefs_JK_pres, file = filename)
    }
    
    #For ar1/exchangeable, save the jackknifed rho's...
    if(corstr_pres!="jackknifed")
    { 
      filename = paste0("Results/03_presence_modelling/", corstr_pres, 
                        "/bootstrapping/age", age, "/rho_pres_JK_age", age, "_b",
                        b,".RData")
      save(rho_pres_JK, file = filename)
      
    }
    
    ##################################################
    ###### Get SD of the jackknifed estimators #######
    ##################################################
    
    JK_SD_pres = data.frame(apply(coefs_JK_pres, MARGIN=1, FUN= compute_JK_SE))
    colnames(JK_SD_pres) = "SD"
    
    filename = paste0("Results/03_presence_modelling/", corstr_pres, "/bootstrapping/age", age, 
                      "/JK_SD_pres_age", age, "_b", b,".RData")
    save(JK_SD_pres, file = filename)
    
    ##################################################
    #########  Standardize the estimators:   #########
    ##################################################
    
    std_coef_pres = coef_pres[,2]/JK_SD_pres
    colnames(std_coef_pres) = "Standardized Estimate"
    
    filename = paste0("Results/03_presence_modelling/", corstr_pres, "/bootstrapping/age", age, 
                      "/std_coefs_pres_age", age, "_b", b,".RData")
    save(std_coef_pres, file = filename)
    
    print(b)
    
    end = Sys.time()
    
  }else{
    print("length(unique(IFS_dat_BS$FRI_Score))<4")  
  }

}else{
  print("Output for this case has already exists.")
}

elapsed = end - start; elapsed