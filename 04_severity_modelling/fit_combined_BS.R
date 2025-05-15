#####################
###### Read b #######
#####################

#Read the bootstrap sample #
args = commandArgs(trailingOnly=TRUE)

b = as.numeric(args[1])
age = as.numeric(args[2])
corstr_pres = args[3]
corstr_sev = args[4]

# Next: need to cluster-based bootstrap this process B=100 times.
start = Sys.time()

######################
### Load packages ####
######################

library(writexl)
library(readxl)

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

filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev,
                  "/bootstrapping/age", age, 
                  "/std_coefs_pres_age", age, "_b", b,".RData")

if(!file.exists(filename))
{
  source("Codes/functions.R")
  
  filename = paste0("Results/02_preprocess_data/preprocessed_IFS_data_age", age, ".xlsx")
  IFS_dat <- read_xlsx(filename)
  
  ###########################################################
  ####### Read in estimates from separate models ############
  ###########################################################
  
  path_pres = paste0(ifelse(inSLURM, "/blue/", "W:/"), 
                     "somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/", 
                     corstr_pres, "/bootstrapping/age",age)
  
  path_sev = paste0(ifelse(inSLURM, "/blue/", "W:/"), 
                    "somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/", 
                    corstr_sev, "/bootstrapping/age", age)
  
  
  obj_pres = load(paste0(path_pres,"/coefs_pres_age", age, "_b", b, ".RData"))
  assign("coef_pres", get(obj_pres))
  
  obj_sev = load(paste0(path_sev,"/coefs_sev_age", age, "_b", b, ".RData"))
  assign("coef_sev", get(obj_sev))
  
  ######################################################################
  ####### Read in jackknifed estimates from separate models ############
  ######################################################################
  
  obj_pres_jk = load(paste0(path_pres,"/coefs_JK_pres_age", age, "_b", b, ".RData"))
  assign("coef_pres_jk", get(obj_pres_jk))
  
  obj_sev_jk = load(paste0(path_sev,"/coefs_JK_sev_age", age, "_b", b, ".RData"))
  assign("coef_sev_jk", get(obj_sev_jk))
  
  
  ##### Bootstrap to get the current dataset:
  
  all_IDs = unique(IFS_dat$SUBJECT_ID)
  
  set.seed(b)
  
  selected_ids = sample(x = all_IDs, size = length(all_IDs), replace = T) #select clusters with replacement
  
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
  
  if(length(unique(IFS_dat_BS$FRI_Score))==4 || any(is.na(c(coef_pres$Estimates,
                                                            coef_sev$Estimates))))
  {
    
    fit = fit_GEE_combined(dat=IFS_dat_BS, corstr_pres = corstr_pres, corstr_sev = corstr_sev,
                           init_pres_coef = coef_pres$Estimates, init_sev_coef = coef_sev$Estimates,
                           maxIter = 100
    )
    
    #Save presence results...
    combined_coef_pres = fit$pres_coefs
    filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev,
                      "/bootstrapping/age", age,
                      "/coefs_pres_b", b, ".RData")
    save(combined_coef_pres, file = filename) 
    
    #Iterations...
    iter=fit$iter
    filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev,
                      "/bootstrapping/age", age, 
                      "/iter_pres_b", b, ".RData")
    save(iter, file = filename) 
    
    #Save severity results...
    combined_coef_sev = fit$sev_coefs
    filename = paste0("Results/05b_combined_severity_modelling/", corstr_pres, ",", corstr_sev,
                      "/bootstrapping/age", age,
                      "/coefs_sev_b", b, ".RData")
    save(combined_coef_sev, file = filename) 
    
    #Save the estimates of rho, if applicable
    
    #Presence
    if(corstr_pres=="ar1" || corstr_pres=="exchangeable")
    {
      rho_pres = fit$rho_pres
      filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev,
                        "/bootstrapping/age", age, 
                        "/rho_pres_b", b, ".RData")
      save(rho_pres, file = filename) 
    }
    
    #Severity
    if(corstr_sev=="ar1" || corstr_sev=="exchangeable")
    {
      rho_sev = fit$rho_sev
      filename = paste0("Results/05b_combined_severity_modelling/", corstr_pres, ",", corstr_sev,
                        "/bootstrapping/age", age,
                        "/rho_sev_b", b, ".RData")
      save(rho_sev, file = filename) 
    }
    
    #######################################
    ###### Get jackknifed estimates #######
    #######################################
    
    all_IDs = unique(IFS_dat_BS$SUBJECT_ID)
    
    combined_coefs_JK_pres = matrix( , nrow=(1+1+2+(ncol(IFS_dat)-2)),
                                     ncol= length(all_IDs) )
    
    combined_coefs_JK_sev = matrix( , nrow=(1+1+2+(ncol(IFS_dat_BS)-2)),
                                    ncol= length(all_IDs) )
    
    rho_pres_JK = c(); rho_sev_JK = c()
    iter_JK = c()
    
    for(id in all_IDs)
    {
      print(paste0("JK leave ",id, " out"))
      
      JK_subdat = IFS_dat_BS %>% filter(IFS_dat_BS$SUBJECT_ID != id)  #leave cluster i out
      
      if(length(unique(JK_subdat$FRI_Score))<4 || any(is.na(c(unlist(coef_pres_jk[,id]), 
                                                              unlist(coef_sev_jk[,id])))))
      {
        next
      }
      
      #for init_pres_coef, and init_sev_coef, can use the JK coefs saved from the separate analyses.
      
      temp_fit = fit_GEE_combined(dat = JK_subdat, corstr_pres = corstr_pres, 
                                  corstr_sev = corstr_sev, 
                                  init_pres_coef = unlist(coef_pres_jk[,id]), 
                                  init_sev_coef = unlist(coef_sev_jk[,id]),
                                  maxIter = 100
      )
      
      ind = which(all_IDs == id)
      
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
                        "/bootstrapping/age", age, "/coefs_JK_pres_b", b, ".RData")
      save(combined_coefs_JK_pres, file = filename)
      
      filename = paste0("Results/05b_combined_severity_modelling/", corstr_pres, ",", corstr_sev, 
                        "/bootstrapping/age", age, "/coefs_JK_sev_b", b, ".RData")
      save(combined_coefs_JK_sev, file = filename)
    }
    
    #Save rho, if applicable...
    
    #Presence:
    if(corstr_pres=="ar1" || corstr_pres=="exchangeable")
    { 
      filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev, 
                        "/bootstrapping/age", age, "/rho_pres_JK_age", age, "_b",
                        b,".RData")
      save(rho_pres_JK, file = filename)
    }
    
    #Severity:
    if(corstr_sev=="ar1" || corstr_sev=="exchangeable")
    { 
      filename = paste0("Results/05b_combined_severity_modelling/", corstr_pres, ",", corstr_sev, 
                        "/bootstrapping/age", age, "/rho_sev_JK_age", age, "_b",
                        b,".RData")
      save(rho_sev_JK, file = filename)
    }
    
    ##################################################
    ###### Get SD of the jackknifed estimators #######
    ##################################################
    
    JK_SD_pres = data.frame(apply(combined_coefs_JK_pres, MARGIN=1, FUN= compute_JK_SE))
    colnames(JK_SD_pres) = "SD"
    
    JK_SD_sev = data.frame(apply(combined_coefs_JK_sev, MARGIN=1, FUN= compute_JK_SE))
    colnames(JK_SD_sev) = "SD"
    
    filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev,
                      "/bootstrapping/age", age, 
                      "/JK_SD_pres_age", age, "_b", b,".RData")
    save(JK_SD_pres, file = filename)
    
    filename = paste0("Results/05b_combined_severity_modelling/", corstr_pres, ",", corstr_sev,
                      "/bootstrapping/age", age, 
                      "/JK_SD_sev_age", age, "_b", b,".RData")
    save(JK_SD_sev, file = filename)
    
    ##################################################
    #########  Standardize the estimators:   #########
    ##################################################
    
    std_coef_pres = combined_coef_pres[,2]/JK_SD_pres
    colnames(std_coef_pres) = "Standardized Estimate"
    
    std_coef_sev = combined_coef_sev[,2]/JK_SD_sev
    colnames(std_coef_sev) = "Standardized Estimate"
    
    filename = paste0("Results/05a_combined_presence_modelling/", corstr_pres, ",", corstr_sev,
                      "/bootstrapping/age", age, 
                      "/std_coefs_pres_age", age, "_b", b,".RData")
    save(std_coef_pres, file = filename)
    
    filename = paste0("Results/05b_combined_severity_modelling/", corstr_pres, ",", corstr_sev,
                      "/bootstrapping/age", age, 
                      "/std_coefs_sev_age", age, "_b", b,".RData")
    save(std_coef_sev, file = filename)
    
    print(b)
    
    end = Sys.time()
    
  }else{
    print("length(unique(IFS_dat_BS$FRI_Score))<4, or estimates from separate models contain NA.")  
  }
  
  elapsed = end - start; elapsed
  
}else{
  print("Output for this case has already exists.")  
}
