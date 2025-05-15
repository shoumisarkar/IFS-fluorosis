assumed_beta = -2

#Read the MC seed #
args = commandArgs(trailingOnly=TRUE)
mc_seed = as.numeric(args[1])
age = as.numeric(args[2])


start = Sys.time()

#Load libraries:
library(readxl)
library(dplyr)
library(MASS)
library(R.utils)

N = 30
B=10000
is_signif_pres = T
is_signif_sev = T
#mc_seed = 1 #change this as needed
#age = 23 #change this as needed

corstr_pres = "independence" #"jackknifed" #"jackknifed" # "exchangeable" #"independence"
corstr_sev = "independence" #"jackknifed" #"jackknifed" # "exchangeable"  #"independence" #"exchangeable"

exch_rho_pres = 0.3 #0.3 #0.6
exch_rho_sev = 0.8 #0.8 #0.6

print(paste0(corstr_pres, ",", corstr_sev, ",", exch_rho_pres, "," , exch_rho_sev))

setwd("/blue/somnath.datta/shoumisarkar/Fluorosis")
source("Codes/functions.R")


print(paste0("Running with MC seed ", mc_seed))

#Create directories

mode = "combined"

dir.create("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy")
dir.create(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy/N", N, " - beta=", assumed_beta, "/"))
dir.create(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy/N", N, " - beta=", assumed_beta, "/", mode))
dir.create(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy/N", N, " - beta=", assumed_beta,  "/", mode, "/",
                  corstr_pres, ",", corstr_sev))
dir.create(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy/N", N, " - beta=", assumed_beta,  "/", mode, "/",
                  corstr_pres, ",", corstr_sev, "/age", age))
path = paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy/N", N, " - beta=", assumed_beta,  "/", mode, "/",
              corstr_pres, ",", corstr_sev, "/age", age, "/coef_pres_MC_", mc_seed, ".Rdata")

corresp_sev_file = paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy/N", N, " - beta=", assumed_beta, "/", "severity", "/",
                          corstr_pres, ",", corstr_sev, "/age", age,"/bs_Avghomeppm_MC_", mc_seed, ".Rdata")

if(file.exists(path))
{
  print("Output already generated")
}else if(!file.exists(corresp_sev_file)){
  print("Matching severity outputs do not exist for this mc_seed")
}else{
  
  #Load the data for the chosen age group:
  dat <- read_excel(paste0("Results/02_preprocess_data/preprocessed_IFS_data_age", age, ".xlsx"))
  
  #Keep selected columns:
  dat <- dat %>% dplyr::select(SUBJECT_ID, FRI_Score, #dental_age, 
                               Avg_homeppm, 
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
  
  true_params[5] = assumed_beta #Set Avg_homeppm to 0
  
  temp = update_MC_dataset(MC_dataset, mc_seed)
  
  #table(temp$FRI_Score)
  
  
  if(length(unique(temp$FRI_Score))!=4)
  {
    print("length(unique(temp$FRI_Score))!=4")
    print(unique(temp$FRI_Score)) 
    
    if(all(c(0,2,3) %in% unique(temp$FRI_Score)))
    {
      copy = temp$FRI_Score
      
      copy[copy==2] = 1
      copy[copy==3] = 2
      
      temp$FRI_Score = copy #in case the only categories present are 0,2,3, rename them to 0,1,2
      
    }
  }
  
  
  
  if(all(c(0,1,2) %in% unique(temp$FRI_Score))){
    
    #########################################
    
    #Load Presence coefficients:
    
    setwd(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy/N", N, " - beta=", assumed_beta,  "/", "presence", "/",
                 corstr_pres, ",", corstr_sev, "/age", age))
    pres_file = paste0("coef_pres_MC_", mc_seed, ".Rdata")    
    pres_obj = load(pres_file)
    assign("pres_obj", get(pres_obj))
    sep_pres_coef = pres_obj$Estimates
    
    #Load Severity coefficients:
    
    setwd(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy/N", N, " - beta=", assumed_beta,  "/", "severity", "/",
                 corstr_pres, ",", corstr_sev, "/age", age))
    sev_file = paste0("coef_sev_MC_", mc_seed, ".Rdata")    
    sev_obj = load(sev_file)
    assign("sev_obj", get(sev_obj))
    sep_sev_coef = sev_obj$Estimates
    
    #########################################
    
    MC_dataset = temp 
    MC_output = fit_GEE_combined(dat=MC_dataset, corstr_pres = corstr_pres,
                                 corstr_sev = corstr_sev,
                                 sep_pres_coef = sep_pres_coef,
                                 sep_sev_coef = sep_sev_coef,
                                 kappa=0.25, maxIter = 10)

    coef_pres = MC_output$pres_coefs
    coef_sev = MC_output$sev_coefs
    
    # rho_pres = 1
    # 
    # JK_corr_mat_pres = NA
    # 
    # if(corstr_pres=="ar1" || corstr_pres=="exchangeable"){  
    #   rho_pres = MC_output$rho_pres
    # }else if(corstr_pres=="jackknifed")
    # {
    #   JK_corr_mat_pres = MC_output$JK_corr_mat_pres
    # }
    
    ##############################
    ## Bootstrapping to get CIs ##
    ##############################
    
    coef_pres_BS = matrix( , nrow=length(true_params), ncol=B)
    rownames(coef_pres_BS) = names(true_params)
    
    coef_sev_BS = coef_pres_BS
    
    #rho_pres_BS = c()
    
    for(b in 1:B)
    {
      #Get the b-th bootstrapped dataset
      set.seed(b)
      bs_IDs = sample(unique(temp$SUBJECT_ID), size =N, replace = T)
      
      bs_dataset = complete_dat[-c(1:nrow(complete_dat)),]
      count=1
      
      for(i in bs_IDs)
      {
        subdat = temp %>% filter(SUBJECT_ID==i)
        subdat$SUBJECT_ID = count; count=count+1
        bs_dataset = rbind(bs_dataset, subdat) 
      }
      
      if(length(unique(bs_dataset$FRI_Score))>=2)
      {
        # bs_output = fit_GEE_presence(dat=bs_dataset, corstr_pres = corstr_pres, kappa=0.25, maxIter = 10)
        # coef_pres_BS[-c(1,3,4),b] = bs_output$coefs$Estimates
        
        # tryCatch({
        #   bs_output <- fit_GEE_combined(dat = bs_dataset, corstr_pres = corstr_pres,
        #                                 corstr_sev = corstr_sev,
        #                                 sep_pres_coef = sep_pres_coef, 
        #                                 sep_sev_coef = sep_sev_coef,
        #                                 kappa = 0.25, maxIter = 10)
        #   coef_pres_BS[, b] <- bs_output$pres_coefs$Estimates
        #   coef_sev_BS[, b] <- bs_output$sev_coefs$Estimates
        #   
        # }, error = function(e) {
        #   message(sprintf("Bootstrap iteration %d failed (combined model): %s", b, e$message))
        #   # Optionally store NA or another placeholder:
        #   coef_pres_BS[, b] <- NA
        #   coef_sev_BS[, b] <- NA
        # })
        
        tryCatch({
          # Fit combined model with timeout protection
          bs_output <- withTimeout({
            fit_GEE_combined(dat = bs_dataset, 
                             corstr_pres = corstr_pres,
                             corstr_sev = corstr_sev,
                             sep_pres_coef = sep_pres_coef, 
                             sep_sev_coef = sep_sev_coef,
                             kappa = 0.25, 
                             maxIter = 10)
          }, timeout = 16, onTimeout = "error")  # 5 seconds timeout
          
          # If successful, store the estimates
          coef_pres_BS[, b] <- bs_output$pres_coefs$Estimates
          coef_sev_BS[, b] <- bs_output$sev_coefs$Estimates
        }, error = function(e) {
          if (inherits(e, "TimeoutException")) {
            message(sprintf("Bootstrap iteration %d failed due to TIMEOUT (combined model).", b))
          } else {
            message(sprintf("Bootstrap iteration %d failed (combined model): %s", b, e$message))
          }
          # Store NA on error
          coef_pres_BS[, b] <- NA
          coef_sev_BS[, b] <- NA
        })
        
      }
      
      print(b)
    }
    
    bs_pres_Avghomeppm = coef_pres_BS[5,]
    bs_sev_Avghomeppm = coef_sev_BS[5,]
    
    lci_pres = quantile(bs_pres_Avghomeppm, probs = 0.025, na.rm = T)
    uci_pres = quantile(bs_pres_Avghomeppm, probs = 0.975, na.rm = T)
    
    lci_sev = quantile(bs_sev_Avghomeppm, probs = 0.025, na.rm = T)
    uci_sev = quantile(bs_sev_Avghomeppm, probs = 0.975, na.rm = T)
    
    if(lci_pres<=0 && 0<=uci_pres)
    {
      is_signif_pres = F
    }
    
    if(lci_sev<=0 && 0<=uci_sev)
    {
      is_signif_sev = F
    }
    
    setwd(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy/N", N, " - beta=", assumed_beta,  "/", mode, "/",
                 corstr_pres, ",", corstr_sev, "/age", age))
    
    save(coef_pres, file = paste0("coef_pres_MC_", mc_seed, ".Rdata") )
    save(coef_sev, file = paste0("coef_sev_MC_", mc_seed, ".Rdata") )
    
    save(is_signif_pres, file = paste0("is_signif_pres_MC_", mc_seed, ".Rdata") )
    save(is_signif_sev, file = paste0("is_signif_sev_MC_", mc_seed, ".Rdata") )
    
    save(bs_pres_Avghomeppm, file = paste0("bs_pres_Avghomeppm_MC_", mc_seed, ".Rdata") )
    save(bs_sev_Avghomeppm, file = paste0("bs_sev_Avghomeppm_MC_", mc_seed, ".Rdata") )
    
  }else{
    print("MC dataset lacks 0/1/2 FRI Scores.")
  }
  
}


end = Sys.time()

elapsed = end - start; elapsed
