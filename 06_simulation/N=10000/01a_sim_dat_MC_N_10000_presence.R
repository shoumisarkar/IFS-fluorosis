
#Read the MC seed #
args = commandArgs(trailingOnly=TRUE)
mc_seed = as.numeric(args[1])
age = as.numeric(args[2])


start = Sys.time()

#Load libraries:
library(readxl)
library(dplyr)
library(MASS)

N = 10000
#mc_seed = 1 #change this as needed
#age = 9 #change this as needed

corstr_pres = "independence"
corstr_sev = "independence" 

exch_rho_pres = 0.3 #0.6
exch_rho_sev = 0.8 #0.6

print(paste0(corstr_pres, ",", corstr_sev, ",", exch_rho_pres, "," , exch_rho_sev))

print(paste0(corstr_pres, ",", corstr_sev, ",", exch_rho_pres, "," , exch_rho_sev))

setwd("/blue/somnath.datta/shoumisarkar/Fluorosis")
source("Codes/functions.R")


print(paste0("Running with MC seed ", mc_seed))

#Create directories

mode = "presence"
dir.create(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N))

dir.create(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode))

dir.create(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode, "/",
                  corstr_pres, ",", corstr_sev))

dir.create(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode, "/",
                  corstr_pres, ",", corstr_sev, "/age", age))


path = paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode, "/",
                    corstr_pres, ",", corstr_sev, "/age", age, "/coef_pres_MC_", mc_seed, ".Rdata")


if(file.exists(path))
{
  print("Output already generated")
}else{
  
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
  
  
  update_MC_dataset2 = function(dataset=MC_dataset, seed=mc_seed)
  {
    #Update the Y in dataset:
    #####################################################
    ##  Set up Sigma, the cluster correlation matrix:  ##
    #####################################################
    
    #Initialize:
    Sigma_P = diag(16) 
    
    # if(corstr_pres=="independence")
    # { 
    #   #With independence cluster correlation structure
    #   Sigma_P = diag(16)
    # }else if(corstr_pres=="exchangeable")
    # { 
      #With exchangeable cluster correlation structure
      Sigma_P = matrix(exch_rho_pres, 16, 16); diag(Sigma_P) = 1
    # }
    # if(corstr_pres=="jackknifed")
    # {
    #   
    #   #Create master matrix of jackknifed correlations:
    #   JK_corr_mat_pres <- data.frame(diag(4 * 4))
    #   tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
    #   tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
    #   rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
    #   
    #   JK_corr_mat_pres <- get_JK_corrmat_pres(dat_pres = dataset, 
    #                                           init_pres_coef = true_params[-c(1:4)])
    #   
    #   rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
    #   
    #   
    #   #Assign the Sigma_P as this master correlation matrix:
    #   Sigma_P = JK_corr_mat_pres
    #   
    # }
    
    set.seed(seed)
    
    #Generate the latent MVN variables:
    Lp = mvrnorm(length(unique(dataset$SUBJECT_ID)), numeric(16), Sigma = Sigma_P) 
    #output is a matrix, where the i-th row has elements for the i-th cluster.
    
    
    ######################################################
    ## Transform Lp --> F(Lp) ~ U(0,1): ##################
    ######################################################
    
    #Initialize:
    Up = Lp 
    
    #Then the marginal CDF for each element in Lp is the N(0,1):
    Up = apply(Lp, MARGIN = 1, pnorm)
    
    ############################################################
    ## ASSIGN PRESENCE CATEGORIES:                            ##
    ############################################################
    ## Calculate the cumulative probs.                        ##
    ## p0 = exp(alpha + X beta) / (1 + exp(alpha + X beta) )  ##
    ## p1 = 1 / exp(alpha + X beta)                           ##
    ## If Up < p0, then assign Y=0. Else Y=1                  ##
    ############################################################
    
    alpha = true_params["alpha"]
    beta_p = true_params[-c(1:4)]
    
    X = dataset[,-c(1,2)]
    
    linpred = alpha + as.matrix(X) %*% matrix(beta_p, nrow= (ncol(dataset)-2), ncol=1)
    
    
    p0 = matrix(as.vector(inv_logit_vec(linpred)), 
                byrow = T, 
                nrow=16, 
                ncol=length(unique(dataset$SUBJECT_ID)))
    
    bool_Wp_pres = (Up <= p0)
    
    sim_Wp_pres = ifelse(bool_Wp_pres, 0, 9) 
    #TRUE values get assigned category 0 (no fluorosis), FALSE values get assigned to dummy category 9 (fluorosis) 
    
    #Unlist sim_Wp_pres elements row-wise and populate the FRI_Score column in dataset. 
    #This is an intermediate step; this will later be re-updated with severity categories
    dataset[,"FRI_Score"] = as.vector(t(sim_Wp_pres)) 
    
    ##################################
    ###### SEVERITY CATEGORIES #######
    ##################################
    
    Sigma_S = diag(16) #initialize
    
    # if(corstr_sev=="independence")
    # { 
    #   #With independence cluster correlation structure
    #   Sigma_S = diag(16)
    #   
    # }else if(corstr_sev=="exchangeable")
    # { 
      #With exchangeable cluster correlation structure
      Sigma_S = matrix(exch_rho_sev, 16, 16); diag(Sigma_S) = 1
      
    # }else if(corstr_sev=="jackknifed")
    # {
    #   #With jackknifed cluster correlation structure
    #   
    #   JK_corr_mat_sev <- data.frame(diag(4 * 4))
    #   tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
    #   tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
    #   rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
    #   
    #   alphas = true_params[c(3,4)]
    #   betas = true_params[-c(1:4)] * true_params[1]
    #   init_sev_coef = c(alphas, betas)
    #   
    #   JK_corr_mat_sev <- get_JK_corrmat_sev(dat_sev = dataset, 
    #                                         init_sev_coef = init_sev_coef)
    #   
    #   rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
    #   
    #   Sigma_S = JK_corr_mat_sev 
    # }
    
    set.seed(seed)
    
    #Generate latent MVN variables:
    Ls = mvrnorm(length(unique(dataset$SUBJECT_ID)), numeric(16), Sigma = Sigma_S) 
    #the output is a matrix; the i-th row has elements for the i-th cluster.
    
    ######################################################
    ## Transform Ls --> F(Ls) ~ U(0,1): ##################
    ######################################################
    
    #Initialize:
    Us = Ls 
    
    #The the marginal CDF for each Lp is of the N(0,1)
    Us = apply(Ls, MARGIN = 1, pnorm)
    
    ############################################################
    ## ASSIGN SEVERITY CATEGORIES:                            ##
    ############################################################
    ## Calculate the cumulative probs.                        ##
    ## p1 = P(Ws<=1)                                          ##
    ## p2 = P(Ws<=2) - P(Ws<=1)                               ##
    ## p3 = 1 - P(Ws<=2)                                      ##
    ## where P(Ws<=l) = inv_logit(alpha_{l|l+1} + beta_S'X)   ##
    ## If  Wp=0: no updation                                  ##
    ##                                                        ##
    ## For Wp=1: if Us in (0, p1), then Y=1                   ##
    ##           if Us in (p1, p1+p2)  then Y=2               ##
    ##           if Us in (p1+p2, 1) then Y=3                 ##
    ############################################################
    
    alpha = true_params[c("alpha_1|2", "alpha_2|3")]
    gamma = true_params["gamma"]; beta_p = true_params[-c(1:4)]; beta_s = gamma*beta_p
    
    
    X = dataset[,-c(1,2)]
    
    linpred12 = alpha[1] + as.matrix(X) %*% matrix(beta_s, nrow= (ncol(dataset)-2), ncol=1)
    linpred23 = alpha[2] + as.matrix(X) %*% matrix(beta_s, nrow= (ncol(dataset)-2), ncol=1)
    
    cusum1 = matrix(as.vector(inv_logit_vec(linpred12)), byrow = T, nrow=16, ncol=length(unique(dataset$SUBJECT_ID)))
    cusum2 = matrix(as.vector(inv_logit_vec(linpred23)), byrow = T, nrow=16, ncol=length(unique(dataset$SUBJECT_ID)))
    
    p1 = cusum1
    p2 = cusum2 - cusum1
    p3 = 1 - cusum2
    
    
    bool_Y1_pres = ( Us <= p1 )
    bool_Y2_pres = ( Us > p1 & Us <= (p1+p2) )
    bool_Y3_pres = ( Us > p1+p2 )
    
    sim_Y_pres = sim_Wp_pres
    
    #If bool_Wp_pres==T i.e. fluorosis present and bool_Y1_pres==T, assign to category 1
    sim_Y_pres = ifelse(!bool_Wp_pres & bool_Y1_pres, 1, sim_Y_pres) 
    
    #If bool_Wp_pres==T i.e. fluorosis present and bool_Y1_pres==T, assign to category 2
    sim_Y_pres = ifelse(!bool_Wp_pres & bool_Y2_pres, 2, sim_Y_pres) 
    
    #If bool_Wp_pres==T i.e. fluorosis present and bool_Y1_pres==T, assign to category 3
    sim_Y_pres = ifelse(!bool_Wp_pres & bool_Y3_pres, 3, sim_Y_pres) 
    
    #Finally, update the dataset with the generated responses  
    dataset[,"FRI_Score"] = as.vector(t(sim_Y_pres)) #unlist rowwise and populate the FRI_Score column in dataset
    
    
    table(dataset$FRI_Score)
    
    return(dataset)
  }
  
  temp = update_MC_dataset2(MC_dataset, mc_seed)
  
  # if(length(unique(temp$FRI_Score))!=4)
  # {
  #   print("length(unique(temp$FRI_Score))!=4")
  #   
  # }
  
  
  if(all(c(0,1) %in% unique(temp$FRI_Score))){
    
    MC_dataset = temp 
    MC_output = fit_GEE_presence(dat=MC_dataset, corstr_pres = corstr_pres, kappa=0.25, maxIter = 10)
    
    coef_pres = MC_output$coefs
    rho_pres = 1
  
    JK_corr_mat_pres = NA
    
    if(corstr_pres=="ar1" || corstr_pres=="exchangeable"){  
      rho_pres = MC_output$rho_pres
      }else if(corstr_pres=="jackknifed")
      {
        JK_corr_mat_pres = MC_output$JK_corr_mat_pres
      }
    
    ##########################################################
    ## Get jackknifed estimates (for SD) and rhos (for CIs) ##
    ##########################################################
    
    coef_pres_JK = matrix( , nrow=length(true_params), ncol=N)
    rownames(coef_pres_JK) = names(true_params)
    rho_pres_JK = c()
    
    count = 1
    
    for(ID in unique(MC_dataset$SUBJECT_ID))
    {
      
      MC_dataset_JK = MC_dataset %>% filter(SUBJECT_ID != ID)
      
      #save only valid outputs : need 0 and 1 FRI scores present
      if(all(c(0,1) %in% unique(MC_dataset_JK$FRI_Score)))
      {
        MC_output_JK = 
          fit_GEE_presence(dat=MC_dataset_JK, corstr_pres = corstr_pres, kappa=0.25, maxIter = 10,
                           reuse = T, rho_pres = rho_pres, JK_corr_mat_pres = JK_corr_mat_pres)
        
        #store JK outputs:
        
        # if(corstr_pres=="ar1" || corstr_pres=="exchangeable"){
        #   rho_pres_JK = c(rho_pres_JK, MC_output_JK$rho_pres)
        # }
        
        coef_pres_JK[ -c(1,3,4), count] = MC_output_JK$coefs$Estimates
      }
      
      count = count+1
      
      print(ID)
    } 
    
    #Get SD:
    JK_SD_pres = data.frame(apply(coef_pres_JK, MARGIN=1, FUN=compute_JK_SE)) 
    colnames(JK_SD_pres) = "SD"
    JK_SD_pres = JK_SD_pres[-c(1,3,4), , drop=F]
    
    #Convert SD to SE and Get standardized estimates:
    #std_coef_pres = coef_pres[,2]/JK_SD_pres
    std_coef_pres = (coef_pres[,2]/JK_SD_pres)/sqrt(N)
    colnames(std_coef_pres) = "Standardized Estimate"
    
    #save this as an RData list...instead of all this
    
    setwd(paste0("/blue/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode, "/",
                 corstr_pres, ",", corstr_sev, "/age", age))
    
    
    save(coef_pres, file = paste0("coef_pres_MC_", mc_seed, ".Rdata") )
    save(JK_SD_pres, file = paste0("JK_SD_MC_", mc_seed, ".Rdata") )
    save(std_coef_pres, file = paste0("std_coef_pres_MC_", mc_seed, ".Rdata") )
    #save(coef_pres_JK, file = paste0("coef_pres_JK_MC_", mc_seed, ".Rdata") )
    
    if(corstr_pres=="ar1" || corstr_pres=="exchangeable")
    { save(rho_pres, file = paste0("rho_pres_MC_", mc_seed, ".Rdata") )
      #save(rho_pres_JK, file = paste0("rho_pres_JK_MC_", mc_seed, ".Rdata") )
    }
    
    
  }else{
    print("MC dataset lacks both 0/1 FRI Scores.")
  }
  
}


end = Sys.time()

elapsed = end - start; elapsed
