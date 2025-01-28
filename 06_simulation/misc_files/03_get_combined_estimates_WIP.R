

mc_seed = 1
age=9

N=30

corstr_pres="exchangeable"
corstr_sev="independence"

exch_rho_pres = 0.6 
exch_rho_sev = 0.6 

########################################################

inSLURM = !is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))

path_prefix = "W:/"

# Check if running in a SLURM environment
if (inSLURM) {
  # If in SLURM environment
  path_prefix = "/blue/"
} 

###########################################################
#### Load the estimates from the piecewise estimations ####
###########################################################

setwd(paste0(path_prefix, "somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation"))

file = paste0("common_MC_estimates_N", N, "_", corstr_pres, "_", corstr_sev, ".RData")

obj = load(file); assign("common_estimates", get(obj))


age_ind = which(c(9,13,17,23) %in% age)


pres_coef = unlist(common_estimates$MC_estimates_presence_list[[age_ind]][mc_seed, -c(1,2,3,5,6)])
pres_rho = unlist(common_estimates$MC_estimates_presence_list[[age_ind]][mc_seed, 2])

sev_coef = unlist(common_estimates$MC_estimates_severity_list[[age_ind]][mc_seed, -c(1,2,3,4)])
sev_rho = unlist(common_estimates$MC_estimates_severity_list[[age_ind]][mc_seed, 2])


###########################################################
############# Create the MC dataset #######################
###########################################################


setwd(paste0(path_prefix, "somnath.datta/shoumisarkar/Fluorosis"))
source("Codes/functions.R")

print(paste0("Running with MC seed ", mc_seed))

#Load the data for the chosen age group:
dat <- read_excel(paste0("Results/02_preprocess_data/preprocessed_IFS_data_age", age, ".xlsx"))

#Keep selected columns:
dat <- dat %>% dplyr::select(SUBJECT_ID, FRI_Score, dental_age, Avg_homeppm, 
                             SugarAddedBeverageOzPerDay, BrushingFrequencyPerDay, 
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

obj = load(paste0(path_prefix, "somnath.datta/shoumisarkar/Fluorosis/Results/05a_combined_presence_modelling/exchangeable,exchangeable/whole_data_based/coefs_pres_age",age, ".RData"))
assign("coef_obj", get(obj))

true_params = coef_obj[coef_obj$Variables %in% sim_varnames,2] #these are the true values we need.
names(true_params) = sim_varnames


update_MC_dataset = function(dataset=MC_dataset, seed=mc_seed)
{
  #Update the Y in dataset:
  #####################################################
  ##  Set up Sigma, the cluster correlation matrix:  ##
  #####################################################
  
  #Initialize:
  Sigma_P = diag(16) 
  
  if(corstr_pres=="independence")
  { 
    #With independence cluster correlation structure
    Sigma_P = diag(16)
  }else if(corstr_pres=="exchangeable")
  { 
    #With exchangeable cluster correlation structure
    Sigma_P = matrix(exch_rho_pres, 16, 16); diag(Sigma_P) = 1
  }
  if(corstr_pres=="jackknifed")
  {
    
    #Create master matrix of jackknifed correlations:
    JK_corr_mat_pres <- dataset.frame(diag(4 * 4))
    tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
    tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
    rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
    
    JK_corr_mat_pres <- get_JK_corrmat_pres(dat_pres = dataset, 
                                            init_pres_coef = true_params[-c(1:4)])
    
    rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
    
    
    #Assign the Sigma_P as this master correlation matrix:
    Sigma_P = JK_corr_mat_pres
    
  }
  
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
  
  linpred = alpha + as.matrix(X) %*% matrix(beta_p, nrow=10, ncol=1)
  
  
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
  
  if(corstr_sev=="independence")
  { 
    #With independence cluster correlation structure
    Sigma_S = diag(16)
    
  }else if(corstr_sev=="exchangeable")
  { 
    #With exchangeable cluster correlation structure
    Sigma_S = matrix(exch_rho_sev, 16, 16); diag(Sigma_P) = 1
    
  }else if(corstr_sev=="jackknifed")
  {
    #With jackknifed cluster correlation structure
    
    JK_corr_mat_sev <- dataset.frame(diag(4 * 4))
    tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
    tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
    rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
    
    alphas = true_params[c(3,4)]
    betas = true_params[-c(1:4)] * true_params[1]
    init_sev_coef = c(alphas, betas)
    
    JK_corr_mat_sev <- get_JK_corrmat_sev(dat_sev = dataset, 
                                          init_sev_coef = init_sev_coef)
    
    rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
    
    Sigma_S = JK_corr_mat_sev 
  }
  
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
  
  linpred12 = alpha[1] + as.matrix(X) %*% matrix(beta_s, nrow=10, ncol=1)
  linpred23 = alpha[2] + as.matrix(X) %*% matrix(beta_s, nrow=10, ncol=1)
  
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

temp = update_MC_dataset(MC_dataset, mc_seed)

if(length(unique(temp$FRI_Score))!=4)
{
  print("length(unique(temp$FRI_Score))!=4")
  
}else{

  MC_dataset = temp
  
  fit = fit_GEE_combined(dat=MC_dataset, corstr_pres = corstr_pres, corstr_sev = corstr_sev,
                         init_pres_coef = pres_coef, init_sev_coef = sev_coef
  )
  
  ## Write the code here to collect coef_pres, coef_sev, rho_pres, rho_sev
  ## Refer to 01a_sim_dat_MC_N_30_presence
  ## carry out jackknifing... init df for pres est, another df for sev est.
  ## also a vector for rho_pres_jk, rho_sev_jk
  ## JK, SD... the starting values are the same as the main estimates' ; since we dont have them stored
  
}