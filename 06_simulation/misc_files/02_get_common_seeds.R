

inSLURM = !is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))

path_prefix = "W:/"

# Check if running in a SLURM environment
if (inSLURM) {
  # If in SLURM environment
  path_prefix = "/blue/"
} 

get_common_MC_estimates = function(N, corstr_pres, corstr_sev)
{
  
  setwd(paste0(path_prefix, "somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation"))
  
  mode = "presence"
  filename = paste0("validMCseeds_", corstr_pres, "_", corstr_sev, "_N", N, "_", mode, ".RData")
  pobj = load(filename); assign("seeds_pres", get(pobj))
  
  mode = "severity"
  filename = paste0("validMCseeds_", corstr_pres, "_", corstr_sev, "_N", N, "_", mode, ".RData")
  sobj = load(filename); assign("seeds_sev", get(sobj))
  
  ages = c(9,13,17,23)
  
  selected_seeds = list()
  
  for(i in 1:length(ages))
  {
    
    #select first 100 seeds that produce complete MC samples
    selected_seeds[[i]] = intersect(seeds_pres[[i]], seeds_sev[[i]])[1:100]
  }
  
  names(selected_seeds) = paste0("age", ages)
  
  #Now: Get the estimates...
  
  ###################################
  ## Extract the variable names... ##
  ###################################
  
  mode = "presence"; age=9
  setwd(paste0(path_prefix, "somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode, "/",
               corstr_pres, ",", corstr_sev, "/age", age))
  mc_seed = 89; file = paste0("std_coef_pres_MC_", mc_seed, ".Rdata")
  obj = load(file); assign("std_coef_pres", get(obj))
  param_names = c("rho", "gamma", "alpha", "alpha_1|2", "alpha_2|3", rownames(std_coef_pres)[-c(1)]) #variable names extracted
  
  #############################################################
  ## Set up lists to store the estimates from the MC samples ##
  #############################################################
  
  MC_estimates_temp_df = data.frame(matrix(NA, nrow=100, ncol=(1+length(param_names))))
  colnames(MC_estimates_temp_df) = c("MC_seed", param_names)
  
  MC_estimates_presence_list = list(age9 = MC_estimates_temp_df, 
                                    age13 = MC_estimates_temp_df, 
                                    age17 = MC_estimates_temp_df, 
                                    age23 = MC_estimates_temp_df)
  
  MC_estimates_severity_list = MC_estimates_presence_list
  
  ############################
  ## Extract the estimates: ##
  ############################
  
  for(i in 1:length(ages))
  {
    age = ages[i]
    
    MC_estimates_presence_list[[i]]$MC_seed = 
      MC_estimates_severity_list[[i]]$MC_seed = selected_seeds[[i]] #fills in the selected seeds into the MC_seed column
    
    for(seed_ind in 1:length(selected_seeds[[i]]))
    {
      mc_seed = selected_seeds[[i]][seed_ind] #choose the current MC seed
      
      ################
      ### Presence ###
      ################
      
      mode="presence"
      setwd(paste0(path_prefix, "somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode, "/",
                   corstr_pres, ",", corstr_sev, "/age", age))
      
      # Load estimates:
      file = paste0("std_coef_pres_MC_", mc_seed, ".Rdata")
      
      obj = load(file); assign("std_coef_obj", get(obj))
      age_ind = which(ages %in% age)
      pres_est = as.vector(std_coef_obj[,1]); names(pres_est) = rownames(std_coef_obj)
      
      #Rho:
      
      pres_rho = 1
      
      if(corstr_pres %in% c("exchangeable", "ar1"))
      {
        file = paste0("rho_pres_MC_", mc_seed, ".Rdata")
        
        obj = load(file); assign("rho_obj", get(obj))
        pres_rho = rho_obj
      }
      
      ################
      ### Severity ###
      ################
      
      mode="severity"
      setwd(paste0(path_prefix, "somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N", N, "/", mode, "/",
                   corstr_pres, ",", corstr_sev, "/age", age))
      
      #Load estimates:
      file = paste0("std_coef_sev_MC_", mc_seed, ".Rdata")
      
      obj = load(file); assign("std_coef_obj", get(obj))
      age_ind = which(ages %in% age)
      sev_est = as.vector(std_coef_obj[,1]); names(sev_est) = rownames(std_coef_obj)
      
      #Rho:
      
      sev_rho = 1
      
      if(corstr_sev %in% c("exchangeable", "ar1"))
      {
        file = paste0("rho_sev_MC_", mc_seed, ".Rdata")
        
        obj = load(file); assign("rho_obj", get(obj))
        sev_rho = rho_obj
      }
      
      
      #Fill in the lists...
      
      # Presence:
      MC_estimates_presence_list[[i]][seed_ind, -c(1,2,3,5,6)] = pres_est #intercept and beta parameters
      MC_estimates_presence_list[[i]][seed_ind, c(2)] = pres_rho #rho
      
      # Severity:
      MC_estimates_severity_list[[i]][seed_ind, -c(1,2,3,4)] = sev_est #intercept and beta parameters 
      MC_estimates_severity_list[[i]][seed_ind, c(2)] = sev_rho #rho
      
    }
    
    print(paste0("Age ", age, " completed."))
  }
  
  output_list = list(MC_estimates_presence_list=MC_estimates_presence_list,
                     MC_estimates_severity_list=MC_estimates_severity_list)
  
  setwd(paste0(path_prefix, "somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation"))
  
  save(output_list, file= paste0("common_MC_estimates_N", N, "_", 
                             corstr_pres, "_", corstr_sev, ".RData"))
  
}

get_common_MC_estimates(N = 30, corstr_pres = "exchangeable", corstr_sev = "independence")
get_common_MC_estimates(N = 50, corstr_pres = "exchangeable", corstr_sev = "independence")


#Get summaries for these sev and pres pieces.


#Turn this into a function that inputs
## N = 30, mode = "presence", corstr_pres = "exchangeable", corstr_sev = "independence"
#and saves it as an RData file.

# then have a mc_seed based script that takes in these inputs one by one, and passes them to the combined script.



#once these are populated...feed them into combined estimation scripts later...
#and get their col means and SD for estimate and SE.
