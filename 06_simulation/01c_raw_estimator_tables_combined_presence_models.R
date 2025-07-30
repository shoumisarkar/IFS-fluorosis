library(dplyr)
library(openxlsx)

path_prefix = "W:/"

# Check if running in a SLURM environment
if (!is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))) {
  # If in SLURM environment
  path_prefix = "/blue/"
}

setwd(paste0(path_prefix, "path/to/Fluorosis/"))
source(file = "Codes/functions.R")

model = "combined"

corstr_pres = "jackknifed"
corstr_sev = "jackknifed"

#corstr_pres = "exchangeable"
#corstr_sev = "exchangeable" 

# corstr_pres = "independence"
# corstr_sev = "exchangeable" 

# corstr_pres = "exchangeable"
# corstr_sev = "independence"

mc_seed_range = 1:100


for(N in c(30, 50, 200))
{
  
  ############################
  ### Get parameter names: ###
  ############################
  
  param_names = c()
  ages = c(9,13,17,23)
  
  
  for(mc_seed in mc_seed_range)
  {
    setwd(paste0(path_prefix, "path/to/Fluorosis/Results/06_simulation/N", N, "/", 
                 model, "/", corstr_pres, ",", corstr_sev, "/age9"))
    
    fn_est = paste0("coef_pres_MC_", mc_seed, ".Rdata")
    
    #if(model=="severity"){
    #  fn_est = paste0("coef_sev_MC_", mc_seed, ".Rdata")
    #}
    
    if(file.exists(fn_est))
    {
      obj = load(fn_est)
      assign("temp_est", get(obj))
      
      param_names = temp_est$Variables
      
      break
    }
  }
  
  #we want to keep all param names, so comment this out
  #param_names = param_names[-c(1:4)]
  
  ###################################################
  ### Extract estimates into a list of dataframes ###
  ###################################################
  
  # Create an empty list of size 4
  est_list <- vector("list", length = 4) ; stdest_list = rho_list = est_list
  
  for(age in ages)
  {
    age_ind = which(ages %in% age)
    
    setwd(paste0(path_prefix, "path/to/Fluorosis/Results/06_simulation/N", N, "/", 
                 model, "/", corstr_pres, ",", corstr_sev, "/age", age))
    
    
    #Initialize data frames to store the outputs
    est_df = data.frame(matrix(NA, nrow=length(param_names), ncol=length(mc_seed_range)))
    rownames(est_df) = param_names
    colnames(est_df) = paste0("MC", mc_seed_range)
    
    stdest_df = est_df[-c(1),]
    
    pres_rho_vec = c()
    sev_rho_vec = c()
    
    for(mc_seed in mc_seed_range)
      
    {
      fn_est = paste0("coef_pres_MC_", mc_seed, ".Rdata")
      fn_stdest = paste0("std_coef_pres_MC_", mc_seed, ".Rdata")
      fn_rho = paste0("rho_pres_MC_", mc_seed, ".Rdata")
      
      if(file.exists(fn_stdest))
      {
      
        # Safely load 'fn_est'
        tryCatch({
          obj <- load(fn_est)
          assign("temp_est", get(obj), envir = .GlobalEnv)
          est_df[, mc_seed] <- temp_est$Estimates  #save all estimates: gamma, alpha, alpha_1|2, alpha_2|3, beta_1,..., beta_q
        }, error = function(e) {
          warning("Failed to load or assign from fn_est: ", e$message)
          est_df[, mc_seed] <- NA  # Assign NA to this column if there's an error
        })
        
        
        # Safely load 'fn_stdest'
        tryCatch({
          obj <- load(fn_stdest)
          assign("temp_stdest", get(obj), envir = .GlobalEnv)
          stdest_df[-c(2,3), mc_seed] <- temp_stdest[,1]  #save for all estimates except (gamma already dropped), alpha_1|2, alpha_2|3: alpha, beta_1,..., beta_q
        }, error = function(e) {
          warning("Failed to load or assign from fn_stdest: ", e$message)
          stdest_df[, mc_seed] <- NA  # Assign NA to this column if there's an error
        })
        
      }
    }
    
    est_list[[age_ind]] = est_df
    stdest_list[[age_ind]] = stdest_df
    
    print(paste0("Age ", age, " completed..."))
  }
  
  names(est_list) = names(stdest_list) = paste0("age", ages)
  
  
  ########################
  ###   J-S shrinkage  ###
  ########################
  
  JS_stdest_list = stdest_list
  
  #define for only #save for only: beta_1,..., beta_q ; drop: (gamma already dropped), alpha, alpha_1|2, alpha_2|3
  JS_stdest_list$age9 = JS_stdest_list$age9[-c(1,2,3),] 
  JS_stdest_list$age13 = JS_stdest_list$age13[-c(1,2,3),]
  JS_stdest_list$age17 = JS_stdest_list$age17[-c(1,2,3),]
  JS_stdest_list$age23 = JS_stdest_list$age23[-c(1,2,3),]
  
  
  #loop over vars, the loop over mc_seed
  
  vars = param_names[-c(1,2,3,4)] #keep only the beta_1,...,beta_q (drop gamma, alpha, alpha_1|2, alpha_2|3)
  
  n_intercept = 3
  

  for(var in vars)
  {
    var_ind = which(vars %in% var) #+ n_intercept 
    
    for(mc_seed in mc_seed_range)
    {
      t <- c(
        stdest_list$age9[(var_ind+n_intercept),  mc_seed],
        stdest_list$age13[(var_ind+n_intercept), mc_seed],
        stdest_list$age17[(var_ind+n_intercept), mc_seed],
        stdest_list$age23[(var_ind+n_intercept), mc_seed]
      )
      
      t = js(t)
      
      JS_stdest_list$age9[var_ind,mc_seed] = t[1]
      JS_stdest_list$age13[var_ind,mc_seed] = t[2]
      JS_stdest_list$age17[var_ind,mc_seed] = t[3]
      JS_stdest_list$age23[var_ind,mc_seed] = t[4]
      
    }
    
  }
  
  # JS_stdest_list$age9[1:n_intercept,] = NA
  # JS_stdest_list$age13[1:n_intercept,] = NA
  # JS_stdest_list$age17[1:n_intercept,] = NA
  # JS_stdest_list$age23[1:n_intercept,] = NA
  
  
  #############################
  ### Compute summary table ###
  #############################
  
  #For each age, want stdest, their SD (SE), JS_stdest, their SD
  
  # Create an empty list of size 4
  summary_list <- vector("list", length = 4) 
  
  for(age in ages)
  {
    age_ind = which(ages %in% age)
    
    est = apply(est_list[[age_ind]], MARGIN = 1, function(x){mean(x, na.rm = T)  })#[-c(1:n_intercept)]
    SE_est = apply(est_list[[age_ind]], MARGIN = 1, function(x){sd(x, na.rm = T)  }) #[-c(1:n_intercept)]

    temp_df = data.frame(Variable = trueTarget$Variables,
         
                         #Pre_JS_Estimate = stdest[-c(1:3)],
                         #Pre_JS_SE = SE_stdest[-c(1:3)],
                         Target = trueTarget$Estimates,
                         Estimate = est[-c(1:4)],
                         SE = SE_est[-c(1:4)])
    
    temp_df$Bias = temp_df$Estimate - temp_df$Target
    temp_df$MSE = (temp_df$Bias)^2 + (temp_df$SE)^2
    
    temp_df = temp_df[,c("Variable", 
                         "Estimate", "Bias", "SE", "MSE")]
    
    temp_df[-1] <- round(as.numeric(as.character(unlist(temp_df[-1]))), 3)  # Round the numeric columns
    
    summary_list[[age_ind]] = temp_df
  }
  
  names(summary_list) = paste0("age", ages)
  
  
  ########################
  ### Save the outputs ###
  ########################
  
  setwd(paste0(path_prefix, "path/to/Fluorosis/Results/06_simulation/N", N, "/", 
               model, "/"))
  
  save(summary_list, file=paste0("combined_presence_summarytable_raw_est_", model,"_", corstr_pres, ",", corstr_sev, ".Rdata"))
  
  write.xlsx(summary_list, file = paste0("summarytable_raw_est_", model,"_presence_", corstr_pres, ",", corstr_sev, "_N_", N, ".xlsx"), 
             sheetNames = names(summary_list))

}
