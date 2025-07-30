
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

# corstr_pres = "jackknifed"
# corstr_sev = "jackknifed"

#corstr_pres = "exchangeable"
#corstr_sev = "exchangeable" 

# corstr_pres = "independence"
# corstr_sev = "exchangeable" 

corstr_pres = "exchangeable"
corstr_sev = "independence"


models = c("presence", "severity")

for(model in models)
{
  print(model)
  
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
      
      if(model=="severity"){
        fn_est = paste0("coef_sev_MC_", mc_seed, ".Rdata")
      }
      
      if(file.exists(fn_est))
      {
        obj = load(fn_est)
        assign("temp_est", get(obj))
        
        param_names = temp_est$Variables
        
        break
      }
    }
    
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
      rownames(est_df) = temp_est$Variables
      colnames(est_df) = paste0("MC", mc_seed_range)
      
      stdest_df = est_df
      
      rho_vec = c()
      
      for(mc_seed in mc_seed_range)
        
      {
        fn_est = paste0("coef_pres_MC_", mc_seed, ".Rdata")
        fn_stdest = paste0("std_coef_pres_MC_", mc_seed, ".Rdata")
        fn_rho = paste0("rho_pres_MC_", mc_seed, ".Rdata")
        
        if(model=="severity")
        {
          fn_est = paste0("coef_sev_MC_", mc_seed, ".Rdata")
          fn_stdest = paste0("std_coef_sev_MC_", mc_seed, ".Rdata")
          fn_rho = paste0("rho_sev_MC_", mc_seed, ".Rdata")
        }
        
        if(file.exists(fn_stdest))
        {
          
          
          # Safely load 'fn_est'
          tryCatch({
            obj <- load(fn_est)
            assign("temp_est", get(obj), envir = .GlobalEnv)
            est_df[, mc_seed] <- temp_est$Estimates  # Assuming est_df and mc_seed are correctly defined
          }, error = function(e) {
            warning("Failed to load or assign from fn_est: ", e$message)
            est_df[, mc_seed] <- NA  # Assign NA to this column if there's an error
          })
          
          
          # Safely load 'fn_stdest'
          tryCatch({
            obj <- load(fn_stdest)
            assign("temp_stdest", get(obj), envir = .GlobalEnv)
            stdest_df[, mc_seed] <- temp_stdest$`Standardized Estimate`  # Assuming stdest_df and mc_seed are correctly defined
          }, error = function(e) {
            warning("Failed to load or assign from fn_stdest: ", e$message)
            stdest_df[, mc_seed] <- NA  # Assign NA to this column if there's an error
          })
          
          
          if(model=="presence")
          {
            if(corstr_pres %in% c("exchangeable", "ar1"))
            {
              obj = load(fn_rho)
              assign("rho", get(obj))
              rho_vec = c(rho_vec, rho)
            }
          }else if(model=="severity")
          {
            if(corstr_sev %in% c("exchangeable", "ar1"))
            {
              obj = load(fn_rho)
              assign("rho", get(obj))
              rho_vec = c(rho_vec, rho)
            }
          }
          
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
    
    #loop over vars, the loop over mc_seed
    
    vars = param_names; n_intercept = 1
    
    if(model=="severity"){
      n_intercept = 2
    }
    
    vars = param_names[-c(1:n_intercept)]
    
    for(var in vars)
    {
      var_ind = n_intercept + which(vars %in% var)
      
      for(mc_seed in mc_seed_range)
      {
        t <- c(
          stdest_list$age9[var_ind,mc_seed],
          stdest_list$age13[var_ind,mc_seed],
          stdest_list$age17[var_ind,mc_seed],
          stdest_list$age23[var_ind,mc_seed]
        )
        
        t = js(t)
        
        JS_stdest_list$age9[var_ind,mc_seed] = t[1]
        JS_stdest_list$age13[var_ind,mc_seed] = t[2]
        JS_stdest_list$age17[var_ind,mc_seed] = t[3]
        JS_stdest_list$age23[var_ind,mc_seed] = t[4]
        
      }
      
    }
    
    JS_stdest_list$age9[1:n_intercept,] = NA
    JS_stdest_list$age13[1:n_intercept,] = NA
    JS_stdest_list$age17[1:n_intercept,] = NA
    JS_stdest_list$age23[1:n_intercept,] = NA
    
    
    #############################
    ### Compute summary table ###
    #############################
    
    # Create an empty list of size 4
    summary_list <- vector("list", length = 4) 
    
    for(age in ages)
    {
      age_ind = which(ages %in% age)
      
      est = apply(est_list[[age_ind]], MARGIN = 1, function(x){mean(x, na.rm = T)  })[-c(1:n_intercept)]
      SE_est = apply(est_list[[age_ind]], MARGIN = 1, function(x){sd(x, na.rm = T)  })[-c(1:n_intercept)]
      
      prefix = ifelse(model=="severity", "sev", "pres")
      
      trueTarget = load(paste0(path_prefix, "path/to/Fluorosis/Results/06_simulation/N10000/",
                               model, "/independence,independence/age", age, "/coef_", prefix, "_MC_2.Rdata"))
      
      assign("trueTarget", get(trueTarget))
      
      if(model=="presence")
      {
        trueTarget = trueTarget[-c(1),]
      }else if(model=="severity")
      {
        trueTarget = trueTarget[-c(1,2),]
      }
      
      temp_df = data.frame(Variable = names(SE_est),
                           Estimate = est,
                           SE = SE_est,

                           # JS_StdEstimate = JS_stdest,
                           # JS_SE = SE_JS_stdest)
                           
                           Target = trueTarget$Estimates
                           
                           #Estimate = JS_stdest,
                           #SE = SE_JS_stdest
      )
      
      temp_df$Bias = temp_df$Estimate - temp_df$Target
      temp_df$MSE = (temp_df$Bias)^2 + (temp_df$SE)^2
      
      temp_df = temp_df[,c("Variable",  
                           "Estimate",
                           "Bias", "SE", "MSE")]
      
      temp_df[-1] <- round(as.numeric(as.character(unlist(temp_df[-1]))), 3)  # Round the numeric columns
      
      summary_list[[age_ind]] = temp_df
    }
    
    names(summary_list) = paste0("age", ages)
    
    
    ########################
    ### Save the outputs ###
    ########################

    setwd(paste0(path_prefix, "path/to/Fluorosis/Results/06_simulation/N", N, "/",
                 model, "/"))

    write.xlsx(summary_list, file = paste0("summarytable_raw_est_", model,"_", corstr_pres, ",", corstr_sev, "_N_", N, ".xlsx"),
               sheetNames = names(summary_list))

  
