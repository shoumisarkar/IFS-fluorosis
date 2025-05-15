library(openxlsx)

# Check if running in a SLURM environment
if (!is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))) {
  # If in SLURM environment
  setwd("/blue/somnath.datta/shoumisarkar/Fluorosis/")
} else {
  # If not in SLURM environment
  setwd("W:/somnath.datta/shoumisarkar/Fluorosis/")
}

source(file = "Codes/functions.R")


# Helper function to compute CI
get_CI <- function(confidence_level = 0.95, BS_dat, vars) {
  lq <- (1 - confidence_level) / 2
  uq <- 1 - lq
  
  data.frame(
    Variable = vars,
    LCL = apply(BS_dat, MARGIN = 1, FUN = function(x) quantile(x, lq, na.rm = TRUE)),
    UCL = apply(BS_dat, MARGIN = 1, FUN = function(x) quantile(x, uq, na.rm = TRUE))
  )
}

# Helper function to write CI results to Excel
write_CI_to_excel_preJS <- function(CI_list, filename, corstr_pres, conf_level) {
  wb <- createWorkbook()
  
  for (name in names(CI_list)) {
    addWorksheet(wb, name)
    writeData(wb, sheet = name, CI_list[[name]])
  }
  
  saveWorkbook(wb, paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/preJS_",
                          conf_level, "_CI_", corstr_pres, ".xlsx"), overwrite = TRUE)
}

# Helper function to write CI results to Excel
write_CI_to_excel_postJS <- function(CI_list, filename, corstr_pres, conf_level) {
  wb <- createWorkbook()
  
  for (name in names(CI_list)) {
    addWorksheet(wb, name)
    writeData(wb, sheet = name, CI_list[[name]])
  }
  
  saveWorkbook(wb, paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/postJS_",
                          conf_level, "_CI_", corstr_pres, ".xlsx"), overwrite = TRUE)
}

# Function to perform bootstrapping and save CI results
process_age_data <- function(age, B, corstr_pres, get_CI, vars) {
  stdcoefs_BS_mat <- matrix(NA, nrow = 13, ncol = B)
  path <- paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/",
                 corstr_pres, "/bootstrapping/age", age, "/")
  
  for (b in 1:B) {
    fn_stdcoef_pres <- paste0(path, "std_coefs_pres_age", age, "_b", b, ".RData")
    if (file.exists(fn_stdcoef_pres)) {
      obj_stdcoef_pres <- load(fn_stdcoef_pres)
      assign("std_coef_pres", get(obj_stdcoef_pres))
      stdcoefs_BS_mat[, b] <- std_coef_pres$`Standardized Estimate`[-1]  # Adjusted to remove first row
    }
  }
  
  return(stdcoefs_BS_mat)
}

# Main logic
corstr_pres <- "independence"
ages <- c(9, 13, 17, 23)
B <- 120
vars <- c("dental_age", "Total_mgF", "SugarAddedBeverageOzPerDay", "BrushingFrequencyPerDay",
          "Avg_homeppm", "Prop_DentAppt", "Prop_FluorideTreatment", "Tooth8", "Tooth9",
          "Tooth10", "ZoneM", "ZoneI", "ZoneO")

stdcoefs_BS_list <- lapply(ages, process_age_data, B = B, corstr_pres = corstr_pres, get_CI = get_CI, vars = vars)
names(stdcoefs_BS_list) <- paste0("age", ages)

# Save the pre-JS standardized coefficients
save(stdcoefs_BS_list, file = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/preJS_stdcoefBS_",
                                     corstr_pres, ".RData"))

# Calculate CIs
preJS_CI_95 <- lapply(stdcoefs_BS_list, get_CI, confidence_level = 0.95, vars = vars)
preJS_CI_90 <- lapply(stdcoefs_BS_list, get_CI, confidence_level = 0.90, vars = vars)

# Write CIs to Excel files
write_CI_to_excel_preJS(preJS_CI_95, "preJS_stdcoefBS", corstr_pres, 95)
write_CI_to_excel_preJS(preJS_CI_90, "preJS_stdcoefBS", corstr_pres, 90)

# Repeat the process for post-James-Stein
stdcoefsJS_BS_list <- stdcoefs_BS_list  # will be updated below, following the initialized structure

#vars = rownames(std_coef_pres)[-c(1)]

for(var in vars)
{
  var_ind = which(vars %in% var)
  
  for(b in 1:B)
  {
    t <- c(
      (stdcoefs_BS_list[[1]])[var_ind, b],
      (stdcoefs_BS_list[[2]])[var_ind, b],
      (stdcoefs_BS_list[[3]])[var_ind, b],
      (stdcoefs_BS_list[[4]])[var_ind, b]
    )
    
    t = js(t)
    
    stdcoefsJS_BS_list[[1]][var_ind,b] = t[1]
    stdcoefsJS_BS_list[[2]][var_ind,b] = t[2]
    stdcoefsJS_BS_list[[3]][var_ind,b] = t[3]
    stdcoefsJS_BS_list[[4]][var_ind,b] = t[4]
    
  }
}


postJS_CI_95 <- lapply(stdcoefsJS_BS_list, get_CI, confidence_level = 0.95, vars = vars)
postJS_CI_90 <- lapply(stdcoefsJS_BS_list, get_CI, confidence_level = 0.90, vars = vars)

# Write post-JS CIs to Excel files
write_CI_to_excel_postJS(postJS_CI_95, "postJS_stdcoefBS", corstr_pres, 95)
write_CI_to_excel_postJS(postJS_CI_90, "postJS_stdcoefBS", corstr_pres, 90)
