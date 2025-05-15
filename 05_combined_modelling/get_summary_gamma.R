

#corstr_pres = "independence"
#corstr_sev = "independence"

library(openxlsx)

# Check if running in a SLURM environment
if (!is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))) {
  # If in SLURM environment
  setwd("/blue/somnath.datta/shoumisarkar/Fluorosis/")
} else {
  # If not in SLURM environment
  setwd("W:/somnath.datta/shoumisarkar/Fluorosis/")
}


ages = c(9, 13, 17, 23)

gamma_df = data.frame(matrix(numeric(16*4),16, 4))
colnames(gamma_df) = c("corstr_pres", "corstr_sev", "age", "Estimate")
corstrs = c("independence", "exchangeable", "ar1", "jackknifed")

count = 1

for(corstr in corstrs)
{
  for(age in ages)
  {
    
    filepath = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/05a_combined_presence_modelling/",
                      corstr, ",", corstr, "/whole_data_based/coefs_pres_age", age, ".RData")
    
    obj_coef <- load(filepath)
    assign("coef_pres", get(obj_coef))
    
    gamma_df[count,1:2] = corstr
    gamma_df[count,3] = age
    gamma_df[count,4] = coef_pres[1,2]
    #gamma_df[count,4] = round(coef_pres[1,2], digits = 2)
    
    count = count+1
  }
}

#View(gamma_df)

library(writexl)

filepath2 = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/05a_combined_presence_modelling/summary_gamma.xlsx")
write_xlsx(gamma_df, path = filepath2)
