#Convert the std_est into 100*std_est

#Presence:

age=9

ages = c(9,13,17,23)



for(age in ages)
{
  setwd(paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N10000/severity/independence,independence/age", age))
  
  obj = load("std_coef_sev_MC_2_ver0.RData")
  assign("std_coef_sev", get(obj))
  
  std_coef_sev$`Standardized Estimate` = std_coef_sev$`Standardized Estimate`/100
  
  save(std_coef_sev, file = "std_coef_sev_MC_2.RData")
}