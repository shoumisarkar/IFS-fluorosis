
corstr_pres = "ar1"
corstr_sev = "ar1"
age = 23

library(writexl)

present_b = c()

for(b in 1:100)
  
{
  path = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/05b_combined_severity_modelling/",
                corstr_pres, ",", corstr_sev, "/bootstrapping/age", age, "/")
  
  fn_stdcoef_pres = paste0(path, "std_coefs_sev_age", age, "_b", b, ".RData")
  
  if(!file.exists(fn_stdcoef_pres))
  {
    next
  }
  
  present_b = c(present_b, b)
}

library(dplyr)
missing_b = setdiff(1:100, present_b); missing_b
length(missing_b)
