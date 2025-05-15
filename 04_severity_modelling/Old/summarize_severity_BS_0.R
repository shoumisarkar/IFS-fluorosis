
corstr_sev = "exchangeable"

library(writexl)

ages = c(9, 13, 17, 23)
B=120
stdcoefs_BS_list = list()

for(age in ages)
{
  stdcoefs_BS_mat = matrix(, nrow=13, ncol=B)
  
  ind_age = which(ages %in% age)
  
  for(b in 1:B)
    
  {
    path = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/",
                  corstr_sev, "/bootstrapping/age", age, "/")
    
    fn_stdcoef_sev = paste0(path, "std_coefs_sev_age", age, "_b", b, ".RData")
    
    if(!file.exists(fn_stdcoef_sev))
    {
      next
    }
    
    obj_stdcoef_sev = load(fn_stdcoef_sev); assign("std_coef_sev", get(obj_stdcoef_sev))
    
    stdcoefs_BS_mat[,b] = std_coef_sev$`Standardized Estimate`[-c(1,2)]
  }
  
  stdcoefs_BS_list[[ind_age]] = stdcoefs_BS_mat
  names(stdcoefs_BS_list)[ind_age] = paste0("age", age)
}

#Standardized coefs before James-Stein updation:
filename_stdcoefs_BS_list = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/preJS_stdcoefBS_",
                                   corstr_sev, ".RData")

save(stdcoefs_BS_list, file = filename_stdcoefs_BS_list)

###################################################################

#What are the missing b?

missing_b_list = list()

for(age in ages)
{
  ind_age = which(ages %in% age)
  stdcoefs_BS_mat = stdcoefs_BS_list[[ind_age]]
  
  print(age)
  #these are the b that are missing:
  missing_b_list[[ind_age]] = which(is.na(colSums(stdcoefs_BS_mat)))
  names(missing_b_list)[ind_age] = paste0("age", age)
}

#Display the b's by age group
missing_b_list

#How many missing b per age group?
lapply(missing_b_list, length)

###################################################################

#Function to compute CI
get_CI = function(confidence_level=0.95, BS_dat)
{
  lq = (1-confidence_level)/2
  uq = 1 - lq
  
  d = data.frame(Variable=rownames(std_coef_sev)[-c(1,2)], 
                 LCL=apply(BS_dat, MARGIN = 1, FUN = function(x){quantile(x, lq, na.rm=T)}),
                 UCL=apply(BS_dat, MARGIN = 1, FUN = function(x){quantile(x, uq, na.rm=T)}))
  
  return(d)
}

preJS_CI_90=list()
preJS_CI_95=list()

#Get pre-JS CIs... 
for(age in ages)
{
  ind_age = which(ages %in% age)
  stdcoefs_BS_mat = stdcoefs_BS_list[[ind_age]]
  
  #95% bootstrap based quantiles:
  preJS_CI_95[[ind_age]] = get_CI(0.95, BS_dat = stdcoefs_BS_mat)
  
  #90% bootstrap based quantiles:
  preJS_CI_90[[ind_age]] = get_CI(0.90, BS_dat = stdcoefs_BS_mat)
  
  names(preJS_CI_95)[ind_age] =
    names(preJS_CI_90)[ind_age] = paste0("age", age)
}

# Write the CI results:

library(openxlsx)

# Create a new workbook for 95% CI
wb_95 <- createWorkbook()

# Loop through the list of dataframes for 95% CI
for (name in names(preJS_CI_95)) {
  addWorksheet(wb_95, name)  # Add a new worksheet with the name of the dataframe
  writeData(wb_95, sheet = name, preJS_CI_95[[name]])  # Write the dataframe to the worksheet
}

# Save the workbook for 95% CI
saveWorkbook(wb_95, paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/95_CI_",
                           corstr_sev, ".xlsx"), overwrite = TRUE)

# Create a new workbook for 90% CI
wb_90 <- createWorkbook()

# Loop through the list of dataframes for 90% CI
for (name in names(preJS_CI_90)) {
  addWorksheet(wb_90, name)  # Add a new worksheet with the name of the dataframe
  writeData(wb_90, sheet = name, preJS_CI_90[[name]])  # Write the dataframe to the worksheet
}

# Save the workbook for 90% CI
saveWorkbook(wb_90, paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/90_CI_",
                           corstr_sev, ".xlsx"), overwrite = TRUE)
