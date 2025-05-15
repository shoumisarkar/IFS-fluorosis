#Update to higher B are we obtain more b (bootstrapped observations)

corstr = "ar1"
modelling_type = "04_severity_modelling" ####
B=400


path = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/", modelling_type, 
              "/", corstr, "/bootstrapping")

all_parameters = c("gamma", "alpha", "alpha_1|2", "alpha_2|3", "dental_age" ,  
                   "Total_mgF" , "SugarAddedBeverageOzPerDay" , "BrushingFrequencyPerDay",
                   "Avg_homeppm" , "Prop_DentAppt" , "Prop_FluorideTreatment" ,
                   "Tooth8" , "Tooth9" , "Tooth10" , "ZoneM" , "ZoneI" , "ZoneO")

sev_parameters = all_parameters[-c(1,2)]

age_folders = suppressWarnings(paste0("age", c(9, 13, 17, 23)))

#Initialize list to store severity BS results in:
df_BS_sev = data.frame(matrix(NA, nrow = length((sev_parameters)), ncol = B))
rownames(df_BS_sev) = sev_parameters

list_df_BS_sev = list(df_BS_sev, df_BS_sev, df_BS_sev, df_BS_sev)
names(list_df_BS_sev) = age_folders

age_folders = supsevsWarnings(paste0("age", c(9,13,17,23)))

for(age_folder in age_folders)
{
  filepath = paste0(path, "/", age_folder); setwd(filepath)
  
  for(b in 1:B)
  {
    filename_sev = paste0("BS_std_estimators_sev_", age_folder, "_b", b, ".RData") ####
    
    if(!file.exists(filename_sev))
    {
      next  
    }
    
    loaded_obj = load(filename_sev)
    assign("BS_std_estimators_sev", get(loaded_obj))
    
    age_folder_ind = which(age_folders %in% age_folder)
    BS_std_estimators_sev = as.vector(BS_std_estimators_sev)
    
    list_df_BS_sev[[age_folder_ind]][,b] = BS_std_estimators_sev
    
  }
  
}


# Tells which b indices are missing...
missing_b_indices_sev <- list(age9 = which(sapply(list_df_BS_sev$age9, function(x) any(is.na(x)))),
                               age13 = which(sapply(list_df_BS_sev$age13, function(x) any(is.na(x)))),
                               age17 = which(sapply(list_df_BS_sev$age17, function(x) any(is.na(x)))),
                               age23 = which(sapply(list_df_BS_sev$age23, function(x) any(is.na(x)))))

missing_b_indices_sev

# Save the workbook to a file
setwd(path)
save(list_df_BS_sev, file = "list_std_BS_estimators_sev_aggregated.RData")
save(missing_b_indices_sev, file = "missing_b_indices_sev.RData")
