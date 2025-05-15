#Update to higher B are we obtain more b (bootstrapped observations)

corstr = "ar1"
modelling_type = "03_presence_modelling" ####
B=400


path = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/", modelling_type, 
              "/", corstr, "/bootstrapping")

all_parameters = c("gamma", "alpha", "alpha_1|2", "alpha_2|3", "dental_age" ,  
                   "Total_mgF" , "SugarAddedBeverageOzPerDay" , "BrushingFrequencyPerDay",
                   "Avg_homeppm" , "Prop_DentAppt" , "Prop_FluorideTreatment" ,
                   "Tooth8" , "Tooth9" , "Tooth10" , "ZoneM" , "ZoneI" , "ZoneO")

pres_parameters = all_parameters[-c(1,3,4)]

age_folders = suppressWarnings(paste0("age", c(9, 13, 17, 23)))

#Initialize list to store presence BS results in:
df_BS_pres = data.frame(matrix(NA, nrow = length((pres_parameters)), ncol = B))
rownames(df_BS_pres) = pres_parameters

list_df_BS_pres = list(df_BS_pres, df_BS_pres, df_BS_pres, df_BS_pres)
names(list_df_BS_pres) = age_folders

age_folders = suppressWarnings(paste0("age", c(9,13,17,23)))

for(age_folder in age_folders)
{
  filepath = paste0(path, "/", age_folder); setwd(filepath)
  
  for(b in 1:B)
  {
    filename_pres = paste0("BS_std_estimators_pres_", age_folder, "_b", b, ".RData") ####
    
    if(!file.exists(filename_pres))
    {
        next  
    }
    
    loaded_obj = load(filename_pres)
    assign("BS_std_estimators_pres", get(loaded_obj))
    
    age_folder_ind = which(age_folders %in% age_folder)
    BS_std_estimators_pres = as.vector(BS_std_estimators_pres)
    
    list_df_BS_pres[[age_folder_ind]][,b] = BS_std_estimators_pres
    
  }
  
}


# Tells which b indices are missing...
missing_b_indices_pres <- list(age9 = which(sapply(list_df_BS_pres$age9, function(x) any(is.na(x)))),
                                age13 = which(sapply(list_df_BS_pres$age13, function(x) any(is.na(x)))),
                                age17 = which(sapply(list_df_BS_pres$age17, function(x) any(is.na(x)))),
                                age23 = which(sapply(list_df_BS_pres$age23, function(x) any(is.na(x)))))

missing_b_indices_pres


# Save the workbook to a file
setwd(path)
save(list_df_BS_pres, file = "list_std_BS_estimators_pres_aggregated.RData")
save(missing_b_indices_pres, file = "missing_b_indices_pres.RData")
