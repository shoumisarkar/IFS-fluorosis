corstr = "ar1"

path_pres = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/", corstr, 
                   "/whole_data_based")

path_sev = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/", corstr,
                  "/whole_data_based")

age_folders = suppressWarnings(paste0("age", c(9,13,17,23)))

df_gamma = data.frame(matrix(nrow=1, ncol=4))
colnames(df_gamma) = age_folders

for(age_folder in age_folders)
{
  obj_pres = load(paste0(path_pres,"/coeffs_pres_", age_folder, ".RData"))
  assign("coef_pres", get(obj_pres))
  
  obj_sev = load(paste0(path_sev,"/coeffs_sev_", age_folder, ".RData"))
  assign("coef_sev", get(obj_sev))
  
  pre_est_gamma = mean(coef_sev$Estimates[-c(1:2)]/coef_pres$Estimates[-c(1)])
  
  age_folder_ind = which(age_folders %in% age_folder)

  df_gamma[1,age_folder_ind] = pre_est_gamma
    
}

output_path_a = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/05a_combined_presence_modelling/",
                      corstr, "/whole_data_based")
                      
output_path_b = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/05b_combined_severity_modelling/",
                       corstr, "/whole_data_based")

save(df_gamma, file=paste0(output_path_a, "/pre_est_gamma.RData"))
save(df_gamma, file=paste0(output_path_b, "/pre_est_gamma.RData"))
