

corstr = "ar1"
B=100

path_pres = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/", corstr, "/bootstrapping")
path_sev = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/", corstr, "/bootstrapping")

age_folders = paste0("age", c(9,13,17,23))

load_pres_obj = load(paste0(path_pres,
                            "/list_std_BS_estimators_pres_aggregated.RData"))
assign("df_BS_pres", get(load_pres_obj))

load_sev_obj = load(paste0(path_sev,
                            "/list_std_BS_estimators_sev_aggregated.RData"))
assign("df_BS_sev", get(load_sev_obj))


#Initialize data frame for gamma
df_gamma = data.frame(matrix(NA, ncol=4, nrow=B))
colnames(df_gamma) = age_folders

for(age_folder in age_folders)
{
  for(b in 1:B)
  {
    age_folder_ind = which(age_folders %in% age_folder)
  
    pres_b = df_BS_pres[[age_folder_ind]][,b]
    sev_b = df_BS_sev[[age_folder_ind]][,b]
    
    if( any(is.na(pres_b), is.na(sev_b)) ) #if there are any missing values in b
    {
      next
    }
    
    gamma_b = mean(sev_b[-c(1:2)]/pres_b[-c(1)])
    
    df_gamma[b, age_folder_ind] = gamma_b
  }
}

path_comb_a = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/05a_combined_presence_modelling/", corstr, "/bootstrapping")
path_comb_b = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/05b_combined_severity_modelling/", corstr, "/bootstrapping")

save(df_gamma, file=paste0(path_comb_a, "/pre_est_gamma_BS.RData"))
save(df_gamma, file=paste0(path_comb_b, "/pre_est_gamma_BS.RData"))
