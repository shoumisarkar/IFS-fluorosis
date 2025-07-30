
#model_type = "03_presence_modelling"
#model_type = "04_severity_modelling"
#model_type = "05a_combined_presence_modelling"
model_type = "05b_combined_severity_modelling"

type = "pres"

if(grepl("severity", model_type, fixed = T))
{
  type = "sev"
}


library(openxlsx)

# Check if running in a SLURM environment
if (!is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))) {
  # If in SLURM environment
  setwd("/inSLURM/path/to/Fluorosis/")
} else {
  # If not in SLURM environment
  setwd("/nonSLURM/path/to/Fluorosis/")
}

ages = c(9, 13, 17, 23)

rho_df = data.frame(matrix(numeric(8*3),8, 3))
colnames(rho_df) = c("corstr", "age", "rho")
corstrs = c("exchangeable", "ar1")

count = 1

for(corstr in corstrs)
{
  for(age in ages)
  {
    corstr2 = corstr
    
    if(grepl("combined", model_type, fixed = T))
    {
      corstr2 = paste0(corstr, ",", corstr)
    }
    
    filepath = paste0("path/to/Fluorosis/Results/", model_type, "/",
                      corstr2, "/whole_data_based/rho_", type, "_age", age, ".RData")
    
    obj_rho <- load(filepath)
    assign("rho", get(obj_rho))
    
    rho_df[count,1] = corstr2
    rho_df[count,2] = age
    rho_df[count,3] = round(rho, digits = 4)
    
    count = count+1
  }
}

#View(rho_df)

filepath2 = paste0("path/to/Fluorosis/Results/",
                   model_type, "/summary_rho.xlsx")
write_xlsx(rho_df, path = filepath2)

