# Subtable with the raw estimates

#############################################
### Choose correlation structure and age: ###
#############################################

corstr_pres = "jackknifed"
corstr_sev = "jackknifed"

# corstr_pres = "exchangeable"
# corstr_sev = "exchangeable"

# corstr_pres = "independence"
# corstr_sev = "exchangeable"

# corstr_pres = "exchangeable"
# corstr_sev = "independence"

corstrs = c("independence", "exchangeable", "jackknifed")
corstr_index_pres = which(corstrs %in% corstr_pres)
corstr_index_sev = which(corstrs %in% corstr_sev)


ages = c(9,13,17,23)

for(age in ages)
{
  ##############################
  ### Other initializations: ###
  ##############################
  
  age_ind = which(ages %in% age)
  
  library(xtable)
  
  # Create a .tex file to store the output
  tex_file <- paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/raw_tables_", corstr_pres,
                     ",", corstr_sev, "_age", age, ".tex")
  sink(tex_file)
  
  #########################
  ### PRESENCE subtable ###
  #########################
  
  #load the presence RData file
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N30/presence/summarytable_raw_est_presence_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table30", get(obj))
  
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N50/presence/summarytable_raw_est_presence_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table50", get(obj))
  
  
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N200/presence/summarytable_raw_est_presence_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table200", get(obj))
  
  table = table30
  
  for(i in ages)
  {
    table[[which(ages %in% age)]] = cbind(table30[[which(ages %in% age)]], table50[[which(ages %in% age)]][,-c(1)])  
    table[[which(ages %in% age)]] = cbind(table[[which(ages %in% age)]], 
                                          table200[[which(ages %in% age)]][,-c(1)])
  }
  
  caption_text <- paste0("Simulation results for age ", age, " data generated with ", corstr_pres, " and ", 
                         corstr_sev, " correlation structures for presence and severity respectively.")
  
  cat("\\begin{table}[ht]\n")
  cat("\\label{ch4:table:sim:raw:", corstr_pres, "}\n", sep="")
  cat("\\centering\n")
  cat("\\caption{", caption_text, "}\n")
  cat("\\begin{threeparttable}\n")
  cat("\\centering\n")
  
  # Add the required LaTeX package
  #cat("\\usepackage{threeparttable}\n\n")
  
  # Loop through the list of data frames (tables)
  
  cat("%Table for presence \n")
  
  cat("\\begin{subtable}{\\linewidth}\n")
  cat("\\centering\n")
  cat("\\caption{Presence model A$_s$.", corstr_index_pres, ".", age_ind, "}\n", sep="")
  cat("\\scalebox{0.65}{\n")
  cat("\\begin{tabular}{l|cccc|cccc|cccc}\n")
  
  cat("{} & \\multicolumn{4}{c|}{\\textbf{N=30}} & \\multicolumn{4}{c|}{\\textbf{N=50}} & \\multicolumn{4}{c}{\\textbf{N=200}} \\\\ \\hline\n")
  
  cat("Variable & Estimate & Bias & SE & MSE & Estimate & Bias & SE & MSE & Estimate & Bias & SE & MSE \\\\ \n")
  
  # Print the table body without column names
  print(xtable(table[[age_ind]], align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c")), 
        include.rownames = FALSE, include.colnames = FALSE, type = "latex", only.contents = TRUE)
  
  cat("\\end{tabular}\n")
  cat("}\n")
  cat("\\end{subtable}\n")
  
  
  
  #########################
  ### SEVERITY subtable ###
  #########################
  
  #load the presence RData file
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N30/severity/summarytable_raw_est_severity_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table30", get(obj))
  
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N50/severity/summarytable_raw_est_severity_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table50", get(obj))
  
  
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N200/severity/summarytable_raw_est_severity_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table200", get(obj))
  
  table = table30
  
  for(i in ages)
  {
    table[[which(ages %in% age)]] = cbind(table30[[which(ages %in% age)]], table50[[which(ages %in% age)]][,-c(1)])  
    table[[which(ages %in% age)]] = cbind(table[[which(ages %in% age)]], 
                                          table200[[which(ages %in% age)]][,-c(1)])
  }
  
  
  cat("%Table for severity \n")
  
  cat("\\begin{subtable}{\\linewidth}\n")
  cat("\\centering\n")
  cat("\\caption{Severity model B$_s$.", corstr_index_sev, ".", age_ind, "}\n", sep="")
  cat("\\scalebox{0.65}{\n")
  cat("\\begin{tabular}{l|cccc|cccc|cccc}\n")
  
  cat("{} & \\multicolumn{4}{c|}{\\textbf{N=30}} & \\multicolumn{4}{c|}{\\textbf{N=50}} & \\multicolumn{4}{c}{\\textbf{N=200}} \\\\ \\hline\n")
  
  cat("Variable & Estimate & Bias & SE & MSE & Estimate & Bias & SE & MSE & Estimate & Bias & SE & MSE \\\\ \n")
  
  # Print the table body without column names
  print(xtable(table[[age_ind]], align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c")), 
        include.rownames = FALSE, include.colnames = FALSE, type = "latex", only.contents = TRUE)
  
  cat("\\end{tabular}\n")
  cat("}\n")
  cat("\\end{subtable}\n")
  
  
  
  ##################################
  ### combined PRESENCE subtable ###
  ##################################
  
  #load the presence RData file
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N30/combined/combined_presence_summarytable_raw_est_combined_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table30", get(obj))
  
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N50/combined/combined_presence_summarytable_raw_est_combined_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table50", get(obj))
  
  
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N200/combined/combined_presence_summarytable_raw_est_combined_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table200", get(obj))
  
  table = table30
  
  for(i in ages)
  {
    table[[which(ages %in% age)]] = cbind(table30[[which(ages %in% age)]], table50[[which(ages %in% age)]][,-c(1)])  
    table[[which(ages %in% age)]] = cbind(table[[which(ages %in% age)]], 
                                          table200[[which(ages %in% age)]][,-c(1)])
  }
  
  
  cat("%Table for severity \n")
  
  cat("\\begin{subtable}{\\linewidth}\n")
  cat("\\centering\n")
  cat("\\caption{Presence piece of model C$_s$.", corstr_index_pres, ".", corstr_index_sev, ".", age_ind, "}\n", sep="")
  cat("\\scalebox{0.65}{\n")
  cat("\\begin{tabular}{l|cccc|cccc|cccc}\n")
  
  cat("{} & \\multicolumn{4}{c|}{\\textbf{N=30}} & \\multicolumn{4}{c|}{\\textbf{N=50}} & \\multicolumn{4}{c}{\\textbf{N=200}} \\\\ \\hline\n")
  
  cat("Variable & Estimate & Bias & SE & MSE & Estimate & Bias & SE & MSE & Estimate & Bias & SE & MSE \\\\ \n")
  
  # Print the table body without column names
  print(xtable(table[[age_ind]], align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c")), 
        include.rownames = FALSE, include.colnames = FALSE, type = "latex", only.contents = TRUE)
  
  cat("\\end{tabular}\n")
  cat("}\n")
  cat("\\end{subtable}\n")
  
  
  ##################################
  ### combined SEVERITY subtable ###
  ##################################
  
  #load the presence RData file
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N30/combined/combined_severity_summarytable_raw_est_combined_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table30", get(obj))
  
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N50/combined/combined_severity_summarytable_raw_est_combined_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table50", get(obj))
  
  
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N200/combined/combined_severity_summarytable_raw_est_combined_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table200", get(obj))
  
  table = table30
  
  for(i in ages)
  {
    table[[which(ages %in% age)]] = cbind(table30[[which(ages %in% age)]], table50[[which(ages %in% age)]][,-c(1)])  
    table[[which(ages %in% age)]] = cbind(table[[which(ages %in% age)]], 
                                          table200[[which(ages %in% age)]][,-c(1)])
  }
  
  
  cat("%Table for severity \n")
  
  cat("\\begin{subtable}{\\linewidth}\n")
  cat("\\centering\n")
  cat("\\caption{Severity piece of model C$_s$.", corstr_index_pres, ".", corstr_index_sev, ".", age_ind, "}\n", sep="")
  cat("\\scalebox{0.65}{\n")
  cat("\\begin{tabular}{l|cccc|cccc|cccc}\n")
  
  cat("{} & \\multicolumn{4}{c|}{\\textbf{N=30}} & \\multicolumn{4}{c|}{\\textbf{N=50}} & \\multicolumn{4}{c}{\\textbf{N=200}} \\\\ \\hline\n")
  
  cat("Variable & Estimate & Bias & SE & MSE & Estimate & Bias & SE & MSE & Estimate & Bias & SE & MSE \\\\ \n")
  
  # Print the table body without column names
  print(xtable(table[[age_ind]], align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c")), 
        include.rownames = FALSE, include.colnames = FALSE, type = "latex", only.contents = TRUE)
  
  cat("\\end{tabular}\n")
  cat("}\n")
  cat("\\end{subtable}\n")
  
  cat("\\end{threeparttable}\n")
  cat("\\end{table}\n\n")
  
  # Stop writing to the .tex file
  sink()
}

  