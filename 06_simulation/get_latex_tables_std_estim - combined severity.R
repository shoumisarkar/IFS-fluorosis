##############################
#### COMBINED SEVERITY #######
##############################

# Subtable with the std estimates

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
  tex_file <- paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/std_tables_combined_severity_", corstr_pres,
                     ",", corstr_sev, "_age", age, ".tex")
  sink(tex_file)
  
  #load the severity RData file
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N30/combined/summarytable_std_estimates_combined_severity_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table30", get(obj))
  
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N50/combined/summarytable_std_estimates_combined_severity_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table50", get(obj))
  
  filename = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N200/combined/summarytable_std_estimates_combined_severity_",
                    corstr_pres, ",", corstr_sev, ".RData")
  obj = load(filename); assign("table200", get(obj))
  
  #per model per age group, have subtables based on cluster size N
  
  caption_text <- paste0("Simulation results for age ", age, " data generated with ", corstr_pres, " and ", 
                         corstr_sev, " correlation structures for presence and severity respectively -",
                         " properties of standardized estimates arising from the severity piece of the combined model C$_s$.", 
                         corstr_index_pres, ".", corstr_index_sev, ".", age_ind, ".")
  
  cat("\\begin{table}[ht]\n")
  cat("\\label{ch4:table:sim:std:comb:sev:", corstr_sev, "}\n", sep="")
  cat("\\centering\n")
  cat("\\caption{", caption_text, "}\n")
  cat("\\begin{threeparttable}\n")
  cat("\\centering\n")
  
  # Add the required LaTeX package
  #cat("\\usepackage{threeparttable}\n\n")
  
  
  #####################
  ### N=30 subtable ###
  #####################
  
  cat("%Table for N=30 \n")
  
  cat("\\begin{subtable}{\\linewidth}\n")
  cat("\\centering\n")
  cat("\\caption{N=30", "}\n", sep="")
  cat("\\scalebox{0.65}{\n")
  cat("\\begin{tabular}{lcccc|cccc}\n")
  
  #cat("{} & \\multicolumn{4}{c|}{\\textbf{N=30}} & \\multicolumn{4}{c|}{\\textbf{N=50}} & \\multicolumn{4}{c}{\\textbf{N=200}} \\\\ \\hline\n")
  
  cat("Variable & \\makecell{Standardized\\\\Estimate} & \\makecell{Bias\\\\(Standardized\\\\Estimate)} & \\makecell{SE\\\\(Standardized\\\\Estimate)} & \\makecell{MSE\\\\(Standardized\\\\Estimate)} & \\makecell{Standardized\\\\James-Stein\\\\Estimate} & \\makecell{Bias\\\\(Standardized\\\\James-Stein\\\\Estimate)} & \\makecell{SE\\\\(Standardized\\\\James-Stein\\\\Estimate)} & \\makecell{MSE\\\\(Standardized\\\\James-Stein\\\\Estimate)} \\\\ \n")
  
  # Print the table body without column names
  print(xtable(table30[[age_ind]], align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c")), 
        include.rownames = FALSE, include.colnames = FALSE, type = "latex", only.contents = TRUE)
  cat("\\end{tabular}\n")
  cat("}\n")
  cat("\\end{subtable}\n")
  
  
  
  
  
  #####################
  ### N=50 subtable ###
  #####################
  
  cat("%Table for N=50 \n")
  
  cat("\\begin{subtable}{\\linewidth}\n")
  cat("\\centering\n")
  cat("\\caption{N=50", "}\n", sep="")
  cat("\\scalebox{0.65}{\n")
  cat("\\begin{tabular}{lcccc|cccc}\n")
  
  #cat("{} & \\multicolumn{4}{c|}{\\textbf{N=30}} & \\multicolumn{4}{c|}{\\textbf{N=50}} & \\multicolumn{4}{c}{\\textbf{N=200}} \\\\ \\hline\n")
  
  cat("Variable & \\makecell{Standardized\\\\Estimate} & \\makecell{Bias\\\\(Standardized\\\\Estimate)} & \\makecell{SE\\\\(Standardized\\\\Estimate)} & \\makecell{MSE\\\\(Standardized\\\\Estimate)} & \\makecell{Standardized\\\\James-Stein\\\\Estimate} & \\makecell{Bias\\\\(Standardized\\\\James-Stein\\\\Estimate)} & \\makecell{SE\\\\(Standardized\\\\James-Stein\\\\Estimate)} & \\makecell{MSE\\\\(Standardized\\\\James-Stein\\\\Estimate)} \\\\ \n")
  
  # Print the table body without column names
  print(xtable(table50[[age_ind]], align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c")), 
        include.rownames = FALSE, include.colnames = FALSE, type = "latex", only.contents = TRUE)
  cat("\\end{tabular}\n")
  cat("}\n")
  cat("\\end{subtable}\n")
  
  
  
  
  #####################
  ### N=200 subtable ###
  #####################
  
  cat("%Table for N=200 \n")
  
  cat("\\begin{subtable}{\\linewidth}\n")
  cat("\\centering\n")
  cat("\\caption{N=200", "}\n", sep="")
  cat("\\scalebox{0.65}{\n")
  cat("\\begin{tabular}{lcccc|cccc}\n")
  
  #cat("{} & \\multicolumn{4}{c|}{\\textbf{N=30}} & \\multicolumn{4}{c|}{\\textbf{N=50}} & \\multicolumn{4}{c}{\\textbf{N=200}} \\\\ \\hline\n")
  
  cat("Variable & \\makecell{Standardized\\\\Estimate} & \\makecell{Bias\\\\(Standardized\\\\Estimate)} & \\makecell{SE\\\\(Standardized\\\\Estimate)} & \\makecell{MSE\\\\(Standardized\\\\Estimate)} & \\makecell{Standardized\\\\James-Stein\\\\Estimate} & \\makecell{Bias\\\\(Standardized\\\\James-Stein\\\\Estimate)} & \\makecell{SE\\\\(Standardized\\\\James-Stein\\\\Estimate)} & \\makecell{MSE\\\\(Standardized\\\\James-Stein\\\\Estimate)} \\\\ \n")
  
  # Print the table body without column names
  print(xtable(table200[[age_ind]], align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c")), 
        include.rownames = FALSE, include.colnames = FALSE, type = "latex", only.contents = TRUE)
  cat("\\end{tabular}\n")
  cat("}\n")
  cat("\\end{subtable}\n")
  

  
  
  
  cat("\\end{threeparttable}\n")
  cat("\\end{table}\n\n")
  
  # Stop writing to the .tex file
  sink()
}

  