library(RankAggreg)
library(readxl)
library(dplyr)

ages = c(9, 13, 17, 23)
corstrs = c("independence", "exchangeable", "ar1", "jackknifed")
corstrs = paste(corstrs, corstrs, sep=",")

model_names = paste("C", 1:4, 1:4, sep=".")

#Get variable names
temp = read_xlsx(paste0("path/to/Fluorosis/Results/05b_combined_severity_modelling/95_table_",
                        "ar1,ar1", ".xlsx"), sheet = 1)
vars = temp$Variable

age_wise_mat = matrix(, nrow=4, ncol=4)
rownames(age_wise_mat) = paste0("Age ",ages)

for(age in ages)
{
  age_ind = which(ages %in% age)
  
  age_spec_mat = matrix(, nrow = length(vars), ncol = 4)
  rownames(age_spec_mat) = vars
  #colnames(age_spec_mat) = corstrs
  
  for(corstr in corstrs)
  {
    temp = read_xlsx(paste0("path/to/Fluorosis/Results/05b_combined_severity_modelling/95_table_",
                            corstr, ".xlsx"), sheet = age_ind)
    
    corstr_ind = which(corstrs %in% corstr)
    
    age_spec_mat[, corstr_ind] = as.numeric(temp$`James-Stein MSE`)
  }
  
  #age_spec_mat = t(apply(age_spec_mat, MARGIN = 1, function(x){corstrs[order(x)]}))
  age_spec_mat = t(apply(age_spec_mat, MARGIN = 1, function(x){model_names[order(x)]}))
  
  result <- RankAggreg(age_spec_mat, 4, method="CE", distance="Spearman", N=100, convIn=5, rho=.1)
  age_wise_mat[age_ind,] = result$top.list
}

overall_result <- RankAggreg(age_wise_mat, 4, method="CE", distance="Spearman", N=100, convIn=5, rho=.1)

all_outputs = rbind(age_wise_mat, Overall = overall_result$top.list); colnames(all_outputs) = paste0("Rank ", 1:4)
all_outputs

################################################
#### Write latex codes for Overleaf tables #####
################################################

library(xtable)

write_table_to_latex <- function(all_outputs) {
  # Create a .tex file to store the output
  tex_file <- paste0("path/to/Fluorosis/Results/05b_combined_severity_modelling/RankAggreg_combined_sev.tex")
  sink(tex_file)
  
  corstr_index = which(corstrs %in% corstr)
  
  #Model name: A.1.1 means A independence age
  caption_text <- paste0("Rank aggregation over different working correlation structures for the severity models from the combined modeling at different ages, and overall across all ages.")
  
  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat("\\caption{", caption_text, "}\n")
  cat("\\label{ch4:table:comb:sev:", "rankAggreg", "}\n", sep="")
  
  # Add the required LaTeX package
  #cat("\\usepackage{threeparttable}\n\n")
  
  
  #cat("\\scalebox{0.65}{\n")
  #cat("\\begin{tabular}{rlcccccl}\n")
  cat("\\begin{tabular}{rccccccl}\n")
  #cat("\\hline\n")
  
  #cat("Variable & Estimate & SE & {95\\% CI} \\\\\n")
  cat("Age & Rank 1 & Rank 2 & Rank 3 & Rank 4 \\\\\n")
  
  # Print the table body without column names
  print(xtable(all_outputs[1:5,], align = c("r", "c", "c", "c", "c")), 
        #print(xtable(my_list[[name]], align = c("r", "l", "c", "c", "c", "c", "c", "l")), 
        #print(xtable(my_list[[name]], align = c("r", "l", "c", "c", "c")), 
        include.rownames = T, include.colnames = FALSE, type = "latex", only.contents = TRUE)
  
  cat("\\end{tabular}\n")
  #cat("}\n")
  
  
  # Add custom notes at the end of the table
  cat("\\begin{tablenotes}\n")
  cat("\\scriptsize\n")
  cat("\\item C.$c_p$.$c_s$ denotes the model. ``C'' indicates the combined model, $c_p$ and $c_s$ denote the working correlation structures for \\\\presence and severity components respectively: $1$ for independence, $2$ for exchangeable, $3$ for AR(1), $4$ for jackknifed. \\\\Lower ranks correspond to lower SE of the James-Stein adjusted standardized estimate.\n")
  cat("\\end{tablenotes}\n")
  
  cat("\\end{table}\n\n")
  
  # Stop writing to the .tex file
  sink()
}


# Write the table to LaTeX
write_table_to_latex(all_outputs)


