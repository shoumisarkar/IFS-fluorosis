library(openxlsx)

# Check if running in a SLURM environment
if (!is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))) {
  # If in SLURM environment
  setwd("/blue/somnath.datta/shoumisarkar/Fluorosis/")
} else {
  # If not in SLURM environment
  setwd("W:/somnath.datta/shoumisarkar/Fluorosis/")
}

source(file = "Codes/functions.R")


# Helper function to compute CI
get_CI <- function(confidence_level = 0.95, BS_dat, vars) {
  lq <- (1 - confidence_level) / 2
  uq <- 1 - lq
  
  data.frame(
    Variable = vars,
    LCL = apply(BS_dat, MARGIN = 1, FUN = function(x) quantile(x, lq, na.rm = TRUE)),
    UCL = apply(BS_dat, MARGIN = 1, FUN = function(x) quantile(x, uq, na.rm = TRUE))
  )
}


# Helper function to write CI results to Excel
write_table_to_excel <- function(my_list, corstr_pres, conf_level) {
  wb <- createWorkbook()
  
  for (name in names(my_list)) {
    addWorksheet(wb, name)
    writeData(wb, sheet = name, my_list[[name]])
  }
  
  saveWorkbook(wb, paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/",
                          conf_level, "_table_", corstr_pres, ".xlsx"), overwrite = TRUE)
}

# Function to perform bootstrapping and save CI results
process_age_data <- function(age, B, corstr_pres, get_CI, vars) {
  stdcoefs_BS_mat <- matrix(NA, nrow = 13, ncol = B)
  path <- paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/",
                 corstr_pres, "/bootstrapping/age", age, "/")
  
  for (b in 1:B) {
    fn_stdcoef_pres <- paste0(path, "std_coefs_pres_age", age, "_b", b, ".RData")
    if (file.exists(fn_stdcoef_pres)) {
      obj_stdcoef_pres <- load(fn_stdcoef_pres)
      assign("std_coef_pres", get(obj_stdcoef_pres))
      stdcoefs_BS_mat[, b] <- std_coef_pres$`Standardized Estimate`[-1]  # Adjusted to remove first row
    }
  }
  
  return(stdcoefs_BS_mat)
}

# Main logic
corstr_pres <- "exchangeable"
ages <- c(9, 13, 17, 23)
B <- 120
vars <- c("dental_age", "Total_mgF", "SugarAddedBeverageOzPerDay", "BrushingFrequencyPerDay",
          "Avg_homeppm", "Prop_DentAppt", "Prop_FluorideTreatment", "Tooth8", "Tooth9",
          "Tooth10", "ZoneM", "ZoneI", "ZoneO")

stdcoefs_BS_list <- lapply(ages, process_age_data, B = B, corstr_pres = corstr_pres, get_CI = get_CI, vars = vars)
names(stdcoefs_BS_list) <- paste0("age", ages)


# Calculate CIs
preJS_CI_95 <- lapply(stdcoefs_BS_list, get_CI, confidence_level = 0.95, vars = vars)
preJS_CI_90 <- lapply(stdcoefs_BS_list, get_CI, confidence_level = 0.90, vars = vars)


# Repeat the process for post-James-Stein
stdcoefsJS_BS_list <- stdcoefs_BS_list  # will be updated below, following the initialized structure

#vars = rownames(std_coef_pres)[-c(1)]

for(var in vars)
{
  var_ind = which(vars %in% var)
  
  for(b in 1:B)
  {
    t <- c(
      (stdcoefs_BS_list[[1]])[var_ind, b],
      (stdcoefs_BS_list[[2]])[var_ind, b],
      (stdcoefs_BS_list[[3]])[var_ind, b],
      (stdcoefs_BS_list[[4]])[var_ind, b]
    )
    
    t = js(t)
    
    stdcoefsJS_BS_list[[1]][var_ind,b] = t[1]
    stdcoefsJS_BS_list[[2]][var_ind,b] = t[2]
    stdcoefsJS_BS_list[[3]][var_ind,b] = t[3]
    stdcoefsJS_BS_list[[4]][var_ind,b] = t[4]
    
  }
}


postJS_CI_95 <- lapply(stdcoefsJS_BS_list, get_CI, confidence_level = 0.95, vars = vars)
postJS_CI_90 <- lapply(stdcoefsJS_BS_list, get_CI, confidence_level = 0.90, vars = vars)


#Obtain the whole data based estimates and SE:

coef_pres_list = list()
se_pres_list = list()

for(age in ages)
{
  path1 = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/", corstr_pres, "/whole_data_based/coefs_pres_age", age, ".RData")
  obj_coef_pres <- load(path1)
  assign("coef_pres", get(obj_coef_pres))
  
  
  path2 = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/", corstr_pres, "/whole_data_based/JK_SD_pres_age", age, ".RData")
  obj_se_pres <- load(path2)
  assign("se_pres", get(obj_se_pres))
  
  
  age_ind = which(ages %in% age)
  
  coef_pres_list[[age_ind]] = coef_pres[-c(1),]
  se_pres_list[[age_ind]] = se_pres[-c(1), , drop=F]
}

names(coef_pres_list) = names(se_pres_list) = paste0("age", ages)


#Get JS estimates based on these whole data based estimates

postJS_stdcoef_pres_list = coef_pres_list #initialize

for(var in vars)
{
  t = c()
  var_ind = which(vars %in% var)
  
  for(age in ages)
  {
    age_ind = which(ages %in% age)
    
    stdcoef = (coef_pres_list[[age_ind]])$Estimates[var_ind]/
      (se_pres_list[[age_ind]])$SD[var_ind]
    
    t = c(t, stdcoef)
  }
  
  js_output = js(t) #James-Stein updates for current variable across the four timepoints
  
  for(age in ages)
  {
    age_ind = which(ages %in% age)
    
    postJS_stdcoef_pres_list[[age_ind]]$Estimates[var_ind] = js_output[age_ind] #update in the JS estimates
  }
  
}


#Create the summary table with 95% CI:

create_summary_table = function(level=95)
{
  
  preJS_CI = preJS_CI_90
  postJS_CI = postJS_CI_90
  
  if(level==95)
  {
    preJS_CI = preJS_CI_95
    postJS_CI = postJS_CI_95
  }
  
  summary_pres = list()
  
  for(age in ages)
  {
    age_ind = which(ages %in% age)
    
    preJS_CI_int =
      paste0( "(",
              format(round(preJS_CI[[age_ind]]$LCL,3), nsmall=3), ", ",
              format(round(preJS_CI[[age_ind]]$UCL,3), nsmall=3),
              ")"
      )
    
    # Determine significance and direction symbols
    ind_preJS_signif <- preJS_CI[[age_ind]]$LCL * preJS_CI[[age_ind]]$UCL > 0
    preJS_signif <- ifelse(ind_preJS_signif, "*", "")
    preJS_direction <- ifelse(coef_pres_list[[age_ind]]$Estimates > 0, "+", "-")
    preJS_symbol <- ifelse(ind_preJS_signif, paste0("$^{", preJS_signif, preJS_direction, "}$"), "")
    preJS_CI_int <- paste0(preJS_CI_int, preJS_symbol)
    
    
    
    postJS_CI_int =
      paste0( "(",
              format(round(postJS_CI[[age_ind]]$LCL,3), nsmall=3), ", ",
              format(round(postJS_CI[[age_ind]]$UCL,3), nsmall=3),
              ")"
      )
    
    
    ind_postJS_signif <- postJS_CI[[age_ind]]$LCL * postJS_CI[[age_ind]]$UCL > 0
    postJS_signif <- ifelse(ind_postJS_signif, "*", "")
    postJS_direction <- ifelse(coef_pres_list[[age_ind]]$Estimates > 0, "+", "-")
    postJS_symbol <- ifelse(ind_postJS_signif, paste0("$^{", postJS_signif, postJS_direction, "}$"), "")
    postJS_CI_int <- paste0(postJS_CI_int, postJS_symbol)
    
    temp_df  = data.frame(Variable=rownames(se_pres_list$age9),
                          Estimate=format(round((coef_pres_list[[age_ind]])$Estimates, 3), nsmall=3),
                          SE=format(round((se_pres_list[[age_ind]])$SD, 3), nsmall=3), 
                          StdEstimate=round((coef_pres_list[[age_ind]])$Estimates/se_pres_list[[age_ind]]$SD, 3),
                          StdCI=preJS_CI_int,
                          JS_StdEstimate = format(round(postJS_stdcoef_pres_list[[age_ind]]$Estimates, 3), nsmall=3),
                          JS_StdCI=postJS_CI_int)
    
    colnames(temp_df) = c("Variable", "Estimate", "Standard Error", "Standarized Estimate",
                          "Standardized CI", "James-Stein Standardized Estimate", "James-Stein Standardized CI")
    
    summary_pres[[age_ind]] = temp_df
  }
  
  names(summary_pres) = paste0("age", ages)
  
  return(summary_pres)
}



# Create summary tables with 95% CI and 90% CI
summary_95 <- create_summary_table(95)
summary_90 <- create_summary_table(90)


write_table_to_excel(summary_95, corstr_pres, 95)
write_table_to_excel(summary_90, corstr_pres, 90)


#Write latex tables:

library(xtable)

write_table_to_latex <- function(my_list, corstr_pres, conf_level) {
  # Create a .tex file to store the output
  tex_file <- paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/03_presence_modelling/", conf_level, "_table_", corstr_pres, ".tex")
  sink(tex_file)
  
  corstr_index = which(c("independence", "exchangeable", "ar1", "jackknifed") %in% corstr_pres)
  
  #Model name: A.1.1 means A independence age
  caption_text <- paste0("Estimates from models A.", corstr_index, ".1-A.", corstr_index, ".4, the separate presence models with the ",
                         corstr_pres, " cluster correlation structure")
  
  # gsub("age", "", name)
  
  #cat("\\section*{", name, "}\n")
  cat("\\begin{table}[ht]\n")
  cat("\\label{table:sep:pres:", corstr_pres, "}\n", sep="")
  cat("\\centering\n")
  cat("\\caption{", caption_text, "}\n")
  cat("\\begin{threeparttable}\n")
  cat("\\centering\n")
  
  # Add the required LaTeX package
  #cat("\\usepackage{threeparttable}\n\n")
  
  # Loop through the list of data frames (tables)
  for (name in names(my_list)) {
      
    age = gsub("age", "", name); age_index = which(ages %in% age)
    
    
    cat("%Table for", name,"\n")
    
    
    cat("\\begin{subtable}{\\linewidth}\n")
    cat("\\centering\n")
    cat("\\caption{Model A.", corstr_index, ".", age_index, " (age ",age,")}\n", sep="")
    cat("\\scalebox{0.85}{\n")
    cat("\\begin{tabular}{rlcccccl}\n")
    
    #cat("\\hline\n")
    cat("Variable & Estimate & SE & \\makecell{Standardized\\\\Estimate} & 95\\% CI & \\makecell{James-Stein\\\\Estimate} & \\makecell{95\\% CI\\\\(James-Stein)} \\\\\n")
    
    # Print the table body without column names
    print(xtable(my_list[[name]], align = c("r", "l", "c", "c", "c", "c", "c", "l")), 
          include.rownames = FALSE, include.colnames = FALSE, type = "latex", only.contents = TRUE)
    
    cat("\\end{tabular}\n")
    cat("}\n")
    cat("\\end{subtable}\n")
    
  }
  
  # Add custom notes at the end of the table
  cat("\\begin{tablenotes}[para,flushleft]\n")
  cat("\\scriptsize\n")
  cat("\\item Superscripts $*+$ and $*-$ denote significant positive and negative effects at the 5\\% significance level, respectively.\n")
  cat("\\end{tablenotes}\n")
  
  cat("\\end{threeparttable}\n")
  cat("\\end{table}\n\n")
  
  # Stop writing to the .tex file
  sink()
}


# Write the 95% CI and 90% CI tables to LaTeX
write_table_to_latex(my_list = summary_95, corstr_pres = corstr_pres, conf_level = 95)
write_table_to_latex(summary_90, corstr_pres, 90)

