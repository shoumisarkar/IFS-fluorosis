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
write_table_to_excel <- function(CI_list, corstr_sev, conf_level) {
  wb <- createWorkbook()
  
  for (name in names(CI_list)) {
    addWorksheet(wb, name)
    writeData(wb, sheet = name, CI_list[[name]])
  }
  
  saveWorkbook(wb, paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/",
                          conf_level, "_table_", corstr_sev, ".xlsx"), overwrite = TRUE)
}

# Function to perform bootstrapping and save CI results
process_age_data <- function(age, B, corstr_sev, get_CI, vars) {
  stdcoefs_BS_mat <- matrix(NA, nrow = 13, ncol = B)
  path <- paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/",
                 corstr_sev, "/bootstrapping/age", age, "/")
  
  for (b in 1:B) {
    fn_stdcoef_sev <- paste0(path, "std_coefs_sev_age", age, "_b", b, ".RData")
    if (file.exists(fn_stdcoef_sev)) {
      obj_stdcoef_sev <- load(fn_stdcoef_sev)
      assign("std_coef_sev", get(obj_stdcoef_sev))
      stdcoefs_BS_mat[, b] <- std_coef_sev$`Standardized Estimate`[-c(1:2)]  # Adjusted to remove first row
    }
  }
  
  return(stdcoefs_BS_mat)
}

# Main logic
corstr_sev <- "jackknifed"
ages <- c(9, 13, 17, 23)
B <- 120
vars <- c("dental_age", "Total_mgF", "SugarAddedBeverageOzPerDay", "BrushingFrequencyPerDay",
          "Avg_homeppm", "Prop_DentAppt", "Prop_FluorideTreatment", "Tooth8", "Tooth9",
          "Tooth10", "ZoneM", "ZoneI", "ZoneO")

stdcoefs_BS_list <- lapply(ages, process_age_data, B = B, corstr_sev = corstr_sev, get_CI = get_CI, vars = vars)
names(stdcoefs_BS_list) <- paste0("age", ages)


# Calculate CIs
preJS_CI_95 <- lapply(stdcoefs_BS_list, get_CI, confidence_level = 0.95, vars = vars)
preJS_CI_90 <- lapply(stdcoefs_BS_list, get_CI, confidence_level = 0.90, vars = vars)

# Repeat the process for post-James-Stein
stdcoefsJS_BS_list <- stdcoefs_BS_list  # will be updated below, following the initialized structure

#vars = rownames(std_coef_sev)[-c(1:2)]

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
    
    t=js(t)
    
    stdcoefsJS_BS_list[[1]][var_ind,b] = t[1]
    stdcoefsJS_BS_list[[2]][var_ind,b] = t[2]
    stdcoefsJS_BS_list[[3]][var_ind,b] = t[3]
    stdcoefsJS_BS_list[[4]][var_ind,b] = t[4]
    
  }
}


postJS_CI_95 <- lapply(stdcoefsJS_BS_list, get_CI, confidence_level = 0.95, vars = vars)
postJS_CI_90 <- lapply(stdcoefsJS_BS_list, get_CI, confidence_level = 0.90, vars = vars)


#Obtain the whole data based estimates and SE:

coef_sev_list = list()
se_sev_list = list()

for(age in ages)
{
  path1 = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/", corstr_sev, "/whole_data_based/coefs_sev_age", age, ".RData")
  obj_coef_sev <- load(path1)
  assign("coef_sev", get(obj_coef_sev))
  
  
  path2 = paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/", corstr_sev, "/whole_data_based/JK_SD_sev_age", age, ".RData")
  obj_se_sev <- load(path2)
  assign("se_sev", get(obj_se_sev))
  
  
  age_ind = which(ages %in% age)
  
  coef_sev_list[[age_ind]] = coef_sev[-c(1,2),]
  se_sev_list[[age_ind]] = se_sev[-c(1,2), , drop=F]
}

names(coef_sev_list) = names(se_sev_list) = paste0("age", ages)


#Get JS estimates based on these whole data based estimates

postJS_stdcoef_sev_list = coef_sev_list #initialize

for(var in vars)
{
  t = c()
  var_ind = which(vars %in% var)
  
  for(age in ages)
  {
    age_ind = which(ages %in% age)
    
    stdcoef = (coef_sev_list[[age_ind]])$Estimates[var_ind]/
      (se_sev_list[[age_ind]])$SD[var_ind]
    
    t = c(t, stdcoef)
  }
  
  js_output = js(t) #James-Stein updates for current variable across the four timepoints
  
  for(age in ages)
  {
    age_ind = which(ages %in% age)
    
    postJS_stdcoef_sev_list[[age_ind]]$Estimates[var_ind] = js_output[age_ind] #update in the JS estimates
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
  
  summary_sev = list()
  
  for(age in ages)
  {
    age_ind = which(ages %in% age)
    
    
    preJS_CI_int =
      paste0( "(",
              format(round(preJS_CI[[age_ind]]$LCL,3), nsmall=3), ", ",
              format(round(preJS_CI[[age_ind]]$UCL,3), nsmall=3),
              ")"
      )
    
    
    postJS_CI_int =
      paste0( "(",
              format(round(postJS_CI[[age_ind]]$LCL,3), nsmall=3), ", ",
              format(round(postJS_CI[[age_ind]]$UCL,3), nsmall=3),
              ")"
      )
    
    
    
    temp_df  = data.frame(Variable=rownames(se_sev_list$age9),
                          Estimate=format(round((coef_sev_list[[age_ind]])$Estimates, 3), nsmall=3),
                          SE=format(round((se_sev_list[[age_ind]])$SD, 3), nsmall=3), 
                          StdEstimate=round((coef_sev_list[[age_ind]])$Estimates/se_sev_list[[age_ind]]$SD, 3),
                          StdCI=preJS_CI_int,
                          JS_StdEstimate = format(round(postJS_stdcoef_sev_list[[age_ind]]$Estimates, 3), nsmall=3),
                          JS_StdCI=postJS_CI_int)
    
    colnames(temp_df) = c("Variable", "Estimate", "Standard Error", "Standarized Estimate",
                          "Standardized CI", "James-Stein Standardized Estimate", "James-Stein Standardized CI")
    
    summary_sev[[age_ind]] = temp_df
  }
  
  names(summary_sev) = paste0("age", ages)
  
  return(summary_sev)
}

write_table_to_excel(create_summary_table(95), corstr_sev, 95)
write_table_to_excel(create_summary_table(90), corstr_sev, 90)

#Write latex tables:

library(xtable)

write_table_to_latex <- function(my_list, corstr_sev, conf_level) {
  # Create a .tex file to store the output
  tex_file <- paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/04_severity_modelling/", conf_level, "_table_", corstr_sev, ".tex")
  sink(tex_file)
  
  # Add the required LaTeX package
  #cat("\\usepackage{threeparttable}\n\n")
  
  # Loop through the list of data frames (tables)
  for (name in names(my_list)) {
    caption_text <- paste("Summary table for age", gsub("age", "", name), "- severity -", corstr_sev)
    
    #cat("\\section*{", name, "}\n")
    cat("\\begin{table}[ht]\n")
    cat("\\centering\n")
    cat("\\begin{threeparttable}\n")
    cat("\\caption{", caption_text, "}\n")
    
    # Custom header
    cat("\\scriptsize\n")
    cat("\\begin{tabular}{rlcccccl}\n")
    cat("\\hline\n")
    cat("Variable & Estimate & SE & \\makecell{Standardized\\\\Estimate} & 95\\% CI & \\makecell{James-Stein\\\\Estimate} & \\makecell{95\\% CI\\\\(James-Stein)} \\\\\n")
    
    # Print the table body without column names
    print(xtable(my_list[[name]], align = c("r", "l", "c", "c", "c", "c", "c", "c")), 
          include.rownames = FALSE, include.colnames = FALSE, type = "latex", only.contents = TRUE)
    
    cat("\\end{tabular}\n")
    
    # Add custom notes at the end of the table
    cat("\\begin{tablenotes}[para,flushleft]\n")
    cat("\\scriptsize\n")
    cat("\\item Superscripts $*+$ and $*-$ denote significant positive and negative effects at the 5\\% significance level, respectively.\n")
    cat("\\end{tablenotes}\n")
    
    cat("\\end{threeparttable}\n")
    cat("\\end{table}\n\n")
  }
  
  # Stop writing to the .tex file
  sink()
}


# Write the 95% CI and 90% CI tables to LaTeX
write_table_to_latex(summary_95, corstr_sev, 95)
write_table_to_latex(summary_90, corstr_sev, 90)

