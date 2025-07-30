library(openxlsx)

# Check if running in a SLURM environment
if (!is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))) {
  # If in SLURM environment
  setwd("inSLURM/path/to/Fluorosis/")
} else {
  # If not in SLURM environment
  setwd("SLURM/path/to/Fluorosis/")
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

get_SE = function(BS_dat, vars){
  data.frame(
    Variable = vars,
    SE = apply(BS_dat, MARGIN = 1, FUN = function(x) sd(x, na.rm = TRUE))
  )
}


# Helper function to write CI results to Excel
write_table_to_excel <- function(CI_list, corstr_pres, corstr_sev, conf_level) {
  wb <- createWorkbook()
  
  for (name in names(CI_list)) {
    addWorksheet(wb, name)
    writeData(wb, sheet = name, CI_list[[name]])
  }
  
  saveWorkbook(wb, paste0("path/to/Fluorosis/Results/05b_combined_severity_modelling/",
                          conf_level, "_table_", corstr_pres, ",", corstr_sev, ".xlsx"), overwrite = TRUE)
}

# Function to perform bootstrapping and save CI results
process_age_data <- function(age, B, corstr_pres, corstr_sev, vars) {
  stdcoefs_BS_mat <- matrix(NA, nrow = length(vars), ncol = B)
  path <- paste0("path/to/Fluorosis/Results/05b_combined_severity_modelling/",
                 corstr_pres, ",", corstr_sev, "/bootstrapping/age", age, "/")
  
  for (b in 1:B) {
    fn_stdcoef_sev <- paste0(path, "std_coefs_sev_age", age, "_b", b, ".RData")
    if (file.exists(fn_stdcoef_sev)) {
      obj_stdcoef_sev <- load(fn_stdcoef_sev)
      assign("std_coef_sev", get(obj_stdcoef_sev))
      stdcoefs_BS_mat[, b] <- std_coef_sev$`Standardized Estimate`[-c(1:4)]  # Adjusted to remove the first 4 rows
    }
  }
  
  return(stdcoefs_BS_mat)
}

# Main logic
corstr_pres <- "independence" #"exchangeable" #"jackknifed" #"ar1"
corstr_sev = "independence" #"exchangeable" #"jackknifed" #"ar1"

ages <- c(9, 13, 17, 23)
B <- 120
vars <- c("dental_age", "Total_mgF", "SugarAddedBeverageOzPerDay", "BrushingFrequencyPerDay",
          "Avg_homeppm", "Prop_DentAppt", "Prop_FluorideTreatment", "Tooth8", "Tooth9",
          "Tooth10", "ZoneM", "ZoneI", "ZoneE"   )

stdcoefs_BS_list <- lapply(ages, process_age_data, B = B, corstr_pres = corstr_pres, corstr_sev = corstr_sev, vars = vars)
names(stdcoefs_BS_list) <- paste0("age", ages)


# Calculate CIs
preJS_CI_95 <- lapply(stdcoefs_BS_list, get_CI, confidence_level = 0.95, vars = vars)
preJS_SE = lapply(stdcoefs_BS_list, get_SE, var=vars) #this is the SE for the stdcoefs_BS_list


# Perform James-Stein shrinkage
stdcoefsJS_BS_list <- stdcoefs_BS_list

for (var in vars) {
  var_ind <- which(vars %in% var)
  
  for (b in 1:B) {
    t <- c(
      (stdcoefs_BS_list[[1]])[var_ind, b],
      (stdcoefs_BS_list[[2]])[var_ind, b],
      (stdcoefs_BS_list[[3]])[var_ind, b],
      (stdcoefs_BS_list[[4]])[var_ind, b]
    )
    
    t = js(t)
    
    # Update the list with James-Stein shrinkage results
    stdcoefsJS_BS_list[[1]][var_ind, b] <- t[1]
    stdcoefsJS_BS_list[[2]][var_ind, b] <- t[2]
    stdcoefsJS_BS_list[[3]][var_ind, b] <- t[3]
    stdcoefsJS_BS_list[[4]][var_ind, b] <- t[4]
  }
}


# Calculate CIs for post-JS coefficients
postJS_CI_95 <- lapply(stdcoefsJS_BS_list, get_CI, confidence_level = 0.95, vars = vars)
postJS_SE = lapply(stdcoefsJS_BS_list, get_SE, var=vars) #this is the SE for the stdcoefs_JS_BS_list



#Obtain the whole data based estimates and SE:

coef_sev_list = list()
se_sev_list = list()

for(age in ages)
{
  path1 = paste0("path/to/Fluorosis/Results/05b_combined_severity_modelling/", corstr_pres,
                 ",", corstr_sev, "/whole_data_based/coefs_sev_age", age, ".RData")
  obj_coef_sev <- load(path1)
  assign("coef_sev", get(obj_coef_sev))
  
  
  path2 = paste0("path/to/Fluorosis/Results/05b_combined_severity_modelling/", corstr_pres, 
                 ",", corstr_sev, "/whole_data_based/JK_SD_sev_age", age, ".RData")
  obj_se_sev <- load(path2)
  assign("se_sev", get(obj_se_sev))
  
  
  age_ind = which(ages %in% age)
  
  coef_sev_list[[age_ind]] = coef_sev[-c(1:4),]
  se_sev_list[[age_ind]] = se_sev[-c(1:4), , drop=F]
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
  
 
  preJS_CI = preJS_CI_95
  postJS_CI = postJS_CI_95
  
  
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
    
    # Determine significance and direction symbols
    ind_preJS_signif <- preJS_CI[[age_ind]]$LCL * preJS_CI[[age_ind]]$UCL > 0
    preJS_signif <- ifelse(ind_preJS_signif, "*", "")
    preJS_direction <- ifelse(coef_sev_list[[age_ind]]$Estimates > 0, "+", "-")
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
    postJS_direction <- ifelse(coef_sev_list[[age_ind]]$Estimates > 0, "+", "-")
    postJS_symbol <- ifelse(ind_postJS_signif, paste0("$^{", postJS_signif, postJS_direction, "}$"), "")
    postJS_CI_int <- paste0(postJS_CI_int, postJS_symbol)
    
    temp_df  = data.frame(Variable=rownames(se_sev_list$age9),
                          Estimate=format(round((coef_sev_list[[age_ind]])$Estimates, 3), nsmall=3),
                          Std_Estimate= format(round((coef_sev_list[[age_ind]])$Estimates/se_sev_list[[age_ind]]$SD, 3), nsmall=3),
                          Std_SE = format(round(preJS_SE[[age_ind]]$SE, 3), nsmall=3),
                          JS_Estimate = format(round(postJS_stdcoef_sev_list[[age_ind]]$Estimates, 3), nsmall=3),
                          JS_Bias = 
                            format(round( apply(stdcoefsJS_BS_list[[age_ind]], MARGIN = 1, FUN = function(x) mean(x, na.rm = TRUE)) -
                                            postJS_stdcoef_sev_list[[age_ind]]$Estimates , 3), nsmall=3),
                          JS_SE = format(round(postJS_SE[[age_ind]]$SE, 3), nsmall=3),
                          JS_MSE = format(round( (apply(stdcoefsJS_BS_list[[age_ind]], MARGIN = 1, FUN = function(x) mean(x, na.rm = TRUE)) -
                                                    postJS_stdcoef_sev_list[[age_ind]]$Estimates)^2 + postJS_SE[[age_ind]]$SE, 3), nsmall=3),
                          JS_CI=postJS_CI_int)
    
    colnames(temp_df) = c("Variable", "Estimate", "Standardized Estimate", "SE", "James-Stein Estimator",
                          "James-Stein Bias", "James-Stein SE", "James-Stein MSE", "James-Stein CI")
    
    summary_sev[[age_ind]] = temp_df
  }
  
  names(summary_sev) = paste0("age", ages)
  
  return(summary_sev)
}

summary_95 = create_summary_table(95)

write_table_to_excel(summary_95, corstr_pres, corstr_sev, 95)

#########################################################
######### Write latex code for Overleaf tables: #########
#########################################################
                                   
library(xtable)

write_table_to_latex <- function(my_list, corstr_pres, corstr_sev, conf_level) {
  # Create a .tex file to store the output
  tex_file <- paste0("path/to/Fluorosis/Results/05b_combined_severity_modelling/", 
                     conf_level, "_table_", corstr_pres, ",", corstr_sev, ".tex")
  sink(tex_file)
  
  corstr_pres_index = which(c("independence", "exchangeable", "ar1", "jackknifed") %in% corstr_pres)
  corstr_sev_index = which(c("independence", "exchangeable", "ar1", "jackknifed") %in% corstr_sev)
  
  #Model name: A.1.1 means A independence age
  caption_text <- paste0("Severity estimates from models C.", corstr_pres_index, ".", corstr_sev_index,
                         ".1-C.", corstr_pres_index, ".", corstr_sev_index,
                         ".4, the combined models with the ",
                         corstr_pres, " and ", corstr_sev, " presence and severity cluster correlation structures respectively.")
  
  # gsub("age", "", name)
  
  #cat("\\section*{", name, "}\n")
  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat("\\caption{", caption_text, "}\n")
  cat("\\label{table:comb:sev:", corstr_sev, "}\n", sep="")
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
    cat("\\caption{Model C.", corstr_pres_index, ".", corstr_sev_index,
        ".", age_index, " (age ",age,")}\n", sep="")
    cat("\\scalebox{0.55}{\n")
    cat("\\begin{tabular}{rrrrrrrrl}\n")
    #cat("\\hline\n")
    
    cat("Variable & Estimate & \\makecell{Standardized\\\\Estimate} & \\makecell{SE\\\\(Standardized\\\\Estimate)} & \\makecell{James-Stein\\\\Estimator} & \\makecell{Bias\\\\(James-Stein\\\\Estimator)} & \\makecell{SE\\\\(James-Stein\\\\Estimator)} & \\makecell{MSE\\\\(James-Stein\\\\Estimator)} & \\makecell{95\\% CI\\\\(James-Stein\\\\Estimator)} \\\\\n")
    
    # Print the table body without column names
    print(xtable(my_list[[name]], align = c("r", "r", "r", "r", "r", "r", "r", "r", "r", "l")), 
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


# Write the 95% CI table to LaTeX
write_table_to_latex(summary_95, corstr_pres, corstr_sev, 95)
