##############################################################################################################
######################################      FUNCTIONS   ######################################################
##############################################################################################################

# List of required packages
required_packages <- c("dplyr", "data.table", "reshape2", "parallel", "MASS", "geepack", "Matrix")

# Function to check and install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)  # Load the package after installation
    }
  }
}

# Load all packages at once
lapply(required_packages, library, character.only = TRUE)



get_corr_submat <- function(master_mat, selected_combinations){ 
  m.x <- melt(master_mat)
  ind.subi <- expand.grid(selected_combinations, selected_combinations)
  n.x <- ind.subi %>% inner_join(m.x) %>% dcast(Var1~Var2, value.var = 'value')
  
  mat = as.matrix(n.x[,-c(1)])
  rownames(mat) = colnames(mat)
  
  return(mat)
}

# Function to create the matrix of multinomial correlations for a given row
get_multinom_corrmat <- function(pi1, pi2, pi3) {
  M <- matrix(0, nrow = 3, ncol = 3)
  
  # Assign the upper triangular part
  M[1,1] <- 1
  M[1,2] <- -pi1 * (1 - pi2) / sqrt( pi1 * (1 - pi1)  *  pi2 * (1 - pi2))
  M[1,3] <- -pi1 * (1 - pi3)  /  sqrt( pi1 * (1 - pi1)  *  pi3 * (1 - pi3))
  
  M[2,2] <- 1
  M[2,3] <- -pi2 * (1 - pi3) / sqrt( pi2 * (1 - pi2)  *  pi3 * (1 - pi3))
  
  M[3,3] <- 1
  
  # Fill in the lower triangular part using symmetry
  M[2,1] <- M[1,2]
  M[3,1] <- M[1,3]
  M[3,2] <- M[2,3]
  
  return(M)
}

get_tooth_index = function(dat = IFS_dat)
{
  Tooth = ifelse(dat$Tooth8==1, 8, ifelse(dat$Tooth9==1, 9, ifelse(dat$Tooth10==1, 10, 7)))
  return(Tooth)
}

get_zone_index = function(dat = IFS_dat)
{
  Zone = ifelse(dat$ZoneM==1, "M", ifelse(dat$ZoneI==1, "I", ifelse(dat$ZoneO==1, "O", "C")))
  return(Zone)
}

#for scalar input
inv_logit <- function(x) {
  ifelse(is.infinite(exp(x)), 1, exp(x) / (1 + exp(x)))
}

inv_logit_vec <- function(x) {
  exp_x <- exp(x)  # Precompute exp(x)
  
  # Handle large x values where exp(x) becomes Inf
  result <- exp_x / (1 + exp_x)
  
  # For very large x, exp(x) becomes Inf, so we set those directly to 1
  result[is.infinite(exp_x)] <- 1
  
  return(result)
}

# Function to repeat columns consecutively
repeat_column <- function(column, times=3) {
  matrix(rep(column, times), ncol = times, byrow = FALSE)
}

#function to create an exchangeable matrix
create_exchangeable_matrix <- function(n, rho) {
  # Create an n x n matrix filled with rho
  mat <- matrix(rho, n, n)
  # Replace the diagonal with 1s
  diag(mat) <- 1
  return(mat)
}

transform_to_vec1 = function(x)
{
  return(c(x, -x, 0))
}

transform_to_vec2 = function(x)
{
  return(c(0, x, -x))
}

#Function to calculate presence model's exchangeable correlation parameter:
get_rho_pres_exch = function(dat_pres = dat_pres,
                             init_pres_coef = pres_coef$Coefficient)  #note: the first element should be alpha and the next, betas 
  #for combined, init_pres_coef is init_pres_coef is loaded from the separate models. 
  #For separate presence model, init_pres_coef is the pres_coef$Coefficient calculated from geeglm.
{
  #Presence:
  
  rho_pres = 1
  epsilon <- .Machine$double.eps  # Smallest positive number such that 1 + epsilon != 1
  
  dat2 = dat_pres
  
  dat2$Wp = dat_pres$FRI_Score
  
  dat2$mu = apply(dat_pres[,-c(1,2)], MARGIN = 1, FUN = function(x){ 
    mu_val <- inv_logit(init_pres_coef[1] + sum(init_pres_coef[-c(1)]*x))
    # Adjust mu_val to avoid it being too close to 0 or 1
    pmax(pmin(mu_val, 1 - epsilon), epsilon)
  })
  
  dat2$var = dat2$mu*(1-dat2$mu)
  dat2$r = (dat2$Wp - dat2$mu)/sqrt(dat2$var)
  
  phi = 1 / ( sum(dat2$r^2/ (nrow(dat2) - ncol(dat_pres[,-c(1,2)]))) )
  
  prods = outer(dat2$r, dat2$r, FUN = function(x, y){x*y})
  
  compare_inds_i = outer(1:nrow(dat_pres), 1:nrow(dat_pres), FUN = function(x, y){x==y})
  compare_inds_ZoneM = outer(dat_pres$ZoneM, dat_pres$ZoneM, FUN = function(x, y){x==y})
  compare_inds_ZoneI = outer(dat_pres$ZoneI, dat_pres$ZoneI, FUN = function(x, y){x==y})
  compare_inds_ZoneO = outer(dat_pres$ZoneO, dat_pres$ZoneO, FUN = function(x, y){x==y})
  compare_inds_Tooth8 = outer(dat_pres$Tooth8, dat_pres$Tooth8, FUN = function(x, y){x==y})
  compare_inds_Tooth9 = outer(dat_pres$Tooth9, dat_pres$Tooth9, FUN = function(x, y){x==y})
  compare_inds_Tooth10 = outer(dat_pres$Tooth10, dat_pres$Tooth10, FUN = function(x, y){x==y})
  
  keep_ind = !(compare_inds_i | 
                 ((compare_inds_ZoneM & compare_inds_ZoneI & compare_inds_ZoneO) &
                    (compare_inds_Tooth8 & compare_inds_Tooth9 & compare_inds_Tooth10)) 
  )
  
  keep_prods = prods[keep_ind]
  
  rho_pres = phi * mean(keep_prods)
  
  #Correct for overflow/underflow
  rho_pres = ifelse(rho_pres< -1, -1, rho_pres)
  rho_pres = ifelse(rho_pres> 1, 1, rho_pres)
  
  return(rho_pres)
}


#Function to calculate severity model's exchangeable correlation parameter:
get_rho_sev_exch = function(dat_sev=dat_sev, 
                            init_sev_coef=init_sev_coef)
{
  rho_sev = 1
  epsilon <- .Machine$double.eps  # Smallest positive number such that 1 + epsilon != 1
  
  dat2 = dat_sev
  
  dat2$Ws = dat2$FRI_Score
  #dat2$zeta = ifelse(dat2$Ws == 1, sev_coef[1], sev_coef[2])
  
  dat2$cusum1 = apply(dat_sev[,-c(1,2)], MARGIN = 1, FUN = function(x){ 
    cusum1_val <- inv_logit(init_sev_coef[1] + sum(init_sev_coef[-c(1,2)] * x))
    #use pmax and pmin to ensure cusum1_val is within a safe range
    pmax(pmin(cusum1_val, 1 - epsilon), epsilon)
  })
  
  dat2$cusum2 = apply(dat_sev[,-c(1,2)], MARGIN = 1, FUN = function(x){ 
    cusum2_val <- inv_logit(init_sev_coef[2] + sum(init_sev_coef[-c(1,2)] * x))
    #use pmax and pmin to ensure cusum2_val is within a safe range
    pmax(pmin(cusum2_val, 1 - epsilon), epsilon)
  })
  
  
  dat2$cusum3 = rep(1, length(dat2$cusum2))
  
  dat2$pi1 = dat2$cusum1
  dat2$pi2 = dat2$cusum2 - dat2$cusum1
  dat2$pi3 = 1 - dat2$cusum2
  
  dat2$Zminuspi1 = ifelse(dat2$Ws==1, 1, 0) - dat2$pi1
  dat2$Zminuspi2 = ifelse(dat2$Ws==2, 1, 0) - dat2$pi2
  dat2$Zminuspi3 = ifelse(dat2$Ws==3, 1, 0) - dat2$pi3
  
  dat2$varZ1 = dat2$pi1*(dat2$pi1)
  dat2$varZ2 = dat2$pi2*(dat2$pi2)
  dat2$varZ3 = dat2$pi3*(dat2$pi3)
  
  
  #If variance is too close to zero, adjust 
  dat2$varZ1 = ifelse(dat2$varZ1==0, .Machine$double.eps, dat2$varZ1)
  dat2$varZ2 = ifelse(dat2$varZ2==0, .Machine$double.eps, dat2$varZ2)
  dat2$varZ3 = ifelse(dat2$varZ3==0, .Machine$double.eps, dat2$varZ3)
  
  
  dat2$std_Zminuspi1 = dat2$Zminuspi1/sqrt(dat2$varZ1)
  dat2$std_Zminuspi2 = dat2$Zminuspi2/sqrt(dat2$varZ2)
  dat2$std_Zminuspi3 = dat2$Zminuspi3/sqrt(dat2$varZ3)
  
  resids1 = dat2[c("SUBJECT_ID", "Ws", "Tooth8", "Tooth9", "Tooth10", "ZoneM", "ZoneI", "ZoneO", "std_Zminuspi1")]
  colnames(resids1)[ncol(resids1)] = "std_Zminuspi"
  resids2 = dat2[c("SUBJECT_ID", "Ws", "Tooth8", "Tooth9", "Tooth10", "ZoneM", "ZoneI", "ZoneO", "std_Zminuspi2")]; colnames(resids2)[ncol(resids2)] = "std_Zminuspi"
  resids3 = dat2[c("SUBJECT_ID", "Ws" ,"Tooth8", "Tooth9", "Tooth10", "ZoneM", "ZoneI", "ZoneO", "std_Zminuspi3")]; colnames(resids3)[ncol(resids3)] = "std_Zminuspi"
  
  stacked_resids = rbind(resids1, resids2, resids3)
  
  phi = 1 / ( sum(stacked_resids$std_Zminuspi^2/ (nrow(stacked_resids) - ncol(dat_sev[,-c(1,2)]))) )
  
  prods = outer(stacked_resids$std_Zminuspi, stacked_resids$std_Zminuspi, FUN = function(x,y){x*y})
  
  compare_ind_Tooth8 = outer(stacked_resids$Tooth8, stacked_resids$Tooth8, FUN = function(x,y){x==y})
  compare_ind_Tooth9 = outer(stacked_resids$Tooth9, stacked_resids$Tooth9, FUN = function(x,y){x==y})
  compare_ind_Tooth10 = outer(stacked_resids$Tooth10, stacked_resids$Tooth10, FUN = function(x,y){x==y})
  compare_categ = outer(stacked_resids$Ws, stacked_resids$Ws, FUN = function(x,y){x==y})
  
  keep_ind = !(compare_ind_Tooth8 & compare_ind_Tooth9 & compare_ind_Tooth10 & compare_categ)
  
  keep_prods = prods[keep_ind]
  
  rho_sev = phi * mean(keep_prods)
  
  #Correct for overflow/underflow
  rho_sev = ifelse(rho_sev< -1, -1, rho_sev)
  rho_sev = ifelse(rho_sev> 1, 1, rho_sev)
  
  return(rho_sev)
}

#Function to calculate presence model's ar1 correlation parameter:
get_rho_pres_ar1 = function(dat_pres=dat_pres, 
                            init_pres_coef=init_pres_coef)
{
  rho_pres = 1
  epsilon <- .Machine$double.eps  # Smallest positive number such that 1 + epsilon != 1
  
  dat2 = dat_pres
  
  dat2$Wp = dat_pres$FRI_Score
  
  dat2$mu = apply(dat_pres[,-c(1,2)], MARGIN = 1, FUN = function(x){ 
    mu_val <- inv_logit(init_pres_coef[1] + sum(init_pres_coef[-c(1)]*x))
    # Adjust mu_val to avoid it being too close to 0 or 1
    pmax(pmin(mu_val, 1 - epsilon), epsilon)
  })
  
  dat2$var = dat2$mu*(1-dat2$mu)
  dat2$r = (dat2$Wp - dat2$mu)/sqrt(dat2$var)
  dat2$Tooth = get_tooth_index(dat_pres) #apply(dat_pres, 1, get_tooth_index)
  
  all_resids = as.matrix(dat2$r, ncol = 1)
  rrprod_mat = (all_resids %*% t(all_resids))
  keep_ind = rrprod_mat>0
  log_rrprod_mat = log(rrprod_mat[keep_ind])
  
  abstoothdiff_mat = outer(dat2$Tooth, dat2$Tooth, FUN = function(x, y) abs(x - y))
  sub_abstoothdiff_mat = abstoothdiff_mat[keep_ind]
  
  ind_abstoothdiff_zero = which(abstoothdiff_mat != 0)
  
  sub_log_mat = log_rrprod_mat[ind_abstoothdiff_zero]
  sub_abstoothdiff_mat = sub_abstoothdiff_mat[ind_abstoothdiff_zero]
  
  all_log = as.vector(sub_log_mat)
  all_abstoothdiff = as.vector(sub_abstoothdiff_mat)
  
  rho_pres = lm(all_log ~ all_abstoothdiff)$coefficients[2]
  
  rho_pres = ifelse(rho_pres< -1, -1, rho_pres)
  rho_pres = ifelse(rho_pres> 1, 1, rho_pres)
  
  return(rho_pres)
}

#Function to calculate severity model's exchangeable correlation parameter:
get_rho_sev_ar1 = function(dat_sev=dat_sev, 
                           init_sev_coef=init_sev_coef)
{
  epsilon <- 1e-6
  
  rho_sev = 1
  
  dat2 = dat_sev
  
  dat2$Ws = dat2$FRI_Score
  #dat2$zeta = ifelse(dat2$Ws == 1, sev_coef[1], sev_coef[2])
  
  dat2$cusum1 = apply(dat_sev[,-c(1,2)], MARGIN = 1, FUN = function(x){ 
    cusum1_val <- inv_logit(init_sev_coef[1] + sum(init_sev_coef[-c(1,2)] * x))
    #use pmax and pmin to ensure cusum1_val is within a safe range
    pmax(pmin(cusum1_val, 1 - epsilon), epsilon)
  })
  
  dat2$cusum2 = apply(dat_sev[,-c(1,2)], MARGIN = 1, FUN = function(x){ 
    cusum2_val <- inv_logit(init_sev_coef[2] + sum(init_sev_coef[-c(1,2)] * x))
    #use pmax and pmin to ensure cusum2_val is within a safe range
    pmax(pmin(cusum2_val, 1 - epsilon), epsilon)
  })
  
  
  dat2$cusum3 = rep(1, length(dat2$cusum2))
  
  dat2$pi1 = dat2$cusum1
  dat2$pi2 = dat2$cusum2 - dat2$cusum1
  dat2$pi3 = 1 - dat2$cusum2
  
  dat2$Zminuspi1 = ifelse(dat2$Ws==1, 1, 0) - dat2$pi1
  dat2$Zminuspi2 = ifelse(dat2$Ws==2, 1, 0) - dat2$pi2
  dat2$Zminuspi3 = ifelse(dat2$Ws==3, 1, 0) - dat2$pi3
  
  dat2$varZ1 = dat2$pi1*(dat2$pi1)
  dat2$varZ2 = dat2$pi2*(dat2$pi2)
  dat2$varZ3 = dat2$pi3*(dat2$pi3)
  
  dat2$std_Zminuspi1 = dat2$Zminuspi1/sqrt(dat2$varZ1)
  dat2$std_Zminuspi2 = dat2$Zminuspi2/sqrt(dat2$varZ2)
  dat2$std_Zminuspi3 = dat2$Zminuspi3/sqrt(dat2$varZ3)
  
  resids1 = dat2[c("SUBJECT_ID", "Ws", "Tooth8", "Tooth9", "Tooth10", "ZoneM", "ZoneI", "ZoneO", "std_Zminuspi1")]
  colnames(resids1)[ncol(resids1)] = "std_Zminuspi"
  resids2 = dat2[c("SUBJECT_ID", "Ws", "Tooth8", "Tooth9", "Tooth10", "ZoneM", "ZoneI", "ZoneO", "std_Zminuspi2")]; colnames(resids2)[ncol(resids2)] = "std_Zminuspi"
  resids3 = dat2[c("SUBJECT_ID", "Ws" ,"Tooth8", "Tooth9", "Tooth10", "ZoneM", "ZoneI", "ZoneO", "std_Zminuspi3")]; colnames(resids3)[ncol(resids3)] = "std_Zminuspi"
  
  stacked_resids = rbind(resids1, resids2, resids3)
  
  prods = outer(stacked_resids$std_Zminuspi, stacked_resids$std_Zminuspi, FUN = function(x,y){x*y})
  log_prods = suppressWarnings(log(prods))
  
  compare_ind_Tooth8 = outer(stacked_resids$Tooth8, stacked_resids$Tooth8, FUN = function(x,y){x==y})
  compare_ind_Tooth9 = outer(stacked_resids$Tooth9, stacked_resids$Tooth9, FUN = function(x,y){x==y})
  compare_ind_Tooth10 = outer(stacked_resids$Tooth10, stacked_resids$Tooth10, FUN = function(x,y){x==y})
  
  keep_ind = !(compare_ind_Tooth8 & compare_ind_Tooth9 & compare_ind_Tooth10)
  
  keep_log_prods = log_prods[keep_ind]
  
  stacked_resids$Tooth = ifelse(stacked_resids$Tooth8==1, 8, ifelse(stacked_resids$Tooth9==1, 9, 10))
  abstoothdiff_mat = outer(stacked_resids$Tooth, stacked_resids$Tooth, FUN = function(x,y){abs(x-y)})
  keep_abstoothdiff_mat = abstoothdiff_mat[keep_ind]
  
  temp_df = na.omit(data.frame(logprods = as.vector(keep_log_prods), absdiff = as.vector(keep_abstoothdiff_mat)))
  temp_df = temp_df[temp_df$logprods!=Inf,] ; temp_df = temp_df[temp_df$logprods!=-Inf,]
  
  rho_sev = lm(logprods ~ absdiff, data = temp_df)$coefficients[2]
  
  rho_sev = ifelse(rho_sev > 1, 1, rho_sev)
  rho_sev = ifelse(rho_sev< -1, -1, rho_sev)
  
  return(rho_sev)
}

#Function to get presence model's jackknifed correlation matrix
get_JK_corrmat_pres <- function(dat_pres, init_pres_coef) {
  epsilon <- 1e-6
  
  N <- length(unique(dat_pres$SUBJECT_ID))
  
  # Pre-allocate the matrix
  JK_corr_mat_pres <- matrix(0, nrow = 16, ncol = 16)
  
  tooth_zone_comb <- expand.grid(Tooth = 7:10, Zone = c("C", "M", "I", "O"))
  tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
  rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
  
  # Set up parallel backend
  num_cores <- detectCores() - 1  # Use one less than total cores
  cl <- makeCluster(num_cores)
  
  # Export required variables and functions to cluster workers
  clusterExport(cl, c("dat_pres", "tooth_zone_comb", "inv_logit", "init_pres_coef", "get_tooth_index", "get_zone_index", "epsilon"), envir = environment())
  
  # Parallelize the jackknife step using parLapply
  results <- parLapply(cl, unique(dat_pres$SUBJECT_ID), function(id, dat_pres_local, init_pres_coef, tooth_zone_comb, epsilon) {
    # Subset the data for jackknife (leaving one subject out)
    dat_pres_JK <- dat_pres_local[dat_pres_local$SUBJECT_ID != id, ]
    
    # Recreate dat2 locally on each worker
    dat2 <- dat_pres_JK
    dat2$Wp <- dat_pres_JK$FRI_Score
    dat2$mu <- pmax(pmin(inv_logit(init_pres_coef[1] + rowSums(t(t(dat_pres_JK[, -c(1, 2)]) * init_pres_coef[-c(1)]))), 1 - epsilon), epsilon)
    dat2$var <- dat2$mu * (1 - dat2$mu)
    dat2$r <- (dat2$Wp - dat2$mu) / sqrt(dat2$var)
    dat2$Tooth <- get_tooth_index(dat_pres_JK)
    dat2$Zone <- get_zone_index(dat_pres_JK)
    
    # Calculate phi
    phi <- 1 / (sum(dat2$r^2) / (nrow(dat2) - ncol(dat_pres_JK[, -c(1, 2)])))
    
    # Precompute outer product for all combinations in one go
    cov_mat <- matrix(0, nrow = 16, ncol = 16)
    
    for (i in 1:nrow(tooth_zone_comb)) {
      for (j in i:nrow(tooth_zone_comb)) {
        tooth1 <- tooth_zone_comb$Tooth[i]
        zone1 <- tooth_zone_comb$Zone[i]
        tooth2 <- tooth_zone_comb$Tooth[j]
        zone2 <- tooth_zone_comb$Zone[j]
        
        if (i != j) {
          ind1 <- dat2$Tooth == tooth1 & dat2$Zone == zone1
          ind2 <- dat2$Tooth == tooth2 & dat2$Zone == zone2
          
          subdat_pres_JK1 <- dat2[ind1, ]
          subdat_pres_JK2 <- dat2[ind2, ]
          
          # Calculate covariance for (tooth1, zone1) and (tooth2, zone2)
          temp <- mean(outer(subdat_pres_JK1$r, subdat_pres_JK2$r, "*"))
          
          cov_mat[i, j] <- cov_mat[j, i] <- phi * temp
        }
      }
    }
    return(cov_mat)
  }, dat_pres_local = dat_pres, init_pres_coef = init_pres_coef, tooth_zone_comb = tooth_zone_comb, epsilon = epsilon)
  
  stopCluster(cl)  # Stop the cluster when done
  
  # Aggregate results
  JK_corr_mat_pres <- Reduce(`+`, results) / N
  
  JK_corr_mat_pres[JK_corr_mat_pres > 1] <- 1
  JK_corr_mat_pres[JK_corr_mat_pres < -1] <- -1
  
  JK_corr_mat_pres <- (JK_corr_mat_pres + t(JK_corr_mat_pres)) / 2  # Fix rounding issues
  diag(JK_corr_mat_pres) <- 1
  
  JK_corr_mat_pres[is.nan(JK_corr_mat_pres)] <- 0
  
  # Print progress
  print("Jackknife completed")
  
  return(JK_corr_mat_pres)
}


get_JK_corrmat_sev <- function(dat_sev = dat_sev, 
                               init_sev_coef = init_sev_coef) {
  
  # Severity:
  dat2 <- dat_sev
  dat2$Ws <- dat2$FRI_Score
  epsilon <- 1e-6
  
  # Vectorize the cusum calculations
  x_matrix <- as.matrix(dat_sev[, -c(1, 2)])  # Extract relevant data as a matrix
  dat2$cusum1 <- pmax(pmin(inv_logit(init_sev_coef[1] + x_matrix %*% init_sev_coef[-c(1, 2)]), 1 - epsilon), epsilon)
  dat2$cusum2 <- pmax(pmin(inv_logit(init_sev_coef[2] + x_matrix %*% init_sev_coef[-c(1, 2)]), 1 - epsilon), epsilon)
  dat2$cusum3 <- rep(1, nrow(dat2))  # Cusum3 is always 1
  
  # Pi calculations
  dat2$pi1 <- dat2$cusum1
  dat2$pi2 <- dat2$cusum2 - dat2$cusum1
  dat2$pi3 <- 1 - dat2$cusum2
  
  # Zminuspi calculations
  dat2$Zminuspi1 <- ifelse(dat2$Ws == 1, 1, 0) - dat2$pi1
  dat2$Zminuspi2 <- ifelse(dat2$Ws == 2, 1, 0) - dat2$pi2
  dat2$Zminuspi3 <- ifelse(dat2$Ws == 3, 1, 0) - dat2$pi3
  
  # Variance calculations
  dat2$varZ1 <- dat2$pi1 * dat2$pi1
  dat2$varZ2 <- dat2$pi2 * dat2$pi2
  dat2$varZ3 <- dat2$pi3 * dat2$pi3
  
  #If variance is too close to zero, adjust 
  dat2$varZ1 = ifelse(dat2$varZ1==0, .Machine$double.eps, dat2$varZ1)
  dat2$varZ2 = ifelse(dat2$varZ2==0, .Machine$double.eps, dat2$varZ2)
  dat2$varZ3 = ifelse(dat2$varZ3==0, .Machine$double.eps, dat2$varZ3)
  
  # Standardized residuals
  dat2$std_Zminuspi1 <- dat2$Zminuspi1 / sqrt(dat2$varZ1)
  dat2$std_Zminuspi2 <- dat2$Zminuspi2 / sqrt(dat2$varZ2)
  dat2$std_Zminuspi3 <- dat2$Zminuspi3 / sqrt(dat2$varZ3)
  
  dat2$Tooth <- get_tooth_index(dat_sev)
  dat2$Zone <- get_zone_index(dat_sev)
  
  # Combine residuals
  resids_list <- lapply(1:3, function(i) {
    resids <- dat2[, c("SUBJECT_ID", "Tooth", "Zone", paste0("std_Zminuspi", i))]
    colnames(resids)[ncol(resids)] <- "std_Zminuspi"
    resids$SUBJECT_ID <- resids$SUBJECT_ID + (i - 1) * (max(dat2$SUBJECT_ID) + 1)
    return(resids)
  })
  stacked_resids <- do.call(rbind, resids_list)
  
  N <- length(unique(stacked_resids$SUBJECT_ID))
  
  JK_corr_mat_sev <- matrix(0, nrow = 16, ncol = 16)  # Pre-allocate the matrix
  
  tooth_zone_comb <- expand.grid(Tooth = 7:10, Zone = c("C", "M", "I", "O"))
  tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
  rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
  
  # Set up parallel backend
  if (!requireNamespace("parallel", quietly = TRUE)) {
    install.packages("parallel")
  }
  num_cores <- detectCores() - 1  # Use one less than total cores
  cl <- makeCluster(num_cores)
  
  # Parallelize the jackknife step and pass stacked_resids as an argument
  results <- parLapply(cl, unique(stacked_resids$SUBJECT_ID), function(id, stacked_resids_local, tooth_zone_comb, dat_sev_local, epsilon) {
    
    stacked_resids_JK <- stacked_resids_local[stacked_resids_local$SUBJECT_ID != id, ]
    
    phi <- 1 / (sum(stacked_resids_JK$std_Zminuspi^2) / (nrow(stacked_resids_JK) - ncol(dat_sev_local[,-c(1,2)])))
    
    cov_mat <- matrix(0, nrow = 16, ncol = 16)
    
    for (i in 1:nrow(tooth_zone_comb)) {
      for (j in i:nrow(tooth_zone_comb)) {
        tooth1 <- tooth_zone_comb$Tooth[i]
        zone1 <- tooth_zone_comb$Zone[i]
        tooth2 <- tooth_zone_comb$Tooth[j]
        zone2 <- tooth_zone_comb$Zone[j]
        
        if (i != j) {
          ind1 <- stacked_resids_JK$Tooth == tooth1 & stacked_resids_JK$Zone == zone1
          ind2 <- stacked_resids_JK$Tooth == tooth2 & stacked_resids_JK$Zone == zone2
          
          substacked_resids_JK_1 <- stacked_resids_JK[ind1, ]
          substacked_resids_JK_2 <- stacked_resids_JK[ind2, ]
          
          temp <- mean(outer(substacked_resids_JK_1$std_Zminuspi, substacked_resids_JK_2$std_Zminuspi, "*"))
          
          cov_mat[i, j] <- cov_mat[j, i] <- phi * temp
        }
      }
    }
    return(cov_mat)
  }, stacked_resids_local = stacked_resids, tooth_zone_comb = tooth_zone_comb, dat_sev_local = dat_sev, epsilon = epsilon)
  
  stopCluster(cl)
  
  # Aggregate results without a loop
  JK_corr_mat_sev <- Reduce(`+`, results) / N
  
  JK_corr_mat_sev[JK_corr_mat_sev > 1] <- 1
  JK_corr_mat_sev[JK_corr_mat_sev < -1] <- -1
  diag(JK_corr_mat_sev) <- 1
  JK_corr_mat_sev <- (JK_corr_mat_sev + t(JK_corr_mat_sev)) / 2  # Fix rounding issues
  
  JK_corr_mat_sev[is.nan(JK_corr_mat_sev)] <- 0
  
  # Print completion message
  print("(Severity) Jackknife completed")
  
  return(JK_corr_mat_sev)
}



#Drop columns in a dataframe which have nonzero number of NAs
extract_non_NA_cols = function(df)
{
  return(df[ , colSums(is.na(df))==0])
}

#Compute jackknifed standard error:
compute_JK_SE = function(x){
  x=na.omit(x)
  T=length(x)
  SE = sqrt(((T-1)/T) * sum( (x-mean(x))^2 ))
  return(SE)
}


#Define function for James-Stein shrinkage:
js = function(x)
{
  G = length(x)
  result = rep(NA, times=G)
  
  if(!any(is.na(x)))
  {
    xbar = mean(x)
    result =  xbar + max(0, (1 - (G-2)/(sum(x^2))) ) * (x - xbar) 
  }
  
  return(result)
}


fit_GEE_presence <- function(dat = IFS_dat, kappa = 0.25, corstr_pres = "independence", 
                             maxIter = 10, tol = 1e-7, n_cores = min(16, detectCores() - 1),
                             reuse = F, rho_pres = NA, JK_corr_mat_pres = NA,
                             init_pres_coef = NA) {
  
  dat_pres <- dat
  dat_pres$FRI_Score <- ifelse(dat_pres$FRI_Score == 0, 1, 0)  # Convert FRI_Score to W_p (0/1)
  
  pres_coef = init_pres_coef
  
  if(is.na(init_pres_coef))
  {
    pres_mod <- geeglm(FRI_Score ~ . -SUBJECT_ID,
                       data = dat_pres, id = SUBJECT_ID, family = binomial(link = "logit"))
    
    pres_coef <- as.vector(pres_mod$coefficients)
  }
  
  
  n_beta <- length(pres_coef) - 1  # Exclude the intercept
  n_params <- n_beta + 1  # alpha + n_beta
  
  # Initialize the parameters matrix
  params <- matrix(NA, nrow = maxIter, ncol = n_params)
  colnames(params) <- c("alpha", names(pres_mod$coefficients[-c(1)]))
  params[1, ] <- pres_coef
  params[, 1] <- pres_coef[1]
  
  # Set correlation structure variables
  
  if(reuse==T) #if rho_pres or JK_corr_mat_pres is provided
  {
    if( !(!is.na(rho_pres) || !any(is.na(JK_corr_mat_pres))) )
    {
      print("reuse=T but correlation parameters were not passed.")
    }
    
  }else{ #if correlation parameters are not provided, estimate them.
    
    rho_pres <- 1
    
    if(corstr_pres %in% c("exchangeable", "ar1")) #need rho_pres
    {
      if (corstr_pres == "exchangeable") {
        rho_pres <- get_rho_pres_exch(dat_pres = dat_pres, init_pres_coef = pres_coef)
      } else if (corstr_pres == "ar1") {
        rho_pres <- get_rho_pres_ar1(dat_pres = dat_pres, init_pres_coef = pres_coef)
      } 
      
      rho_pres <- pmax(pmin(rho_pres, 1), -1)
      
    }else if(corstr_pres=="jackknifed") #need JK_corr_mat_pres
    {
      JK_corr_mat_pres <- data.frame(diag(4 * 4))
      tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
      tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
      rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
      
      JK_corr_mat_pres <- get_JK_corrmat_pres(dat_pres = dat_pres, init_pres_coef = pres_coef)
      rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
      
    }
    
    
  }
  
  
  # Set up parallel cluster
  cl <- makeCluster(n_cores)
  
  # Load the necessary libraries on each worker using a single call
  clusterEvalQ(cl, {
    library(MASS)
    library(reshape2)
    library(data.table)
    library(dplyr)
  })
  
  # Export necessary objects and functions to the parallel environment from the function's environment
  clusterExport(cl, varlist = c("dat_pres", "get_tooth_index", "get_zone_index", 
                                "inv_logit_vec", "rho_pres", "JK_corr_mat_pres", 
                                "corstr_pres", "get_corr_submat"), envir = environment())
  
  # Start the GEE iterations
  for (iter in 1:(maxIter - 1)) {
    #psi_pres <- numeric(n_params)
    #hessian_pres <- matrix(0, nrow = n_params, ncol = n_params)
    
    # Capture the current row of params to pass to the workers
    params_iter <- params[iter, ]
    
    # Re-export params_iter for each iteration
    clusterExport(cl, "params_iter", envir = environment())
    
    # Parallel loop over each subject
    results <- parLapply(cl, unique(dat_pres$SUBJECT_ID), function(id) {
      subdat_pres <- subset(dat_pres, SUBJECT_ID == id)
      y_pres <- subdat_pres$FRI_Score
      x_pres <- as.matrix(subdat_pres[, -c(1, 2)])  # Exclude SUBJECT_ID and FRI_Score
      
      # Calculate linear predictor and mu using the passed-in params_iter
      lin_pred_pres <- as.numeric(params_iter["alpha"]) + x_pres %*% params_iter[-1]
      mu_pres <- inv_logit_vec(lin_pred_pres)
      
      # Construct D_i_pres
      v1 <- mu_pres * (1 - mu_pres)
      D_i_pres <- cbind(v1, sweep(x_pres, 1, v1, "*"))
      
      # Construct R_i_pres based on the correlation structure
      R_i_pres <- diag(nrow(subdat_pres))  # Default: independence
      if (corstr_pres == "exchangeable") {
        R_i_pres <- matrix(rho_pres, nrow = length(y_pres), ncol = length(y_pres))
        diag(R_i_pres) <- 1
      } else if (corstr_pres == "ar1") {
        tooth <- get_tooth_index(subdat_pres)
        absdiff_toothPos <- abs(outer(tooth, tooth, "-"))
        R_i_pres <- rho_pres ^ absdiff_toothPos
      } else if (corstr_pres == "jackknifed") {
        diag(R_i_pres) <- 1
        selected_comb <- paste0(get_tooth_index(subdat_pres), ",", get_zone_index(subdat_pres))
        R_i_pres <- get_corr_submat(master_mat = JK_corr_mat_pres, selected_combinations = selected_comb)
      }
      
      # Calculate the variance matrix V_i_pres
      sqrt_A_i_pres <- diag(as.vector(sqrt(v1)))
      V_i_pres <- sqrt_A_i_pres %*% R_i_pres %*% sqrt_A_i_pres
      
      # Compute psi and Hessian contributions
      psi_i_pres <- t(D_i_pres) %*% ginv(V_i_pres) %*% (y_pres - mu_pres)
      hessian_i_pres <- t(D_i_pres) %*% ginv(V_i_pres) %*% D_i_pres
      
      list(psi_i_pres = psi_i_pres, hessian_i_pres = hessian_i_pres)
    })
    
    # Aggregate psi and Hessian results
    psi_pres <- Reduce("+", lapply(results, function(res) res$psi_i_pres))
    hessian_pres <- Reduce("+", lapply(results, function(res) res$hessian_i_pres))
    
    # Update beta_1, ..., beta_q
    psi_betas <- psi_pres[-1]
    params[iter + 1, -1] <- params[iter, -1] + kappa * ginv(hessian_pres[-1, -1] + psi_betas %*% t(psi_betas)) %*% psi_betas
    
    # Check convergence
    diff <- sum(abs(params[iter + 1, ] - params[iter, ]))
    if (diff < tol) {
      params <- params[1:(iter + 1), , drop = FALSE]
      break
    }
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  coefs <- data.frame(Variables = colnames(params), Estimates = unlist(params[nrow(params), ]))
  
  if (corstr_pres == "jackknifed") {
    results <- list(coefs = coefs, JK_corr_mat_pres = JK_corr_mat_pres, iter = iter)
  } else {
    results <- list(coefs = coefs, rho_pres = rho_pres, iter = iter)
  }
  
  return(results)
}


#########################################
## Improved functions for simulations: ##
#########################################

fit_GEE_combined <- 
  function(dat = IFS_dat, kappa = 1, corstr_pres = "independence", 
           corstr_sev = "independence", maxIter = 10, tol = 1e-7, 
           sep_pres_coef = NA, sep_sev_coef = NA, 
           n_cores = min(16, detectCores() - 1),
           reuse = F, 
           rho_pres = NA, rho_sev = NA, 
           JK_corr_mat_pres = NA, JK_corr_mat_sev = NA,
           init_pres_coef = NA, init_sev_coef = NA) {
    
    #Number of unique severity responses
    n_sev_categ = length(unique((dat[dat$FRI_Score > 0,])$FRI_Score))                  
    
    dat_pres <- dat
    dat_sev <- dat
    
    # Convert FRI_Score to presence responses...FRI=0 is coded as 1 (success) and FRI=1 is coded as 0 (failure or reference) 
    dat_pres$FRI_Score <- ifelse(dat_pres$FRI_Score == 0, 1, 0)
    
    # Filter out severity responses from FRI_Score (1,2 and/or 3)
    dat_sev <- dat_sev %>% filter(FRI_Score > 0)
    
    #######################################
    #Preprocess the responses in dat_sev: #
    #######################################
    
    if(n_sev_categ == 2) 
    {
      #If there are only two categories for severity (1 and 2):
      #prepare a dichotomous response for a logistic regression type GEE model
      #The higher score FRI=2 is coded as 0 ("failure" or "reference"), the lower score FRI=1 is coded as 1 ("success") 
      dat_sev$FRI_Score <- ifelse(dat_sev$FRI_Score == 2, 0, 1)
      
    }else if(n_sev_categ == 3){ 
      
      #If there are three severity categories FRI=1,2,3:
      #make the response factor type to use in a subsequent proportional odds GEE model
      dat_sev$FRI_Score <- as.factor(dat_sev$FRI_Score)
      
    }
    
    #sep_pres_coef <- init_pres_coef
    #sep_sev_coef <- init_sev_coef
    
    if(n_sev_categ==2){
      sep_sev_coef = sep_sev_coef[-c(2)] #get rid of alpha_2|3 term
    } 
    
    flag_missing = F
    
    # Check for missing coefficients
    if (any(is.na(sep_pres_coef))) {
      print("Missing: separate presence coefficients"); flag_missing=T}
    if (any(is.na(sep_sev_coef))) {
      print("Missing: separate severity coefficients"); flag_missing=T}
    
    if(flag_missing == T)
    {
      stop("Function stopped due to missing inputs")
    }
    
    
    # Initialize parameters
    n_beta <- (length(sep_pres_coef) - 1)
    n_params <- (4 + n_beta) # gamma, alpha, alpha_1|2, alpha_2|3, beta_presence
    
    params <- data.frame(matrix(NA, ncol = n_params, nrow = maxIter)) 
    
    colnames(params) <- c("gamma", "alpha", "alpha_1|2", "alpha_2|3", 
                          colnames(dat[,-c(1:2)]))
    
    #Pre-estimate gamma with the ratio of (sum of severity beta coefs):(sum of presence beta coefs)
    gamma = sum(sep_sev_coef[(length(sep_sev_coef) - n_beta + 1) : length(sep_sev_coef)]) / sum(sep_pres_coef[-c(1)]) 
    params$gamma <- gamma
    
    #Get starting values for the iterations
    
    if(is.na(init_pres_coef))
    {
      
      # Get presence estimates' starting values from geeglm
      pres_mod <- geeglm(FRI_Score ~ . -SUBJECT_ID, 
                         data = dat_pres, id = SUBJECT_ID, family = binomial(link = "logit"))
      init_pres_coef  = pres_mod$coefficients
    }
    
    if(is.na(init_sev_coef))
    {
      sev_mod_coef = rep(NA, times = n_params) #initialize with dummy value
      
      if(n_sev_categ==2)
      {
        sev_mod <- geeglm(FRI_Score ~ . -SUBJECT_ID, 
                          data = dat_sev, id = SUBJECT_ID, family = binomial(link = "logit"))
        init_sev_coef  = c(sev_mod$coefficients[1], sev_mod$coefficients[1], sev_mod$coefficients[-c(1)]) #alpha_2|3 gets assigned alpha_1_2
        
      }else if(n_sev_categ==3)
      {
        # Get severity estimates' starting values from polr
        sev_mod <- polr(FRI_Score ~ ., data = dat_sev[,-c(1)]) 
        init_sev_coef = c(sev_mod$zeta, -sev_mod$coefficients) #we will only utilize the intercept estimates
      }
    }
    
    
    params$alpha <- init_pres_coef[1]
    params$`alpha_1|2` <- init_sev_coef[1]
    params$`alpha_2|3` <- init_sev_coef[2]
    params[1, -c(1:4)] <- init_pres_coef[-c(1)] # beta_pres
    
    #####################################
    #### Get correlation parameters #####
    #####################################
    
    
    ##################
    #### Presence ####
    ##################
    
    if(reuse==T) #if rho_pres or JK_corr_mat_pres is provided
    {
      if( !(!is.na(rho_pres) || !any(is.na(JK_corr_mat_pres))) )
      {
        print("reuse=T but correlation parameters were not passed.")
      }
      
    }else{ #if correlation parameters are not provided, estimate them.
      
      rho_pres = 1
      
      # Get correlations
      if (corstr_pres == "exchangeable") {
        rho_pres <- get_rho_pres_exch(dat_pres = dat_pres, init_pres_coef = sep_pres_coef)
      } else if (corstr_pres == "ar1") {
        rho_pres <- get_rho_pres_ar1(dat_pres = dat_pres, init_pres_coef = sep_pres_coef)
      } else if (corstr_pres == "jackknifed") {
        
        JK_corr_mat_pres <- data.frame(diag(4 * 4))
        tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
        tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
        rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
        
        JK_corr_mat_pres <- get_JK_corrmat_pres(dat_pres = dat_pres, init_pres_coef = sep_pres_coef)
        rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
      }
      
      rho_pres <- pmax(pmin(rho_pres, 1), -1)
    }
    
    
    ##################
    #### Severity ####
    ##################
    
    if(reuse==T) #if rho_sev or JK_corr_mat_sev is provided
    {
      if( !( !is.na(rho_sev) || !any(is.na(JK_corr_mat_sev))) )
      {
        print("reuse=T but correlation parameters were not passed.")
      }
      
    }else{ #if correlation parameters are not provided, estimate them.
      
      rho_sev <- 1
      
      if(n_sev_categ==3)
      {
        
        # Get correlations
        if (corstr_sev == "exchangeable") {
          rho_sev <- get_rho_sev_exch(dat_sev = dat_sev, init_sev_coef = sep_sev_coef)
        } else if (corstr_sev == "ar1") {
          rho_sev <- get_rho_sev_ar1(dat_sev = dat_sev, init_sev_coef = sep_sev_coef)
        } else if (corstr_sev == "jackknifed") {
          
          JK_corr_mat_sev <- data.frame(diag(4 * 4))
          tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
          tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
          
          JK_corr_mat_sev <- get_JK_corrmat_sev(dat_sev = dat_sev, init_sev_coef = sep_sev_coef)
          rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
        }
        
      }else if(n_sev_categ==2)
      {
        
        # Get correlations
        if (corstr_sev == "exchangeable") {
          rho_sev <- get_rho_pres_exch(dat_pres = dat_sev, init_pres_coef = sep_sev_coef[-c(2)])
        } else if (corstr_sev == "ar1") {
          rho_sev <- get_rho_pres_ar1(dat_pres = dat_sev, init_pres_coef = sep_sev_coef[-c(2)])
        } else if (corstr_sev == "jackknifed") {
          
          JK_corr_mat_sev <- data.frame(diag(4 * 4))
          tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
          tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
          
          JK_corr_mat_sev <- get_JK_corrmat_pres(dat_pres = dat_sev, init_pres_coef = sep_sev_coef[-c(2)])
          rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
        } 
      }
      
      rho_sev <- pmax(pmin(rho_sev, 1), -1)
      
    }
    
    
    # Set up parallel cluster outside the loop
    cl <- makeCluster(n_cores)
    
    # Load necessary libraries and export objects once before the loop
    clusterEvalQ(cl, {
      library(MASS)
      library(reshape2)
      library(data.table)
      library(dplyr)
      library(Matrix) # For bdiag
    })
    
    clusterExport(cl, varlist = c("dat", "n_beta", "n_params", "get_tooth_index", "get_zone_index", 
                                  "inv_logit_vec", "rho_sev", "JK_corr_mat_sev", "corstr_sev", 
                                  "get_corr_submat", "repeat_column", "transform_to_vec1", "transform_to_vec2", 
                                  "get_multinom_corrmat", "dat_pres", "rho_pres", "JK_corr_mat_pres", 
                                  "corstr_pres"), envir = environment())
    
    for (iter in 1:(maxIter - 1)) {
      
      # Capture the current row of params
      params_iter <- params[iter, ]
      
      clusterExport(cl, "params_iter", envir = environment())
      
      # Parallel loop for each subject using load-balanced parallelization
      results <- parLapplyLB(cl, unique(dat_pres$SUBJECT_ID), function(id) {
        
        subdat <- dat %>% filter(SUBJECT_ID == id)
        
        ###################
        #  Presence model #
        ###################
        
        subdat_pres <- subdat
        subdat_pres$FRI_Score <- ifelse(subdat_pres$FRI_Score == 0, 1, 0)
        y_pres <- subdat_pres$FRI_Score
        x_pres <- as.matrix(subdat_pres[, -c(1:2)])
        
        lin_pred_pres <- as.numeric(params_iter["alpha"]) + x_pres %*% matrix(unlist(params_iter[-c(1:4)]), n_beta, 1)
        mu_pres <- inv_logit_vec(lin_pred_pres)
        
        v1 <- mu_pres * (1 - mu_pres)
        D_i_pres <- cbind(v1, sweep(x_pres, 1, v1, "*"))
        
        # Construct correlation structure for presence
        R_i_pres <- diag(nrow(subdat_pres))
        if (corstr_pres == "exchangeable") {
          R_i_pres <- matrix(rho_pres, nrow = length(y_pres), ncol = length(y_pres))
          diag(R_i_pres) <- 1
        } else if (corstr_pres == "ar1") {
          tooth <- get_tooth_index(subdat_pres)
          absdiff_toothPos <- abs(outer(tooth, tooth, "-"))
          R_i_pres <- rho_pres ^ absdiff_toothPos
        } else if (corstr_pres == "jackknifed") {
          diag(R_i_pres) <- 1
          selected_comb <- paste0(get_tooth_index(subdat_pres), ",", get_zone_index(subdat_pres))
          R_i_pres <- get_corr_submat(master_mat = JK_corr_mat_pres, selected_combinations = selected_comb)
        }
        
        sqrt_A_i_pres <- diag(as.vector(sqrt(v1)))
        V_i_pres <- sqrt_A_i_pres %*% R_i_pres %*% sqrt_A_i_pres
        psi_i_pres <- t(D_i_pres) %*% ginv(V_i_pres) %*% (y_pres - mu_pres)
        hessian_i_pres <- t(D_i_pres) %*% ginv(V_i_pres) %*% D_i_pres
        
        ####################
        ## Severity model ##
        ####################
        
        #initialize with default values:
        psi_i_sev = matrix(numeric((n_beta+2)), ncol = 1) 
        hessian_i_sev = matrix(0, (n_beta+2), (n_beta+2))
        
        if(any(subdat$FRI_Score > 0)) {
          
          subdat_sev <- subdat %>% filter(FRI_Score != 0)
          
          if( all(c(1,2,3) %in% unique(subdat_sev$FRI_Score)) ) #if unique FRI 1,2,3 are present in subdat_sev
          {
            lp_12 <- as.numeric(params_iter["alpha_1|2"]) + 
              as.numeric(params_iter["gamma"]) * 
              as.matrix(subdat_sev[, -c(1:2)]) %*% matrix(unlist(params_iter[-c(1:4)]), n_beta, 1)
            lp_23 <- as.numeric(params_iter["alpha_2|3"]) + as.numeric(params_iter["gamma"]) * 
              as.matrix(subdat_sev[, -c(1:2)]) %*% matrix(unlist(params_iter[-c(1:4)]), n_beta, 1)
            
            cusum_1 <- inv_logit_vec(lp_12)
            cusum_2 <- inv_logit_vec(lp_23)
            cusum_3 <- rep(1, times = length(cusum_1))
            
            pi_1 <- cusum_1
            pi_2 <- cusum_2 - cusum_1
            pi_3 <- 1 - cusum_2
            
            # Variance terms
            v_1 <- pi_1 * (1 - pi_1)
            v_2 <- pi_2 * (1 - pi_2)
            v_3 <- pi_3 * (1 - pi_3)
            
            D_i_T_sev <- matrix(NA, nrow = (n_beta + 2), ncol = 3 * nrow(subdat_sev))
            t1 <- c(sapply(v_1, FUN = transform_to_vec1))
            t2 <- c(sapply(v_2, FUN = transform_to_vec2))
            D_i_T_sev[1, ] <- t1
            D_i_T_sev[2, ] <- t2
            subdat_x <- t(as.matrix(subdat_sev[, -c(1, 2)]))
            repeated_x <- do.call(cbind, lapply(1:ncol(subdat_x), function(i) repeat_column(subdat_x[, i])))
            D_i_T_sev[-c(1, 2), ] <- t1 + t2
            D_i_T_sev[-c(1, 2), ] <- repeated_x * D_i_T_sev[-c(1, 2), ] * as.numeric(params_iter["gamma"]) #as derivative is w.r.t. beta_presence
            
            # Construct R_i_sev based on the correlation structure
            R_i_sev <- diag(nrow(subdat_sev))
            if (corstr_sev == "exchangeable") {
              R_i_sev <- matrix(rho_sev, nrow = length(subdat_sev$FRI_Score), ncol = length(subdat_sev$FRI_Score))
              diag(R_i_sev) <- 1
            } else if (corstr_sev == "ar1") {
              tooth <- get_tooth_index(subdat_sev)
              absdiff_toothPos <- abs(outer(tooth, tooth, "-"))
              R_i_sev <- rho_sev ^ absdiff_toothPos
            } else if (corstr_sev == "jackknifed") {
              diag(R_i_sev) <- 1
              selected_comb <- paste0(get_tooth_index(subdat_sev), ",", get_zone_index(subdat_sev))
              R_i_sev <- get_corr_submat(master_mat = JK_corr_mat_sev, selected_combinations = selected_comb)
            }
            
            list_of_Ms <- lapply(1:length(pi_1), function(i) {
              pi1 <- pi_1[i]
              pi2 <- pi_2[i]
              pi3 <- pi_3[i]
              
              get_multinom_corrmat(pi1, pi2, pi3)
            })
            
            # Block diagonal construction using bdiag from Matrix package
            block_diag_matrix <- bdiag(list_of_Ms)
            block_diag_matrix <- as.matrix(block_diag_matrix)
            block_diag_matrix[block_diag_matrix == 0] <- 1
            
            temp_R_i <- kronecker(R_i_sev, matrix(1, 3, 3))
            R_i_sev2 <- temp_R_i * block_diag_matrix
            R_i_sev2[is.infinite(R_i_sev2)] <- 1e10  # Replace Inf with a large number
            R_i_sev2[is.nan(R_i_sev2)] <- 0  # Replace NaN with 0
            
            sqrt_A_i_sev <- diag(sqrt(c(v_1, v_2, v_3)))
            V_i_sev <- sqrt_A_i_sev %*% R_i_sev2 %*% sqrt_A_i_sev
            
            z1 <- ifelse(subdat_sev$FRI_Score == 1, 1, 0)
            z2 <- ifelse(subdat_sev$FRI_Score == 2, 1, 0)
            z3 <- ifelse(subdat_sev$FRI_Score == 3, 1, 0)
            
            resid_sev <- c((z1 - pi_1), (z2 - pi_2), (z3 - pi_3))
            psi_i_sev <- D_i_T_sev %*% ginv(V_i_sev) %*% as.matrix(resid_sev)
            hessian_i_sev <- D_i_T_sev %*% ginv(V_i_sev) %*% t(D_i_T_sev)
            
          }else if( all(c(1,2) %in% unique(subdat_sev$FRI_Score)) ) #If the only unique FRI severity scores are 1 & 2
          {
            #Code "FRI=1" as 1 (success) and "FRI=2" as 0 (failure/reference level):
            subdat_sev$FRI_Score <- ifelse(subdat_sev$FRI_Score == 2, 0, 1)
            
            y_sev <- subdat_sev$FRI_Score
            x_sev <- as.matrix(subdat_sev[, -c(1, 2)])  # Exclude SUBJECT_ID and FRI_Score
            
            # Calculate linear predictor and mu using the passed-in params_iter
            lin_pred_sev <- as.numeric(params_iter$`alpha_1|2`) + x_sev %*% matrix(unlist(params_iter[-c(1,2,3,4)]), n_beta, 1)
            mu_sev <- inv_logit_vec(lin_pred_sev)
            
            # Construct D_i_sev
            v1 <- mu_sev * (1 - mu_sev)
            D_i_sev <- cbind(v1, sweep(x_sev, 1, v1, "*"))
            
            # Construct R_i_sev based on the correlation structure
            R_i_sev <- diag(nrow(subdat_sev))  # Default: independence
            if (corstr_sev == "exchangeable") {
              R_i_sev <- matrix(rho_sev, nrow = length(y_sev), ncol = length(y_sev))
              diag(R_i_sev) <- 1
            } else if (corstr_sev == "ar1") {
              tooth <- get_tooth_index(subdat_sev)
              absdiff_toothPos <- abs(outer(tooth, tooth, "-"))
              R_i_sev <- rho_sev ^ absdiff_toothPos
            } else if (corstr_sev == "jackknifed") {
              diag(R_i_sev) <- 1
              selected_comb <- paste0(get_tooth_index(subdat_sev), ",", get_zone_index(subdat_sev))
              R_i_sev <- get_corr_submat(master_mat = JK_corr_mat_sev, selected_combinations = selected_comb)
            }
            
            # Calculate the variance matrix V_i_pres
            sqrt_A_i_sev <- diag(as.vector(sqrt(v1)))
            V_i_sev <- sqrt_A_i_sev %*% R_i_sev %*% sqrt_A_i_sev
            
            # Compute psi and Hessian contributions
            psi_i_sev <- t(D_i_sev) %*% ginv(V_i_sev) %*% (y_sev - mu_sev)
            hessian_i_sev <- t(D_i_sev) %*% ginv(V_i_sev) %*% D_i_sev
            
            #Adjustments for the alpha_2|3 term:
            
            #Add elements for alpha_1|2 and set as zero
            psi_i_sev = matrix(c(psi_i_sev[1], 0, psi_i_sev[-c(1)]), ncol=1)
            
            #insert matrix row and column of zeroes for alpha_2|3
            hessian_i_updated = matrix(0, (n_beta+2), (n_beta+2))
            hessian_i_updated[1:2,1:2] = hessian_i_sev[1:2,1:2]
            hessian_i_updated[4:(n_beta+2), 1:2] = hessian_i_sev[3:(n_beta+1), 1:2]
            hessian_i_updated[1:2, 4:(n_beta+2)] = hessian_i_sev[1:2, 3:(n_beta+1)]
            hessian_i_updated[4:(n_beta+2),4:(n_beta+2)] = hessian_i_sev[3:(n_beta+1), 3:(n_beta+1)]
            
            hessian_i_sev = hessian_i_updated
            
          }
          
        }
        
        output_list = list(psi_i_pres = psi_i_pres, 
                           hessian_i_pres = hessian_i_pres,
                           psi_i_sev = psi_i_sev, 
                           hessian_i_sev = hessian_i_sev)
        
        return(output_list)
      })
      
      # Aggregate results
      psi_pres <- Reduce("+", lapply(results, function(res) res$psi_i_pres))
      hessian_pres <- Reduce("+", lapply(results, function(res) res$hessian_i_pres))
      psi_sev <- Reduce("+", lapply(results, function(res) res$psi_i_sev))
      hessian_sev <- Reduce("+", lapply(results, function(res) res$hessian_i_sev))
      
      # Update beta coefficients
      psi_betas <- psi_pres[-1] + psi_sev[-c(1:2)]
      params[iter + 1, -c(1:4)] <- params[iter, -c(1:4)] + 
        kappa * ginv(hessian_pres[-1, -1] + hessian_sev[-c(1:2), -c(1:2)] + psi_betas %*% t(psi_betas)) %*% psi_betas
      
      
      diff <- sum(abs(params[(iter + 1), ] - params[iter, ]))
      
      # diff <- sum(abs(params[(iter + 1), -c(4)] - params[iter, -c(4)]))
      # 
      # if(n_sev_categ==3)
      # {
      #   diff <- sum(abs(params[(iter + 1), ] - params[iter, ]))
      # }
      
      if (diff < tol) {
        params <- params[1:(iter + 1), ]
        break
      }
      
    }
    
    # Stop the cluster after all iterations
    stopCluster(cl)
    
    pres_coefs_final <- data.frame(Variables = colnames(params), Estimates = as.numeric(params[nrow(params), ]))
    sev_coefs_final <- pres_coefs_final
    sev_coefs_final[-c(1:4), 2] <- as.numeric(pres_coefs_final$Estimates[1]) * pres_coefs_final$Estimates[-c(1:4)]
    
    results <- list(pres_coefs = pres_coefs_final, rho_pres = rho_pres, 
                    sev_coefs = sev_coefs_final, rho_sev = rho_sev, iter=iter)
    if(corstr_sev == "jackknifed") {
      results <- list(pres_coefs = pres_coefs_final, JK_corrmat_pres = JK_corr_mat_pres, 
                      sev_coefs = sev_coefs_final, JK_corrmat_sev = JK_corr_mat_sev, iter=iter)
    }
    
    return(results)
    
  }



fit_GEE_severity <- function(dat = IFS_dat, kappa = 0.25, corstr_sev = "independence", 
                             maxIter = 10, tol = 1e-7, n_cores = min(16, detectCores() - 1),
                             reuse=F, rho_sev=NA, JK_corr_mat_sev=NA,
                             init_sev_coef=NA
) {
  
  #Number of unique severity responses
  n_sev_categ = length(unique((dat[dat$FRI_Score > 0,])$FRI_Score))   
  
  dat_sev <- dat
  dat_sev <- dat_sev %>% filter(FRI_Score > 0) # Convert FRI_Score to W_s (1, 2, 3)
  
  if(n_sev_categ == 2) 
  {
    #If there are only two categories for severity (1 and 2):
    #prepare a dichotomous response for a logistic regression type GEE model
    #The higher score FRI=2 is coded as 0 ("failure" or "reference"), the lower score FRI=1 is coded as 1 ("success") 
    dat_sev$FRI_Score <- ifelse(dat_sev$FRI_Score == 2, 0, 1)
    
  }else if(n_sev_categ == 3){ 
    
    #If there are three severity categories FRI=1,2,3:
    #make the response factor type to use in a subsequent proportional odds GEE model
    dat_sev$FRI_Score <- as.factor(dat_sev$FRI_Score)
  }
  
  
  sev_coef = init_sev_coef
  
  if(is.na(init_sev_coef))
  {
    if(n_sev_categ==2)
    {
      sev_mod <- geeglm(FRI_Score ~ . -SUBJECT_ID, 
                        data = dat_sev, id = SUBJECT_ID, family = binomial(link = "logit"))
      sev_coef  = c(sev_mod$coefficients[1], sev_mod$coefficients[1], sev_mod$coefficients[-c(1)]) #alpha_2|3 gets assigned 0
      
    }else if(n_sev_categ==3)
    {
      # Get severity estimates' starting values from polr
      
      # Try to fit the model with automatic starting values
      sev_mod <- tryCatch({
        polr(FRI_Score ~ ., data = dat_sev[, -c(1)])
      }, error = function(e) {
        # Check if the error is due to failure of finding starting values
        if(grepl("attempt to find suitable starting values failed", e$message)) {
          message("Default starting values failed. Trying with specified starting values...")
          # Retry with specified starting values
          temp_obj = polr(FRI_Score ~ ., data = dat_sev[, -c(1)], 
                          start = c(-2, -0.38, -0.13, -0.08, 0.18, -0.34, -0.74, -0.96, 1, 2))
          message("Converged with specified starting values.")
          return(temp_obj)
        } else {
          # Re-throw the error if it is not related to starting values
          stop(e)
        }
      })
       
      sev_coef = c(sev_mod$zeta, -sev_mod$coefficients) 
    }
  }
  
  
  # Initialize the parameters matrix
  n_beta <- (length(sev_coef) - 2)
  n_params <- (2 + n_beta)
  
  
  params <- data.frame(matrix(NA, ncol = n_params, nrow = maxIter)) 
  colnames(params) <- c("alpha_1|2", "alpha_2|3", names(sev_coef)[-c(1:2)])
  
  params$`alpha_1|2` <- sev_coef[1]
  params$`alpha_2|3` <- sev_coef[2]
  params[1, -c(1:2)] <- sev_coef[-c(1:2)]
  
  ########################################
  ###### Get correlation parameters ######
  ########################################
  
  if(reuse==T) #if rho_sev or JK_corr_mat_sev is provided
  {
    if( !( !is.na(rho_sev) || !any(is.na(JK_corr_mat_sev))) )
    {
      print("reuse=T but correlation parameters were not passed.")
    }
    
  }else{ #if correlation parameters are not provided, estimate them.
    
    rho_sev <- 1
    
    if(corstr_sev %in% c("ar1", "exchangeable"))
    {
      # Initialize rho_sev based on correlation structure
      
      if(n_sev_categ==3)
      {
        # Get correlations
        if (corstr_sev == "exchangeable") {
          rho_sev <- get_rho_sev_exch(dat_sev = dat_sev, init_sev_coef = sev_coef)
        } else if (corstr_sev == "ar1") {
          rho_sev <- get_rho_sev_ar1(dat_sev = dat_sev, init_sev_coef = sev_coef)
        } 
        
      }else if(n_sev_categ==2)
      {
        
        # Get correlations
        if (corstr_sev == "exchangeable") {
          rho_sev <- get_rho_pres_exch(dat_pres = dat_sev, init_pres_coef = sev_coef[-c(2)])
        } else if (corstr_sev == "ar1") {
          rho_sev <- get_rho_pres_ar1(dat_pres = dat_sev, init_pres_coef = sev_coef[-c(2)])
        } 
      }
      
      rho_sev <- pmax(pmin(rho_sev, 1), -1)
      
    }else if(corstr_sev=="jackknifed") #need JK_corr_mat_sev
    {
      
      JK_corr_mat_sev <- data.frame(diag(4 * 4))
      tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
      tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
      rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
      
      if(n_sev_categ==3)
      {
        JK_corr_mat_sev <- get_JK_corrmat_sev(dat_sev = dat_sev, init_sev_coef = sev_coef)
        rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
        
      }else if(n_sev_categ==2)
      {
        JK_corr_mat_sev <- get_JK_corrmat_pres(dat_pres = dat_sev, init_pres_coef = sev_coef[-c(2)])
        rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
      }
    }
  }
  
  # Set up parallel cluster
  cl <- makeCluster(n_cores)
  
  # Load the necessary libraries on each worker
  clusterEvalQ(cl, {
    library(MASS)
    library(reshape2)
    library(data.table)
    library(dplyr)
    library(Matrix)  # For bdiag
  })
  
  # Export necessary objects and functions to the parallel environment
  clusterExport(cl, varlist = c("dat", "n_beta", "n_params", "get_tooth_index", "get_zone_index", 
                                "inv_logit_vec", "rho_sev", "JK_corr_mat_sev", 
                                "corstr_sev", "get_corr_submat", "repeat_column", 
                                "transform_to_vec1", "transform_to_vec2", "get_multinom_corrmat"), 
                envir = environment())
  
  for(iter in 1:(maxIter - 1)) {
    #psi_sev <- numeric((2 + n_beta))
    #hessian_sev <- matrix(0, (2 + n_beta), (2 + n_beta))
    
    # Export params for the current iteration
    params_iter <- params[iter, ]
    clusterExport(cl, "params_iter", envir = environment())
    
    # Parallel loop over each subject
    results <- parLapply(cl, unique(dat_sev$SUBJECT_ID), function(id) {
      
      subdat <- dat %>% filter(SUBJECT_ID == id)
      subdat_sev <- subdat %>% filter(FRI_Score != 0)
      
      #initialize with default values:
      psi_i_sev = matrix(numeric((n_beta+2)), ncol = 1) 
      hessian_i_sev = matrix(0, (n_beta+2), (n_beta+2))
      
      has_nonzeroFRI <- any(subdat_sev$FRI_Score > 0)
      
      if (has_nonzeroFRI) {
        
        if( all(c(1,2,3) %in% unique(subdat_sev$FRI_Score)) ) #if unique FRI 1,2,3 are present in subdat_sev
        {
          lp_12 <- as.numeric(params_iter["alpha_1|2"]) + 
            as.matrix(subdat_sev[, -c(1:2)]) %*% matrix(unlist(params_iter[-c(1:2)]), n_beta, 1)
          lp_23 <- as.numeric(params_iter["alpha_2|3"]) + 
            as.matrix(subdat_sev[, -c(1:2)]) %*% matrix(unlist(params_iter[-c(1:2)]), n_beta, 1)
          
          cusum_1 <- inv_logit_vec(lp_12)
          cusum_2 <- inv_logit_vec(lp_23)
          cusum_3 <- rep(1, times = length(cusum_1))
          
          pi_1 <- cusum_1
          pi_2 <- cusum_2 - cusum_1
          pi_3 <- 1 - cusum_2
          
          # Variance terms
          v_1 <- pi_1 * (1 - pi_1)
          v_2 <- pi_2 * (1 - pi_2)
          v_3 <- pi_3 * (1 - pi_3)
          
          D_i_T_sev <- matrix(NA, nrow = (n_beta + 2), ncol = 3 * nrow(subdat_sev))
          t1 <- c(sapply(v_1, FUN = transform_to_vec1))
          t2 <- c(sapply(v_2, FUN = transform_to_vec2))
          D_i_T_sev[1, ] <- t1
          D_i_T_sev[2, ] <- t2
          subdat_x <- t(as.matrix(subdat_sev[, -c(1, 2)]))
          repeated_x <- do.call(cbind, lapply(1:ncol(subdat_x), function(i) repeat_column(subdat_x[, i])))
          D_i_T_sev[-c(1, 2), ] <- t1 + t2
          D_i_T_sev[-c(1, 2), ] <- repeated_x * D_i_T_sev[-c(1, 2), ] #as derivative is w.r.t. beta_presence
          
          # Construct R_i_sev based on the correlation structure
          R_i_sev <- diag(nrow(subdat_sev))
          if (corstr_sev == "exchangeable") {
            R_i_sev <- matrix(rho_sev, nrow = length(subdat_sev$FRI_Score), ncol = length(subdat_sev$FRI_Score))
            diag(R_i_sev) <- 1
          } else if (corstr_sev == "ar1") {
            tooth <- get_tooth_index(subdat_sev)
            absdiff_toothPos <- abs(outer(tooth, tooth, "-"))
            R_i_sev <- rho_sev ^ absdiff_toothPos
          } else if (corstr_sev == "jackknifed") {
            diag(R_i_sev) <- 1
            selected_comb <- paste0(get_tooth_index(subdat_sev), ",", get_zone_index(subdat_sev))
            R_i_sev <- get_corr_submat(master_mat = JK_corr_mat_sev, selected_combinations = selected_comb)
          }
          
          list_of_Ms <- lapply(1:length(pi_1), function(i) {
            pi1 <- pi_1[i]
            pi2 <- pi_2[i]
            pi3 <- pi_3[i]
            
            get_multinom_corrmat(pi1, pi2, pi3)
          })
          
          # Block diagonal construction using bdiag from Matrix package
          block_diag_matrix <- bdiag(list_of_Ms)
          block_diag_matrix <- as.matrix(block_diag_matrix)
          block_diag_matrix[block_diag_matrix == 0] <- 1
          
          temp_R_i <- kronecker(R_i_sev, matrix(1, 3, 3))
          R_i_sev2 <- temp_R_i * block_diag_matrix
          R_i_sev2[is.infinite(R_i_sev2)] <- 1e10  # Replace Inf with a large number
          R_i_sev2[is.nan(R_i_sev2)] <- 0  # Replace NaN with 0
          
          sqrt_A_i_sev <- diag(sqrt(c(v_1, v_2, v_3)))
          V_i_sev <- sqrt_A_i_sev %*% R_i_sev2 %*% sqrt_A_i_sev
          
          z1 <- ifelse(subdat_sev$FRI_Score == 1, 1, 0)
          z2 <- ifelse(subdat_sev$FRI_Score == 2, 1, 0)
          z3 <- ifelse(subdat_sev$FRI_Score == 3, 1, 0)
          
          resid_sev <- c((z1 - pi_1), (z2 - pi_2), (z3 - pi_3))
          psi_i_sev <- D_i_T_sev %*% ginv(V_i_sev) %*% as.matrix(resid_sev)
          hessian_i_sev <- D_i_T_sev %*% ginv(V_i_sev) %*% t(D_i_T_sev)
          
        }else if( all(c(1,2) %in% unique(subdat_sev$FRI_Score)) ) #If the only unique FRI severity scores are 1 & 2
        {
          #Code "FRI=1" as 1 (success) and "FRI=2" as 0 (failure/reference level):
          subdat_sev$FRI_Score <- ifelse(subdat_sev$FRI_Score == 2, 0, 1)
          
          y_sev <- subdat_sev$FRI_Score
          x_sev <- as.matrix(subdat_sev[, -c(1, 2)])  # Exclude SUBJECT_ID and FRI_Score
          
          # Calculate linear predictor and mu using the passed-in params_iter
          lin_pred_sev <- as.numeric(params_iter$`alpha_1|2`) + x_sev %*% matrix(unlist(params_iter[-c(1,2)]), n_beta, 1)
          mu_sev <- inv_logit_vec(lin_pred_sev)
          
          # Construct D_i_sev
          v1 <- mu_sev * (1 - mu_sev)
          D_i_sev <- cbind(v1, sweep(x_sev, 1, v1, "*"))
          
          # Construct R_i_sev based on the correlation structure
          R_i_sev <- diag(nrow(subdat_sev))  # Default: independence
          if (corstr_sev == "exchangeable") {
            R_i_sev <- matrix(rho_sev, nrow = length(y_sev), ncol = length(y_sev))
            diag(R_i_sev) <- 1
          } else if (corstr_sev == "ar1") {
            tooth <- get_tooth_index(subdat_sev)
            absdiff_toothPos <- abs(outer(tooth, tooth, "-"))
            R_i_sev <- rho_sev ^ absdiff_toothPos
          } else if (corstr_sev == "jackknifed") {
            diag(R_i_sev) <- 1
            selected_comb <- paste0(get_tooth_index(subdat_sev), ",", get_zone_index(subdat_sev))
            R_i_sev <- get_corr_submat(master_mat = JK_corr_mat_sev, selected_combinations = selected_comb)
          }
          
          # Calculate the variance matrix V_i_pres
          sqrt_A_i_sev <- diag(as.vector(sqrt(v1)))
          V_i_sev <- sqrt_A_i_sev %*% R_i_sev %*% sqrt_A_i_sev
          
          # Compute psi and Hessian contributions
          psi_i_sev <- t(D_i_sev) %*% ginv(V_i_sev) %*% (y_sev - mu_sev)
          hessian_i_sev <- t(D_i_sev) %*% ginv(V_i_sev) %*% D_i_sev
          
          #Adjustments for the alpha_2|3 term:
          
          #Add elements for alpha_1|2 and set as zero
          psi_i_sev = matrix(c(psi_i_sev[1], 0, psi_i_sev[-c(1)]), ncol=1)
          
          #insert matrix row and column of zeroes for alpha_2|3
          hessian_i_updated = matrix(0, (n_beta+2), (n_beta+2))
          hessian_i_updated[1:2,1:2] = hessian_i_sev[1:2,1:2]
          hessian_i_updated[4:(n_beta+2), 1:2] = hessian_i_sev[3:(n_beta+1), 1:2]
          hessian_i_updated[1:2, 4:(n_beta+2)] = hessian_i_sev[1:2, 3:(n_beta+1)]
          hessian_i_updated[4:(n_beta+2),4:(n_beta+2)] = hessian_i_sev[3:(n_beta+1), 3:(n_beta+1)]
          
          hessian_i_sev = hessian_i_updated
          
        }
        
        output_list = list(psi_i_sev = psi_i_sev, 
                           hessian_i_sev = hessian_i_sev)
        
        return(output_list)
        
      }
    })
    
    # Aggregate psi and Hessian results
    psi_sev <- Reduce("+", lapply(results, function(res) res$psi_i_sev))
    hessian_sev <- Reduce("+", lapply(results, function(res) res$hessian_i_sev))
    
    # Update parameters
    psi_betas <- psi_sev[-c(1:2)]
    params[iter + 1, -c(1:2)] <- params[iter, -c(1:2)] +
      kappa * ginv(hessian_sev[-c(1:2), -c(1:2)] + psi_betas %*% t(psi_betas)) %*% psi_betas
    
    # Check convergence
    diff <- sum(abs(params[(iter + 1), ] - params[iter, ]))
    if (diff < tol) {
      params <- params[1:(iter + 1), ]
      break
    }
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  coefs <- data.frame(Variables = colnames(params), Estimates = unlist(params[nrow(params), ]))
  
  if (corstr_sev == "jackknifed") {
    results <- list(coefs = coefs, JK_corr_mat_sev = JK_corr_mat_sev, iter = iter)
  } else {
    results <- list(coefs = coefs, rho_sev = rho_sev, iter = iter)
  }
  
  return(results)
}


################################
#### Get simulated dataset: ####
################################

update_MC_dataset = function(dataset=MC_dataset, seed=mc_seed)
{
  #Update the Y in dataset:
  #####################################################
  ##  Set up Sigma, the cluster correlation matrix:  ##
  #####################################################
  
  #Initialize:
  Sigma_P = diag(16) 
  
  if(corstr_pres=="independence")
  { 
    #With independence cluster correlation structure
    Sigma_P = diag(16)
  }else if(corstr_pres=="exchangeable")
  { 
    #With exchangeable cluster correlation structure
    Sigma_P = matrix(exch_rho_pres, 16, 16); diag(Sigma_P) = 1
  }
  if(corstr_pres=="jackknifed")
  {
    
    #Create master matrix of jackknifed correlations:
    JK_corr_mat_pres <- data.frame(diag(4 * 4))
    tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
    tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
    rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
    
    JK_corr_mat_pres <- get_JK_corrmat_pres(dat_pres = dataset, 
                                            init_pres_coef = true_params[-c(1:4)])
    
    rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
    
    
    #Assign the Sigma_P as this master correlation matrix:
    Sigma_P = JK_corr_mat_pres
    
  }
  
  set.seed(seed)
  
  #Generate the latent MVN variables:
  Lp = mvrnorm(length(unique(dataset$SUBJECT_ID)), numeric(16), Sigma = Sigma_P) 
  #output is a matrix, where the i-th row has elements for the i-th cluster.
  
  
  ######################################################
  ## Transform Lp --> F(Lp) ~ U(0,1): ##################
  ######################################################
  
  #Initialize:
  Up = Lp 
  
  #Then the marginal CDF for each element in Lp is the N(0,1):
  Up = apply(Lp, MARGIN = 1, pnorm)
  
  ############################################################
  ## ASSIGN PRESENCE CATEGORIES:                            ##
  ############################################################
  ## Calculate the cumulative probs.                        ##
  ## p0 = exp(alpha + X beta) / (1 + exp(alpha + X beta) )  ##
  ## p1 = 1 / exp(alpha + X beta)                           ##
  ## If Up < p0, then assign Y=0. Else Y=1                  ##
  ############################################################
  
  alpha = true_params["alpha"]
  beta_p = true_params[-c(1:4)]
  
  X = dataset[,-c(1,2)]
  
  linpred = alpha + as.matrix(X) %*% matrix(beta_p, nrow= (ncol(dataset)-2), ncol=1)
  
  
  p0 = matrix(as.vector(inv_logit_vec(linpred)), 
              byrow = T, 
              nrow=16, 
              ncol=length(unique(dataset$SUBJECT_ID)))
  
  bool_Wp_pres = (Up <= p0)
  
  sim_Wp_pres = ifelse(bool_Wp_pres, 0, 9) 
  #TRUE values get assigned category 0 (no fluorosis), FALSE values get assigned to dummy category 9 (fluorosis) 
  
  #Unlist sim_Wp_pres elements row-wise and populate the FRI_Score column in dataset. 
  #This is an intermediate step; this will later be re-updated with severity categories
  dataset[,"FRI_Score"] = as.vector(t(sim_Wp_pres)) 
  
  ##################################
  ###### SEVERITY CATEGORIES #######
  ##################################
  
  Sigma_S = diag(16) #initialize
  
  if(corstr_sev=="independence")
  { 
    #With independence cluster correlation structure
    Sigma_S = diag(16)
    
  }else if(corstr_sev=="exchangeable")
  { 
    #With exchangeable cluster correlation structure
    Sigma_S = matrix(exch_rho_sev, 16, 16); diag(Sigma_S) = 1
    
  }else if(corstr_sev=="jackknifed")
  {
    #With jackknifed cluster correlation structure
    
    JK_corr_mat_sev <- data.frame(diag(4 * 4))
    tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
    tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
    rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
    
    alphas = true_params[c(3,4)]
    betas = true_params[-c(1:4)] * true_params[1]
    init_sev_coef = c(alphas, betas)
    
    JK_corr_mat_sev <- get_JK_corrmat_sev(dat_sev = dataset, 
                                          init_sev_coef = init_sev_coef)
    
    rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
    
    Sigma_S = JK_corr_mat_sev 
  }
  
  set.seed(seed)
  
  #Generate latent MVN variables:
  Ls = mvrnorm(length(unique(dataset$SUBJECT_ID)), numeric(16), Sigma = Sigma_S) 
  #the output is a matrix; the i-th row has elements for the i-th cluster.
  
  ######################################################
  ## Transform Ls --> F(Ls) ~ U(0,1): ##################
  ######################################################
  
  #Initialize:
  Us = Ls 
  
  #The the marginal CDF for each Lp is of the N(0,1)
  Us = apply(Ls, MARGIN = 1, pnorm)
  
  ############################################################
  ## ASSIGN SEVERITY CATEGORIES:                            ##
  ############################################################
  ## Calculate the cumulative probs.                        ##
  ## p1 = P(Ws<=1)                                          ##
  ## p2 = P(Ws<=2) - P(Ws<=1)                               ##
  ## p3 = 1 - P(Ws<=2)                                      ##
  ## where P(Ws<=l) = inv_logit(alpha_{l|l+1} + beta_S'X)   ##
  ## If  Wp=0: no updation                                  ##
  ##                                                        ##
  ## For Wp=1: if Us in (0, p1), then Y=1                   ##
  ##           if Us in (p1, p1+p2)  then Y=2               ##
  ##           if Us in (p1+p2, 1) then Y=3                 ##
  ############################################################
  
  alpha = true_params[c("alpha_1|2", "alpha_2|3")]
  gamma = true_params["gamma"]; beta_p = true_params[-c(1:4)]; beta_s = gamma*beta_p
  
  
  X = dataset[,-c(1,2)]
  
  linpred12 = alpha[1] + as.matrix(X) %*% matrix(beta_s, nrow= (ncol(dataset)-2), ncol=1)
  linpred23 = alpha[2] + as.matrix(X) %*% matrix(beta_s, nrow= (ncol(dataset)-2), ncol=1)
  
  cusum1 = matrix(as.vector(inv_logit_vec(linpred12)), byrow = T, nrow=16, ncol=length(unique(dataset$SUBJECT_ID)))
  cusum2 = matrix(as.vector(inv_logit_vec(linpred23)), byrow = T, nrow=16, ncol=length(unique(dataset$SUBJECT_ID)))
  
  p1 = cusum1
  p2 = cusum2 - cusum1
  p3 = 1 - cusum2
  
  
  bool_Y1_pres = ( Us <= p1 )
  bool_Y2_pres = ( Us > p1 & Us <= (p1+p2) )
  bool_Y3_pres = ( Us > p1+p2 )
  
  sim_Y_pres = sim_Wp_pres
  
  #If bool_Wp_pres==T i.e. fluorosis present and bool_Y1_pres==T, assign to category 1
  sim_Y_pres = ifelse(!bool_Wp_pres & bool_Y1_pres, 1, sim_Y_pres) 
  
  #If bool_Wp_pres==T i.e. fluorosis present and bool_Y1_pres==T, assign to category 2
  sim_Y_pres = ifelse(!bool_Wp_pres & bool_Y2_pres, 2, sim_Y_pres) 
  
  #If bool_Wp_pres==T i.e. fluorosis present and bool_Y1_pres==T, assign to category 3
  sim_Y_pres = ifelse(!bool_Wp_pres & bool_Y3_pres, 3, sim_Y_pres) 
  
  #Finally, update the dataset with the generated responses  
  dataset[,"FRI_Score"] = as.vector(t(sim_Y_pres)) #unlist rowwise and populate the FRI_Score column in dataset
  
  
  table(dataset$FRI_Score)
  
  return(dataset)
}

update_MC_dataset2 = function(dataset=MC_dataset, seed=mc_seed)
{
  #Update the Y in dataset:
  #####################################################
  ##  Set up Sigma, the cluster correlation matrix:  ##
  #####################################################
  
  #Initialize:
  Sigma_P = diag(16) 
  
  # if(corstr_pres=="independence")
  # { 
  #   #With independence cluster correlation structure
  #   Sigma_P = diag(16)
  # }else if(corstr_pres=="exchangeable")
  # { 
  #With exchangeable cluster correlation structure
  Sigma_P = matrix(exch_rho_pres, 16, 16); diag(Sigma_P) = 1
  # }
  # if(corstr_pres=="jackknifed")
  # {
  #   
  #   #Create master matrix of jackknifed correlations:
  #   JK_corr_mat_pres <- data.frame(diag(4 * 4))
  #   tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
  #   tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
  #   rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
  #   
  #   JK_corr_mat_pres <- get_JK_corrmat_pres(dat_pres = dataset, 
  #                                           init_pres_coef = true_params[-c(1:4)])
  #   
  #   rownames(JK_corr_mat_pres) <- colnames(JK_corr_mat_pres) <- tooth_zone_comb$combination
  #   
  #   
  #   #Assign the Sigma_P as this master correlation matrix:
  #   Sigma_P = JK_corr_mat_pres
  #   
  # }
  
  set.seed(seed)
  
  #Generate the latent MVN variables:
  Lp = mvrnorm(length(unique(dataset$SUBJECT_ID)), numeric(16), Sigma = Sigma_P) 
  #output is a matrix, where the i-th row has elements for the i-th cluster.
  
  
  ######################################################
  ## Transform Lp --> F(Lp) ~ U(0,1): ##################
  ######################################################
  
  #Initialize:
  Up = Lp 
  
  #Then the marginal CDF for each element in Lp is the N(0,1):
  Up = apply(Lp, MARGIN = 1, pnorm)
  
  ############################################################
  ## ASSIGN PRESENCE CATEGORIES:                            ##
  ############################################################
  ## Calculate the cumulative probs.                        ##
  ## p0 = exp(alpha + X beta) / (1 + exp(alpha + X beta) )  ##
  ## p1 = 1 / exp(alpha + X beta)                           ##
  ## If Up < p0, then assign Y=0. Else Y=1                  ##
  ############################################################
  
  alpha = true_params["alpha"]
  beta_p = true_params[-c(1:4)]
  
  X = dataset[,-c(1,2)]
  
  linpred = alpha + as.matrix(X) %*% matrix(beta_p, nrow= (ncol(dataset)-2), ncol=1)
  
  
  p0 = matrix(as.vector(inv_logit_vec(linpred)), 
              byrow = T, 
              nrow=16, 
              ncol=length(unique(dataset$SUBJECT_ID)))
  
  bool_Wp_pres = (Up <= p0)
  
  sim_Wp_pres = ifelse(bool_Wp_pres, 0, 9) 
  #TRUE values get assigned category 0 (no fluorosis), FALSE values get assigned to dummy category 9 (fluorosis) 
  
  #Unlist sim_Wp_pres elements row-wise and populate the FRI_Score column in dataset. 
  #This is an intermediate step; this will later be re-updated with severity categories
  dataset[,"FRI_Score"] = as.vector(t(sim_Wp_pres)) 
  
  ##################################
  ###### SEVERITY CATEGORIES #######
  ##################################
  
  Sigma_S = diag(16) #initialize
  
  # if(corstr_sev=="independence")
  # { 
  #   #With independence cluster correlation structure
  #   Sigma_S = diag(16)
  #   
  # }else if(corstr_sev=="exchangeable")
  # { 
  #With exchangeable cluster correlation structure
  Sigma_S = matrix(exch_rho_sev, 16, 16); diag(Sigma_S) = 1
  
  # }else if(corstr_sev=="jackknifed")
  # {
  #   #With jackknifed cluster correlation structure
  #   
  #   JK_corr_mat_sev <- data.frame(diag(4 * 4))
  #   tooth_zone_comb <- expand.grid(Tooth = c(7:10), Zone = c("C", "M", "I", "O"))
  #   tooth_zone_comb$combination <- paste(tooth_zone_comb$Tooth, tooth_zone_comb$Zone, sep = ",")
  #   rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
  #   
  #   alphas = true_params[c(3,4)]
  #   betas = true_params[-c(1:4)] * true_params[1]
  #   init_sev_coef = c(alphas, betas)
  #   
  #   JK_corr_mat_sev <- get_JK_corrmat_sev(dat_sev = dataset, 
  #                                         init_sev_coef = init_sev_coef)
  #   
  #   rownames(JK_corr_mat_sev) <- colnames(JK_corr_mat_sev) <- tooth_zone_comb$combination
  #   
  #   Sigma_S = JK_corr_mat_sev 
  # }
  
  set.seed(seed)
  
  #Generate latent MVN variables:
  Ls = mvrnorm(length(unique(dataset$SUBJECT_ID)), numeric(16), Sigma = Sigma_S) 
  #the output is a matrix; the i-th row has elements for the i-th cluster.
  
  ######################################################
  ## Transform Ls --> F(Ls) ~ U(0,1): ##################
  ######################################################
  
  #Initialize:
  Us = Ls 
  
  #The the marginal CDF for each Lp is of the N(0,1)
  Us = apply(Ls, MARGIN = 1, pnorm)
  
  ############################################################
  ## ASSIGN SEVERITY CATEGORIES:                            ##
  ############################################################
  ## Calculate the cumulative probs.                        ##
  ## p1 = P(Ws<=1)                                          ##
  ## p2 = P(Ws<=2) - P(Ws<=1)                               ##
  ## p3 = 1 - P(Ws<=2)                                      ##
  ## where P(Ws<=l) = inv_logit(alpha_{l|l+1} + beta_S'X)   ##
  ## If  Wp=0: no updation                                  ##
  ##                                                        ##
  ## For Wp=1: if Us in (0, p1), then Y=1                   ##
  ##           if Us in (p1, p1+p2)  then Y=2               ##
  ##           if Us in (p1+p2, 1) then Y=3                 ##
  ############################################################
  
  alpha = true_params[c("alpha_1|2", "alpha_2|3")]
  gamma = true_params["gamma"]; beta_p = true_params[-c(1:4)]; beta_s = gamma*beta_p
  
  
  X = dataset[,-c(1,2)]
  
  linpred12 = alpha[1] + as.matrix(X) %*% matrix(beta_s, nrow= (ncol(dataset)-2), ncol=1)
  linpred23 = alpha[2] + as.matrix(X) %*% matrix(beta_s, nrow= (ncol(dataset)-2), ncol=1)
  
  cusum1 = matrix(as.vector(inv_logit_vec(linpred12)), byrow = T, nrow=16, ncol=length(unique(dataset$SUBJECT_ID)))
  cusum2 = matrix(as.vector(inv_logit_vec(linpred23)), byrow = T, nrow=16, ncol=length(unique(dataset$SUBJECT_ID)))
  
  p1 = cusum1
  p2 = cusum2 - cusum1
  p3 = 1 - cusum2
  
  
  bool_Y1_pres = ( Us <= p1 )
  bool_Y2_pres = ( Us > p1 & Us <= (p1+p2) )
  bool_Y3_pres = ( Us > p1+p2 )
  
  sim_Y_pres = sim_Wp_pres
  
  #If bool_Wp_pres==T i.e. fluorosis present and bool_Y1_pres==T, assign to category 1
  sim_Y_pres = ifelse(!bool_Wp_pres & bool_Y1_pres, 1, sim_Y_pres) 
  
  #If bool_Wp_pres==T i.e. fluorosis present and bool_Y1_pres==T, assign to category 2
  sim_Y_pres = ifelse(!bool_Wp_pres & bool_Y2_pres, 2, sim_Y_pres) 
  
  #If bool_Wp_pres==T i.e. fluorosis present and bool_Y1_pres==T, assign to category 3
  sim_Y_pres = ifelse(!bool_Wp_pres & bool_Y3_pres, 3, sim_Y_pres) 
  
  #Finally, update the dataset with the generated responses  
  dataset[,"FRI_Score"] = as.vector(t(sim_Y_pres)) #unlist rowwise and populate the FRI_Score column in dataset
  
  
  table(dataset$FRI_Score)
  
  return(dataset)
}