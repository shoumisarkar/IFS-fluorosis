
#combined presence

age=9

setwd(paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/N10000/combined/independence,independence/age", age))

obj = load("coef_pres_JK_MC_2_chunk_1.Rdata")
assign("coef_pres_JK", get(obj))

#coef_pres_JK = data.frame(coef_pres_JK)

# Assume 'data' is your matrix
# Check which columns have any non-NA values
valid_columns <- colSums(!is.na(coef_pres_JK)) > 0

# Subset the matrix to keep only columns with data
dat_JK <- coef_pres_JK[, valid_columns]
dat_JK = coef_pres_JK[,1:3000]


for(chunk in 3:9)
{
  
  #load the chunk of 1000 leave-one-cluster-out estimates
  obj1 = load(paste0("coef_pres_JK_MC_2_chunk_", chunk ,".Rdata"))
  assign("obj1", get(obj1))
  
  #coef_pres_JK = data.frame(coef_pres_JK)
  
  # Assume 'data' is your matrix
  # Check which columns have any non-NA values
  valid_columns <- colSums(!is.na(obj1)) > 0
  
  # Subset the matrix to keep only columns with data
  obj1 <- obj1[, valid_columns]
  obj1 = obj1[,-c(1)] #removed the 1st column as it overlaps with the contents of coef_pres_JK
  
  print(chunk)
  print(ncol(obj1))
  
  dat_JK = cbind(dat_JK, obj1)
}


valid_columns <- colSums(!is.na(dat_JK)) > 0

# Subset the matrix to keep only columns with data
dat_JK <- dat_JK[, valid_columns]
dat_JK = dat_JK[,-c(1)] #removed the 1st column as it overlaps with the contents of coef_pres_JK


#coef_pres_MC_2
#JK_SD_MC_2
#std_coef_pres_MC2

JK_SD_MC_2 = apply(dat_JK, MARGIN = 1, FUN = function(x){sd(x, na.rm = T)})

#Need to get coef_pres_MC_2

#std_coef_pres_MC2 = 