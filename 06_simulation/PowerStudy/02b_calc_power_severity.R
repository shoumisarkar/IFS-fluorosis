
empirical_power = function(target_alpha = 0.05, M=1000, beta=-2)
{
  
  setwd(paste0("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy/N30 - beta=",
               beta, "/severity/independence,independence/age23"))
  
  std_estimates = c()
  signif_vec = c()
  
  # count = 0
  
  for(mc_seed in 1:M) 
  {
    filename = paste0("bs_Avghomeppm_MC_", mc_seed,".RData")
    
    alpha = target_alpha
    
    if(!file.exists(filename))
    {
      next
    }
    
    # count = count + 1
    
    obj = load(filename)
    assign("bs_Avghomeppm", get(obj))
    
    lci = quantile(bs_Avghomeppm, probs = alpha/2, na.rm = T)
    uci = quantile(bs_Avghomeppm, probs = (1 - alpha/2), na.rm = T)
    
    mean = mean(bs_Avghomeppm, na.rm = T)
    SD = sd(bs_Avghomeppm, na.rm = T)
    std_estimate = mean/SD
    
    is_signif = T
    
    if(lci<=0 && 0<=uci)
    {
      is_signif = F
    }
    
    signif_vec = c(signif_vec, is_signif)
    std_estimates = c(std_estimates, std_estimate)
    
    # if(count == 100)
    # {
    #   break
    # }
  }
  
  power = length(which(signif_vec))/length((signif_vec))
  effect_size = mean(std_estimates, na.rm=T)
  
  return(list(power=power, effect_size=effect_size))
}

#empirical_power(target_alpha = 0.05)

#Get power curve:

power = c()
effect_size = c()

for(beta in c(-2, -1, 0, 1, 2))
{
  obj = empirical_power(beta=beta)
  
  power = c(power, obj$power) 
  effect_size = c(effect_size, obj$effect_size)
}

plot(x = effect_size[order(effect_size)], y = power[order(effect_size)], type = "b",
     main = "Separate severity model")


data.frame(beta = betas,
           effect_size = effect_size,
           power = power)



# plot(x = effect_size[order(effect_size)], y = power[order(effect_size)], type = "l")
# c(-2, -1, 0, 1, 2)[order(effect_size)]


# # setwd("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy/N30/")
# # obj = load("valid_MC_indices_combined.RData")
# # assign("comb_ind", get(obj))
# 
# 
# M = 300
# setwd("W:/somnath.datta/shoumisarkar/Fluorosis/Results/06_simulation/PowerStudy/N30 - zero beta - level/severity/independence,independence/age23")
# 
# # sev_ind = c()
# # signif_vec_0 = c()
# # 
# # for(mc_seed in 1:M)
# # {
# # 
# #   filename = paste0("bs_Avghomeppm_MC_", mc_seed,".RData")
# #   
# #   if(!file.exists(filename))
# #   {
# #     next
# #   }
# #   
# #   sev_ind = c(sev_ind, mc_seed)
# #   
# #   obj = load(filename)
# #   assign("bs_Avghomeppm", get(obj))
# #   
# #   alpha=0.05
# #   
# #   lci = quantile(bs_Avghomeppm, probs = alpha/2, na.rm = T)
# #   uci = quantile(bs_Avghomeppm, probs = (1 - alpha/2), na.rm = T)
# #   
# #   is_signif = T
# #   
# #   if(lci<=0 && 0<=uci)
# #   {
# #     is_signif = F
# #   }
# #   
# #   signif_vec_0 = c(signif_vec_0, is_signif)
# #   
# # }
# # 
# # library(dplyr)
# # true_ind_sev = sev_ind[which(signif_vec_0)][!(sev_ind[which(signif_vec_0)] %in% comb_ind)][1:3]
# # 
# # use_these_ind = sort(c(comb_ind, true_ind_sev ))
# 
# 
# 
# empirical_power = function(target_alpha = 0.05)
# {
#   
#   #print(target_alpha) 
#   
#   signif_vec = c()
#   
#   for(mc_seed in 1:M) #use_these_ind)#comb_ind)#all_indices)#res$indices)#1:M)
#   {
#     filename = paste0("bs_Avghomeppm_MC_", mc_seed,".RData")
#     
#     alpha = target_alpha
#     #alpha = target_alpha/2  ##
#     
#     if(!file.exists(filename))
#     {
#       next
#     }
#     
#     obj = load(filename)
#     assign("bs_Avghomeppm", get(obj))
#     
#     #bs_Avghomeppm = bs_Avghomeppm[1:1000]
#     
#     # lci = quantile(bs_Avghomeppm, probs = 0.025, na.rm = T)
#     # uci = quantile(bs_Avghomeppm, probs = 0.975, na.rm = T)
#     
#     lci = quantile(bs_Avghomeppm, probs = alpha/2, na.rm = T)
#     uci = quantile(bs_Avghomeppm, probs = (1 - alpha/2), na.rm = T)
#     
#     is_signif = T
#     
#     if(lci<=0 && 0<=uci)
#     {
#       is_signif = F
#     }
#     
#     signif_vec = c(signif_vec, is_signif)
#     
#     #print(mc_seed)
#   }
#   
#   #print(length(which(signif_vec)))
#   #print(length((signif_vec)))
#   #print(length(which(signif_vec))/length((signif_vec)))
#   
#   return(length(which(signif_vec))/length((signif_vec)))
# }
# 
# empirical_power(target_alpha = 0.05)
# 
