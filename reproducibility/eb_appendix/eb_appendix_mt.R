library(ggplot2)
library(DACT)
library(HDMT)
library(cowplot)
library(adaFilter)

set.seed(21)
# First read in the arguments listed at the command line
args=(commandArgs(TRUE))# args is now a list of character vectors
# First check to see if arguments are passed.
# Then cycle through each element of the list and evaluate the expressions.
arguments = c('mag', 'alpha', 'M', 'power_iters')

if (length (args) == 0) {
  print("No arguments supplied.")
  ## supply default values
  mag = 2
  alpha = 0.1
  M = 5000 #change this to get results for different M
  power_iters = 1000
} else {
  for (i in 1:length(args)) {
    eval (parse (text = paste(arguments[i]," = ",args[[i]],sep="") ))
  }
}

dense_null = c(0.6, 0.2, 0.2, 0)
sparse_null = c(0.9, 0.05, 0.05, 0)
complete_null =c(1, 0, 0, 0)
sparse_alt = c(0.88, 0.05, 0.05, 0.02)
dense_alt =  c(0.4, 0.2, 0.2, 0.2)

null_cases = list(dense_null, sparse_null, complete_null)
alt_cases = list(dense_alt, sparse_alt)


#calculates the cPCH p-value for r=m=2
calc_cpch_pvals = function(x, true_mu, oracle = FALSE){
  x1 = x[1]
  x2 = x[2]
  min_val = min(abs(x1),abs(x2))
  max_val = max(abs(x1),abs(x2))
  
  if (oracle == TRUE){
    min_mu = min(abs(true_mu[1]),abs(true_mu[2]))
    max_mu = max(abs(true_mu[1]),abs(true_mu[2]))
    if (max_mu == abs(true_mu[2])){est = c(0, true_mu[2])} else {est = c(true_mu[1], 0)}
  }else if (max_val == abs(x2)){
    est = c(0, x2)
  } else {est = c(x1, 0)}
  
  if (max_val == abs(x2)){
    t2 = x2
  } else {   
    t2 = x1
  }
  
  
  if (est[1] >= 0){
    mix_comp_1 = (-pnorm(min_val,est[1],1,T)+pnorm(max_val,est[1],1,T) +
                    pnorm(-min_val,est[1],1,T)-pnorm(-max_val,est[1],1,T))/(pnorm(max_val,est[1],1, T)-pnorm(-max_val,est[1],1, T))
    
    mw_1 = (dnorm(t2, est[2], 1))*(pnorm(max_val-est[1]) - pnorm(-max_val-est[1]))
    
  }
  else if (est[1] < 0){
    mix_comp_1 =  (pnorm(min_val,est[1],1,F)-pnorm(max_val,est[1],1,F) -
                     pnorm(-min_val,est[1],1,F)+pnorm(-max_val,est[1],1,F))/(-pnorm(max_val,est[1],1,F)+pnorm(-max_val,est[1],1, F))
    mw_1 = (dnorm(t2,est[2], 1))*(-pnorm(max_val-est[1], 0, 1, F) + pnorm(-max_val-est[1], 0, 1, F))
    
  }
  
  if (est[2] >= 0){
    mix_comp_2 = (-pnorm(min_val,est[2],1,T)+pnorm(max_val,est[2],1,T) +
                    pnorm(-min_val,est[2],1,T)-pnorm(-max_val,est[2],1,T))/(pnorm(max_val,est[2],1, T)-pnorm(-max_val,est[2],1, T))
    mw_2 =  (dnorm(t2,est[1], 1))*(pnorm(max_val-est[2], 0, 1, T) - pnorm(-max_val-est[2], 0, 1, T))
    
  }
  else if (est[2] < 0){
    mix_comp_2 =  (pnorm(min_val,est[2],1,F)-pnorm(max_val,est[2],1,F) -
                     pnorm(-min_val,est[2],1,F)+pnorm(-max_val,est[2],1,F))/(-pnorm(max_val,est[2],1,F)+pnorm(-max_val,est[2],1, F))
    mw_2 = (dnorm(t2,est[1], 1))*(-pnorm(max_val-est[2], 0, 1, F) + pnorm(-max_val-est[2], 0, 1, F))
    
  }
  
  Ps = mw_2 + mw_1
  
  pval = mix_comp_2 * (mw_2/Ps) +mix_comp_1 * (mw_1/Ps)
  
  return(pval)
}



generate_XX = function(x, case_config, M, mus){
  num_cases = (case_config * M)
  mus_mat =  matrix(c(c(0, 0), c(mus[1], 0), c(0, mus[2]), mus), nrow = 4, byrow = TRUE)
  
  generate_data_reps = function(i){
    if (num_cases[i] > 0){
      return (matrix(rnorm(2*num_cases[i], mean = mus_mat[i,], 1), nrow = num_cases[i], byrow = TRUE))
    }
  }
  #use lapply to generate num_cases[i] number of draws for each case
  XX = lapply(seq(1, 4, 1), generate_data_reps)
  #rbind together all results
  XX = do.call(rbind, XX)
  null_bool = rep(c(T, T, T, F), num_cases)
  XX = cbind(null_bool, XX)
  return (XX)
}


get_rejections = function(XX, alpha){
  null_bool = as.logical(XX[,1])
  XX_pvals = 2*(pnorm(-abs(as.matrix(XX[, c(2, 3)]))))
  
  dact_pvals = DACT(p_a=XX_pvals[,1], p_b=XX_pvals[, 2], correction = "JC")
  dact_bh = p.adjust(dact_pvals, method ='BH')
  
  ada_decisions = adaFilter(XX_pvals, r = 2, alpha = alpha)$decision
  
  alphas= null_estimation(XX_pvals, lambda = 0.5)
  alpha00 = alphas$alpha00
  alpha01 = alphas$alpha01
  alpha10 = alphas$alpha10
  alpha1 =  alphas$alpha1
  alpha2 =  alphas$alpha2
  hdmt_loc = fdr_est(alpha00, alpha01, alpha10, alpha1, alpha2, XX_pvals, exact = 1)
  
  cpch_pvals = apply(XX[,c(2, 3)], 1, calc_cpch_pvals)
  cpch_bh = p.adjust(cpch_pvals, method = 'BH')
  
  #maxp with BH
  maxp_decisions = ClassicalMTP(XX_pvals, r = 2, alpha = alpha, method = "Bonferroni")$decision
  
  ##### FDR and Power Calculations #####
  cpch_bh_fdr = sum(cpch_bh[null_bool] <= alpha, na.rm = TRUE)/max(1, sum(cpch_bh <= alpha, na.rm = TRUE))
  
  dact_fdr = sum(dact_bh[null_bool] <= alpha, na.rm = TRUE)/max(1, sum(dact_bh <= alpha, na.rm = TRUE))
  
  hdmt_fdr = sum(hdmt_loc[null_bool] <= alpha)/max(1, sum(hdmt_loc <= alpha))
  
  ada_fdr = sum(ada_decisions[null_bool])/max(1, sum(ada_decisions))
  
  maxp_fdr = sum(maxp_decisions[null_bool])/max(1, sum(maxp_decisions))
  
  if (sum(null_bool) < length(null_bool)){
    cpch_bh_power = mean(cpch_bh[!null_bool] <= alpha)
    dact_power = mean(dact_bh[!null_bool] <= alpha)
    hdmt_power = mean(hdmt_loc[!null_bool] <= alpha)
    ada_power = mean(ada_decisions[!null_bool])
    maxp_power = mean(maxp_decisions[!null_bool])
    return (list(cpch_bh_fdr,
                 ada_fdr,
                 dact_fdr,
                 hdmt_fdr,
                 maxp_fdr, cpch_bh_power, 
                 ada_power, dact_power, hdmt_power, maxp_power
    ))
  }
  else{return (list(cpch_bh_fdr, 
                    ada_fdr, dact_fdr, hdmt_fdr, maxp_fdr))
  }
}


######################
##### Null Cases #####
######################

calc_null_cases = function(case, M, mus, alpha, power_iters){
  print(case)
  XX_list =lapply(seq(1, power_iters, 1), generate_XX, case, M, mus)
  
  dh =sapply(XX_list, get_rejections, alpha)
  power_mat = t(as.matrix(dh))
  pm = matrix(unlist(power_mat), nrow = power_iters,  byrow = FALSE)
  fdr = colMeans(pm, na.rm = TRUE)
  fdr_ses = apply(pm, 2, sd, na.rm = FALSE)/(power_iters ** 0.5) 
  print(cat(
    paste('cPCH-BH & ',
          toString(round(fdr[1], digits = 6)),
          ' \\\\ \n & adaFilter  &',
          toString(round(fdr[2], digits = 6)),
          ' \\\\ \n & DACT & ',
          toString(round(fdr[3], digits = 6)),
          ' \\\\ \n & HDMT & ',
          toString(round(fdr[4], digits = 6)),
          ' \\\\ \n & MaxP-BH & ',
          toString(round(fdr[5], digits = 6)),
          "\\\\")
  )
  )
  return (list(fdr, fdr_ses))
}

#####################
##### Alt Cases #####
#####################
calc_cases = function(case, M, mus, alpha, power_iters){
  print(case)
  XX_list =lapply(seq(1, power_iters, 1), generate_XX, case, M, mus)
  dh =sapply(XX_list, get_rejections, alpha)
  power_mat = t(as.matrix(dh))
  pm = matrix(unlist(power_mat), nrow = power_iters,  byrow = FALSE)
  t1_power = colMeans(pm, na.rm = TRUE)
  t1_power_ses = apply(pm, 2, sd, na.rm = FALSE)/(power_iters ** 0.5) 
  print(cat(paste('cPCH-BH & ',
                  toString(round(t1_power[1], digits = 4)), ' & ',
                  toString(round(t1_power[6], digits = 4)),
                  '\\\\ \n & adaFilter  & ',
                  toString(round(t1_power[2], digits = 4)), ' & ',
                  toString(round(t1_power[7], digits = 4)),
                  '\\\\ \n & DACT  & ',
                  toString(round(t1_power[3], digits = 4)), ' & ',
                  toString(round(t1_power[8], digits = 4)),
                  '\\\\ \n & HDMT  & ',
                  toString(round(t1_power[4], digits = 4)), ' & ',
                  toString(round(t1_power[9], digits = 4)),
                  '\\\\ \n & MaxP-BH  & ',
                  toString(round(t1_power[5], digits = 4)), ' & ',
                  toString(round(t1_power[10], digits = 4)),
                  '\\\\')))
  return (list(t1_power, t1_power_ses))
}

print(M)
alpha = 0.1
alt_res = lapply(alt_cases, calc_cases, M, c(mag, mag), alpha, power_iters)
null_res = lapply(null_cases, calc_null_cases, M, c(mag, mag), alpha, power_iters)

#set file directory to where you want the file to be saved
save(alt_res, null_res, file = paste0("eb_appendix_",M,".rda"))


