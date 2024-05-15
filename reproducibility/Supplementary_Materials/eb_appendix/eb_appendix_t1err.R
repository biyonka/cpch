library(HDMT)
library(DACT)
library(adaFilter)
library(locfdr)
library(ggplot2)

set.seed(211)
# First read in the arguments listed at the command line
args=(commandArgs(TRUE))# args is now a list of character vectors
# First check to see if arguments are passed.
# Then cycle through each element of the list and evaluate the expressions.
arguments = c('max_null_mu', 'max_alt_mu', 'alpha', 'repl', 'power_iters')

if (length (args) == 0) {
  print("No arguments supplied.")
  ## supply default values
  max_null_mu = 4
  max_alt_mu = 4
  alpha = 0.1
  repl = 5000 #change this to vary M
  power_iters = 100
} else {
  for (i in 1:length(args)) {
    eval (parse (text = paste(arguments[i]," = ",args[[i]],sep="") ))
  }
}

print(c(max_null_mu, max_alt_mu, alpha, repl, power_iters))
alt_step = 0.05
null_step = 1
null_mus = lapply(seq(0, max_null_mu, null_step), function(x){return (c(0, x))})

##########################
### T1 Error ###
##########################


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


generate_XX = function(x, repl, mus){
  XX = matrix(rnorm(2*repl, mean = mus, 1), nrow = repl, byrow = TRUE)
  return (XX)
}


get_pvals = function(XX, alpha){
  XX_pvals = 2*(pnorm(-abs(XX)))
  
  dact_pvals = DACT(p_a=XX_pvals[,1], p_b=XX_pvals[, 2], correction = 'JC')
  cpch_pvals = apply(XX, 1, calc_cpch_pvals)
  ada_pvals = adaFilter(XX_pvals, r = 2, alpha = alpha)$adjusted.p
  #ada_decisions = adaFilter(XX_pvals, r = 2, alpha = alpha)$decision
  maxp_pvals = apply(XX_pvals, 1, max, na.rm = TRUE)
  
  dact_p_reject = mean(dact_pvals <= alpha, na.rm = TRUE)
  cpch_p_reject = mean(cpch_pvals <= alpha)
  ada_p_reject = mean(ada_pvals <= alpha)
  maxp_p_reject = mean(maxp_pvals <= alpha)
  #ada_decision = mean(ada_decisions)
  return (list(cpch_p_reject, ada_p_reject, dact_p_reject, maxp_p_reject))
}


get_p_reject = matrix(rep(0, length(null_mus)*4), ncol = 4)
for (i in seq(1, length(null_mus), 1)){
  XX_list =lapply(seq(1, power_iters, 1), generate_XX, repl, null_mus[[i]])
  dh = sapply(XX_list, get_pvals, alpha)
  pm = matrix(unlist(dh), nrow = power_iters,  byrow = TRUE)
  p_reject = colMeans(pm, na.rm = TRUE)
  get_p_reject[i,] = p_reject
}

t1_error = data.frame(get_p_reject)
t1_error['mus'] = seq(0, max_null_mu, null_step)
colnames(t1_error) = c( 'cPCH', 'adaFilter', 'DACT', 'Max-P', 'mus')
t1_error_l = reshape(t1_error, direction = 'long', varying = list(1:4), times = c('cPCH', 'adaFilter', 'DACT', 'Max-P'), v.names = 'p_reject')
t1_error_l = t1_error_l[, c(1, 2, 3)]
colnames(t1_error_l) = c( 'nonzero_mu', 'method', 't1_error')
t1_error_l$method <- factor(t1_error_l$method, levels = c('cPCH', 'adaFilter', 'DACT', 'Max-P'))

#change path to where you want data to be saved
write.csv(t1_error_l, file = sprintf("t1_error_mu_%1.0f_a_%.2f_M_%1.0f.csv", max_null_mu, alpha, repl),  row.names = FALSE)


