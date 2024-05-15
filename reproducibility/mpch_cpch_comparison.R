setwd("~/Documents/Research/cpch") #set directory to where this file is located locally
set.seed(21)
library(mediation.test)
library(ggplot2)

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


calc_mpch_pvals = function(t, true_mu, oracle = FALSE){
  t1 = t[1]
  t2 = t[2]
  t1_order = min(abs(t1),abs(t2))
  t2_order = max(abs(t1),abs(t2))
  
  if (oracle == TRUE){
    min_mu = min(abs(true_mu[1]),abs(true_mu[2]))
    max_mu = max(abs(true_mu[1]),abs(true_mu[2]))
    if (max_mu == abs(true_mu[2])){mu_star = true_mu[2]} else {mu_star = true_mu[1]}
  }else if (t2_order == abs(t2)){mu_star = t2} else {mu_star = t1}
  
  pval = (2*pnorm(t1_order, 0, 1, lower.tail = FALSE))*(pnorm(-t1_order, mu_star, 1, lower.tail = TRUE) + pnorm(t1_order, mu_star, 1, lower.tail = FALSE))
  
  return(pval)
}

##############
#Type I Error
##############
it = 100000
alpha = 0.05

calc_pvals = function(nonzero_mu) {
  XX = matrix(c(rnorm(it, 0, 1), rnorm(it, nonzero_mu, 1)), ncol = 2, byrow = FALSE)
  XX_pvals = 2*(1-pnorm(abs(XX)))
  maxp_decisions = apply(XX_pvals, 1, max) <= alpha 
  cpch_decisions =  apply(XX, 1, calc_cpch_pvals, true_mu = c(0, nonzero_mu), oracle = F) <= 0.0425 #using value from lookup table
  mpch_decisions =  apply(XX, 1, calc_mpch_pvals,  true_mu = c(0, nonzero_mu),  oracle = F) <= alpha
  mediation_test = mediation_test(XX, alpha)$decision
  return (c(sum(cpch_decisions)/it, 
            sum(mpch_decisions)/it, 
            sum(maxp_decisions)/it,
            sum(mediation_test)/it,
            sd(cpch_decisions)/sqrt(it), 
            sd(mpch_decisions)/sqrt(it),
            sd(maxp_decisions)/sqrt(it),
            sd(mediation_test)/sqrt(it)
  ))
}

mu_vec = seq(0, 5, 0.5)
pval_ses = sapply(mu_vec, calc_pvals)

pvals_df = data.frame('nonzero_mu' = rep(mu_vec, 4), 
                      'Method' = c(rep('cPCH', length(mu_vec)), rep('mPCH', length(mu_vec)), rep('Max-P', length(mu_vec)), rep('MM Optimal', length(mu_vec))),
                      't1_error' = c(pval_ses[1,], pval_ses[2,], pval_ses[3,], pval_ses[4,]),
                      'ses' = c(pval_ses[5,], pval_ses[6,], pval_ses[7, ], pval_ses[8, ]))



pvals_df$Method = factor(pvals_df$Method, c('Max-P', 'mPCH', 'cPCH', 'MM Optimal'))


t1err = ggplot(pvals_df, aes(x = nonzero_mu, y = t1_error, color = Method)) +
  geom_line() + ylim(c(0, 0.25)) +
  theme_minimal() +
  geom_errorbar(aes(ymin= t1_error-2*ses, ymax=t1_error+2*ses), width=.1, show.legend = FALSE)  +
  geom_hline(aes(yintercept = 0.05), color = 'black', linetype = 'dotted') + 
  #geom_segment(aes(x=0, xend=5, y=0.05, yend=0.05), color = 'black', linetype = 'dashed') + 
  xlab(expression(paste(theta[(2)]))) + ylab('Type I Error') + 
  theme(aspect.ratio = 1, 
        axis.title =element_text(size=12), axis.text = element_text(size = 12),
        legend.title=element_text(size=12), legend.text=element_text(size=12))




##############
#Power 
##############
it = 100000
alpha = 0.05

calc_pvals_power = function(nonzero_mu) {
  XX = matrix(c(rnorm(it, nonzero_mu/2, 1), rnorm(it, nonzero_mu, 1)), ncol = 2, byrow = FALSE)
  XX_pvals = 2*(1-pnorm(abs(XX)))
  maxp_decisions = apply(XX_pvals, 1, max) <= alpha 
  cpch_decisions =  apply(XX, 1, calc_cpch_pvals, true_mu = c(nonzero_mu/2, nonzero_mu), oracle = F) <= 0.043 
  mpch_decisions =  apply(XX, 1, calc_mpch_pvals,  true_mu = c(nonzero_mu/2, nonzero_mu),  oracle = F) <= alpha 
  mediation_test = minimax_test(XX, alpha)$decision
  return (c(sum(cpch_decisions)/it, 
            sum(mpch_decisions)/it, 
            sum(maxp_decisions)/it,
            sum(mediation_test)/it,
            sd(cpch_decisions)/sqrt(it), 
            sd(mpch_decisions)/sqrt(it),
            sd(maxp_decisions)/sqrt(it),
            sd(mediation_test)/sqrt(it)
  ))
}


mu_vec = seq(0, 5, 0.5)
power_pval_ses = sapply(mu_vec, calc_pvals_power)

pvals_power_df = data.frame('nonzero_mu' = rep(mu_vec, 4),
                            'Method' = c(rep('cPCH', length(mu_vec)), rep('mPCH', length(mu_vec)), rep('Max-P', length(mu_vec)),  rep('MM Optimal', length(mu_vec))),
                            'power' = c(power_pval_ses[1,], power_pval_ses[2,], power_pval_ses[3,],  power_pval_ses[4,]),
                            'ses' = c(power_pval_ses[5,], power_pval_ses[6,], power_pval_ses[7, ], power_pval_ses[8, ]))


pvals_power_df$Method = factor(pvals_power_df$Method, c('Max-P','mPCH', 'cPCH', 'MM Optimal'))

ggplot(pvals_power_dfl, aes(x = nonzero_mu, y = power, color = Method)) +
  geom_line(size = 1.25) + 
  theme_minimal() +
  geom_errorbar(aes(ymin= power-2*ses, ymax=power+2*ses), width=.1)  +
  geom_hline(aes(yintercept = 0.05), color = 'black', linetype = 'dotted') +
  xlab(expression(paste(theta[(2)]))) + ylab('Power') + 
  theme(aspect.ratio = 1, 
        axis.title =element_text(size=30), axis.text = element_text(size = 29),
        legend.title=element_text(size=29), legend.text=element_text(size=28))


