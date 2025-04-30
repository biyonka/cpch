setwd("~/Documents/Research/cpch") #set directory to where this file is located locally
set.seed(100)
#library(mediation.test)
library(ggplot2)



calc_cpch_pvals_bivariate = function(x, rho, true_mu, oracle = FALSE){
  x1 = x[1]
  x2 = x[2]
  f_obs = min(abs(x1),abs(x2))
  abs_t_2 = max(abs(x1),abs(x2))
  
  if (oracle == TRUE){
    min_mu = min(abs(true_mu[1]),abs(true_mu[2]))
    max_mu = max(abs(true_mu[1]),abs(true_mu[2]))
    if (max_mu == abs(true_mu[2])){est = c(0, true_mu[2])} else {est = c(true_mu[1], 0)}
  }else if (x1^2 <= x2^2){
    est = c(0, x2-rho*x1)
  } else {est = c(x1-rho*x2, 0)}
  
  if (abs_t_2 == abs(x2)){
    t_2 = x2
  } else {   
    t_2 = x1
  }
  
  
 denom_b1 = dnorm(t_2-est[2]) *(  pnorm(abs_t_2, mean =est[1]+rho*(t_2-est[2]), sd = (1-rho**2)**0.5) -
                                    pnorm(-abs_t_2, mean =est[1]+rho*(t_2-est[2]), sd = (1-rho**2)**0.5) )
 
  
  denom_b2 = dnorm(t_2-est[1]) * (pnorm((abs_t_2 - (est[2]+rho*(t_2-est[1])))/(1-rho**2)**0.5) -
                                    pnorm((-abs_t_2 - (est[2]+rho*(t_2-est[1])))/(1-rho**2)**0.5) )
  
  
  num_2 =  dnorm(t_2-est[1]) * (pnorm((abs_t_2 - (est[2]+rho*(t_2-est[1])))/(1-rho**2)**0.5) - pnorm((f_obs - (est[2]+rho*(t_2-est[1])))/(1-rho**2)**0.5) +
                                  pnorm((-f_obs - (est[2]+rho*(t_2-est[1])))/(1-rho**2)**0.5) - pnorm((-abs_t_2 - (est[2]+rho*(t_2-est[1])))/(1-rho**2)**0.5) 
                                )
  
 num_1 =  dnorm(t_2-est[2]) *(pnorm((abs_t_2 - (est[1]+rho*(t_2-est[2])))/(1-rho**2)**0.5) -pnorm((f_obs - (est[1]+rho*(t_2-est[2])))/(1-rho**2)**0.5)+
                                                   pnorm((-f_obs - (est[1]+rho*(t_2-est[2])))/(1-rho**2)**0.5)-  pnorm((-abs_t_2 - (est[1]+rho*(t_2-est[2])))/(1-rho**2)**0.5)
                                                 )
   
  pval = (num_1+num_2)/(denom_b1 + denom_b2) #+ num_2/(denom_b1 + denom_b2)

  return(pval)
}


calc_cpch_pvals_bivariate(c(4.99, 5), 0.90)

#calc_cpch_pvals(x)


it=400000
alpha = 0.05


#type I error is tiny when rho is large, and that makes sense! If I observe something sigificantly bigger than 0
#and I know that the two RVs are very closely tied in value, then I basically have a lot of evidence against the null
calc_pvals = function(nonzero_mu, rho=0.99) {
  x1 = rnorm(it, 0, 1)
  x2 = rho*x1 + ((1-rho**2)^0.5)*rnorm(it, 0, 1) + nonzero_mu
  XX = matrix(c(x1, x2), ncol = 2, byrow = FALSE)
 # XX_pvals = 2*(1-pnorm(abs(XX)))
  #maxp_decisions = apply(XX_pvals, 1, max) <= alpha 
  cpch_decisions =  apply(XX, 1, calc_cpch_pvals_bivariate, rho=rho, true_mu = c(0, nonzero_mu), oracle = F) <= alpha #using value from lookup table
 # mpch_decisions =  apply(XX, 1, calc_mpch_pvals,  true_mu = c(0, nonzero_mu),  oracle = F) <= alpha
 # mediation_test = mediation_test(XX, alpha)$decision
  return (c(sum(cpch_decisions)/it, 
          #  sum(mpch_decisions)/it, 
          #  sum(maxp_decisions)/it,
          #  sum(mediation_test)/it,
            sd(cpch_decisions)/sqrt(it)
           # sd(mpch_decisions)/sqrt(it),
            #sd(maxp_decisions)/sqrt(it),
            #sd(mediation_test)/sqrt(it)
  ))
}



pvals_over_rho = lapply(c(0.01, 0.1, 0.3, 0.5), function(rho){
  mu_vec = seq(0, 8, 0.5)
  pval_ses = sapply(mu_vec, calc_pvals, rho)
  pvals_df = data.frame('nonzero_mu' = rep(mu_vec, 1), 'rho' = rho,
                        'Method' = c(rep('cPCH', length(mu_vec))),#, rep('mPCH', length(mu_vec)), rep('Max-P', length(mu_vec)), rep('MM Optimal', length(mu_vec))),
                        't1_error' = c(pval_ses[1,]),#, pval_ses[2,], pval_ses[3,], pval_ses[4,]),
                        'ses' = c(pval_ses[2,])#, pval_ses[6,], pval_ses[7, ], pval_ses[8, ])
  )
  
  return(pvals_df)
})


pvals_df = bind_rows(pvals_over_rho, .id = "column_label")



custom_labeller <- labeller(
  .multi_line = FALSE,
  rho = function(x) {
    # Convert rho to the symbol expression
    paste0("ρ=", x)
  }
)


#pvals_df$Method = factor(pvals_df$Method, c('Max-P', 'mPCH', 'cPCH', 'MM Optimal'))


pvals_df = read.csv('~/Documents/Research/cpch/Data/bivariate_cpch_t1_err_it_400000_no_adjustment.csv')
t1err = ggplot(pvals_df, aes(x = nonzero_mu, y = t1_error)) +
  geom_line() + ylim(c(0, 0.2)) +
  theme_minimal() + facet_wrap(~rho,  nrow = 1, dir = 'h', labeller = custom_labeller)+
  geom_errorbar(aes(ymin= t1_error-2*ses, ymax=t1_error+2*ses), width=.1, show.legend = FALSE)  +
  geom_hline(aes(yintercept = 0.05), color = 'black', linetype = 'dotted') + 
  #geom_segment(aes(x=0, xend=5, y=0.05, yend=0.05), color = 'black', linetype = 'dashed') + 
  xlab(expression(paste(theta[(2)]))) + ylab('Type I error') + 
  theme(aspect.ratio = 1, 
        strip.text.x = element_text(size = 20), legend.position="none",
        axis.title =element_text(size=25), axis.text = element_text(size = 20),
        # legend.title=element_text(size=15), legend.text=element_text(size=12)
  ) 



t1err


#write.csv(pvals_df, '~/Documents/Research/cpch/Data/bivariate_cpch_t1_err_it_400000.csv')



#change directory to where you want plot to be saved
ggsave(filename = paste0("bivariate_cpch_t1err.eps"),
       plot = t1err, path = '~/Documents/Research/cpch/Plots/', bg = 'white',
       height = 5, width = 12)



##############
#Power 
##############
it = 10000
alpha = 0.05

#and I know that the two RVs are very closely tied in value, then I basically have a lot of evidence against the null
calc_pvals_power = function(nonzero_mu, rho=0.99) {
  x1 = rnorm(it, 0, 1) + nonzero_mu/2
  x2 = rho*(x1-nonzero_mu/2) + ((1-rho**2)^0.5)*rnorm(it, 0, 1) + nonzero_mu
  XX = matrix(c(x1, x2), ncol = 2, byrow = FALSE)
  # XX_pvals = 2*(1-pnorm(abs(XX)))
  #maxp_decisions = apply(XX_pvals, 1, max) <= alpha 
  cpch_decisions =  apply(XX, 1, calc_cpch_pvals_bivariate, rho=rho, true_mu = c(0, nonzero_mu), oracle = F) <= alpha #using value from lookup table
  # mpch_decisions =  apply(XX, 1, calc_mpch_pvals,  true_mu = c(0, nonzero_mu),  oracle = F) <= alpha
  # mediation_test = mediation_test(XX, alpha)$decision
  return (c(sum(cpch_decisions)/it, 
            #  sum(mpch_decisions)/it, 
            #  sum(maxp_decisions)/it,
            #  sum(mediation_test)/it,
            sd(cpch_decisions)/sqrt(it)
            # sd(mpch_decisions)/sqrt(it),
            #sd(maxp_decisions)/sqrt(it),
            #sd(mediation_test)/sqrt(it)
  ))
}



power_pvals_over_rho = lapply(c(0.01, 0.1, 0.3, 0.5), function(rho){
  mu_vec = seq(0, 8, 0.5)
  power_pval_ses = sapply(mu_vec, calc_pvals_power, rho)
  pvals_df = data.frame('nonzero_mu' = rep(mu_vec, 1), 'rho' = rho,
                        'Method' = c(rep('cPCH', length(mu_vec))),#, rep('mPCH', length(mu_vec)), rep('Max-P', length(mu_vec)), rep('MM Optimal', length(mu_vec))),
                        'power' = c(power_pval_ses[1,]),#, pval_ses[2,], pval_ses[3,], pval_ses[4,]),
                        'ses' = c(power_pval_ses[2,])#, pval_ses[6,], pval_ses[7, ], pval_ses[8, ])
  )
  
  return(pvals_df)
})


power_pvals_df = bind_rows(power_pvals_over_rho, .id = "column_label")




# pvals_power_df = data.frame('nonzero_mu' = rep(mu_vec, 4),
#                             'Method' = c(rep('cPCH', length(mu_vec)), rep('mPCH', length(mu_vec)), rep('Max-P', length(mu_vec)),  rep('MM Optimal', length(mu_vec))),
#                             'power' = c(power_pval_ses[1,], power_pval_ses[2,], power_pval_ses[3,],  power_pval_ses[4,]),
#                             'ses' = c(power_pval_ses[5,], power_pval_ses[6,], power_pval_ses[7, ], power_pval_ses[8, ]))
# 

#pvals_power_df$Method = factor(pvals_power_df$Method, c('Max-P','mPCH', 'cPCH', 'MM Optimal'))

power_pvals_df = read.csv('~/Documents/Research/cpch/Data/bivariate_cpch_power_it_10000.csv')
power_over_rho = ggplot(power_pvals_df, aes(x = nonzero_mu, y = power)) +
  geom_line(size = 0.75) + 
  theme_minimal() + facet_wrap(~rho,  nrow = 1, dir = 'h', labeller = custom_labeller)+
  geom_errorbar(aes(ymin= power-2*ses, ymax=power+2*ses), width=.1)  +
  geom_hline(aes(yintercept = 0.05), color = 'black', linetype = 'dotted') + 
  #geom_segment(aes(x=0, xend=5, y=0.05, yend=0.05), color = 'black', linetype = 'dashed') + 
  xlab(expression(paste(theta[(2)]))) + ylab('Power') + 
  theme(aspect.ratio = 1, 
        strip.text.x = element_text(size = 20), legend.position="none",
        axis.title =element_text(size=25), axis.text = element_text(size = 20),
        # legend.title=element_text(size=15), legend.text=element_text(size=12)
  )


power_over_rho


#write.csv(power_pvals_df, '~/Documents/Research/cpch/Data/bivariate_cpch_power_it_10000.csv')

#change directory to where you want plot to be saved
ggsave(filename = paste0("bivariate_cpch_pow.eps"),
       plot = power_over_rho, path = '~/Documents/Research/cpch/Plots/', bg = 'white',
       height = 5, width = 12)
