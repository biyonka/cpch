library(AdaPTGMM)
library(ggplot2)
library(adaFilter)
iteration=1
set.seed(iteration)
q = 0.1


library(dplyr)
library(tidyr)

storey = function(p, fdr.level){
  qvals_out <- p
  rm_na <- !is.na(p)
  p = p[rm_na]
  pi0hat = min(1, (1+length(p) - sum(p <= 0.5))/(length(p)/2))
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  pi0s = list()
  pi0s$pi0 = pi0hat
  qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m /i ))[ro]
  qvals_out[rm_na] <- qvals
  return (list(qvalues = qvals_out,significant = (qvals <= fdr.level)))
}

# 
new_dmd = read.csv('PATH/dmd.csv')
bhqvals_cpch = p.adjust(new_dmd$pvalue.2, method = 'BH', n = 1871)
rejected_bh_cpch = 1*(bhqvals_cpch <= q) 
sum(rejected_bh_cpch)
new_dmd['rej_first'] = rejected_bh_cpch

#run get_pch_pvals.ipynb to generate the below data
# #read in classical pvals data frame
cpch_pvals_dmd = read.csv('PATH/cpch_fisher_pvals_dmd_followup_screen.csv')
# #read in cpch pvals data frame
fisher_pvals_dmd = read.csv('PATH/fisher_pvals_dmd_followup_screen.csv')



pi1fx = function(x){
  if (x >= 0.95){pi1 = 1}else{pi1 = 0.1}
  return(pi1)
}
q = 0.1
alphas = c(q)
formulas = paste("splines::ns(x,df=", c(2,3), ")")

  x <- data.frame(x = new_dmd$rej_first) 
  #adapt cpch, r=3
  cpch_adapt_3_3 <- adapt_gmm(x=x, pvals=cpch_pvals_dmd$cpch_3_3, alphas=alphas,
                         beta_formulas = formulas, model_type = "neural", nclasses= c(5))
  rejected_adapt_cpch_3_3 = 1*(cpch_adapt_3_3$qvals <= q)
  sum(rejected_adapt_cpch_3_3)
  #adapt cpch, r=2
  cpch_adapt_3_2 <- adapt_gmm(x=x,pvals=cpch_pvals_dmd$cpch_3_2, alphas=alphas,
                        beta_formulas = formulas, model_type = "neural", nclasses= c(5))
  rejected_adapt_cpch_3_2 = 1*(cpch_adapt_3_2$qvals <= q)
  sum(rejected_adapt_cpch_3_2)
  
  #cpch_storey
  rejected_storey_cpch_3_2 = 1*(storey(cpch_pvals_dmd$cpch_3_2, q)$significant)
  sum(rejected_storey_cpch_3_2)
  rejected_storey_cpch_3_3 = 1*(storey(cpch_pvals_dmd$cpch_3_3, q)$significant)
  sum(rejected_storey_cpch_3_3)

  #cpch BH
  bhqvals_cpch_3_2 = p.adjust(cpch_pvals_dmd$cpch_3_2 , method = 'BH', n = nrow(new_dmd))
  rejected_bh_cpch_3_2 = 1*(bhqvals_cpch_3_2 <= q)
  sum(rejected_bh_cpch_3_2)
  bhqvals_cpch_3_3 = p.adjust(cpch_pvals_dmd$cpch_3_3 , method = 'BH', n = nrow(new_dmd))
  rejected_bh_cpch_3_3 = 1*(bhqvals_cpch_3_3 <= q)
  sum(rejected_bh_cpch_3_3)
  
  #adapt fisher classical, r=3
  fisher_adapt_3_3 <- adapt_gmm(x=x, pvals=fisher_pvals_dmd$fisher_pvals_3_3, alphas=alphas,
                              beta_formulas = formulas, model_type = "neural", nclasses= c(5))
  rejected_adapt_fisher_3_3 = 1*(fisher_adapt_3_3$qvals <= q)
  sum(rejected_adapt_fisher_3_3)
  #adapt fisher classical, r=2
  fisher_adapt_3_2 <- adapt_gmm(x=x,pvals=fisher_pvals_dmd$fisher_pvals_3_2, alphas=alphas,
                              beta_formulas = formulas, model_type = "neural", nclasses= c(5))
  rejected_adapt_fisher_3_2 = 1*(fisher_adapt_3_2$qvals <= q)
  sum(rejected_adapt_fisher_3_2)
  
  #fisher classical storey
  rejected_storey_fisher_3_2 = 1*(storey(fisher_pvals_dmd$fisher_pvals_3_2, q)$significant)
  sum(rejected_storey_fisher_3_2)
  rejected_storey_fisher_3_3 = 1*(storey(fisher_pvals_dmd$fisher_pvals_3_3, q)$significant)
  sum(rejected_storey_fisher_3_3)
  
  
  #fisher classical BH
  bhqvals_fisher_3_2 = p.adjust(fisher_pvals_dmd$fisher_pvals_3_2 , method = 'BH', n = nrow(new_dmd))
  rejected_bh_fisher_3_2 = 1*(bhqvals_fisher_3_2 <= q)
  sum(rejected_bh_fisher_3_2)
  bhqvals_fisher_3_3 = p.adjust(fisher_pvals_dmd$fisher_pvals_3_3, method = 'BH', n = nrow(new_dmd))
  rejected_bh_fisher_3_3 = 1*(bhqvals_fisher_3_3 <= q)
  sum(rejected_bh_fisher_3_3)
  
  #adafilter
  rejections_ada_3_2 = adaFilter(as.matrix(dmd[,c('pvalue.1', 'pvalue.4', 'pvalue.3')]), r=2, alpha = q)$decision
  sum(rejections_ada_3_2)
  rejections_ada_3_3 = adaFilter(as.matrix(dmd[,c('pvalue.1', 'pvalue.4', 'pvalue.3')]), r=3, alpha = q)$decision
  sum(rejections_ada_3_3)
  

  dhh_pvals_3_2 = (10)*fisher_pvals_dmd$fisher_pvals_3_2[fisher_pvals_dmd$fisher_pvals_3_2 <= 0.1]
  dhh_x_3_2 = x[fisher_pvals_dmd$fisher_pvals_3_2 <= 0.1,]
  dhh_pvals_3_3 = (10)*fisher_pvals_dmd$fisher_pvals_3_3[fisher_pvals_dmd$fisher_pvals_3_3 <= 0.1]
  dhh_x_3_3 = x[fisher_pvals_dmd$fisher_pvals_3_3 <= 0.1,]
  #adapt dhh
  if (length(dhh_pvals_3_2) >= 1){
    dhh_x_3_2 = data.frame(x = dhh_x_3_2)
    ada_dhh_3_2 <- adapt_gmm(x=dhh_x_3_2, pvals=dhh_pvals_3_2, alphas=alphas,
                         beta_formulas = formulas, model_type = "neural", nclasses= c(5))
    rejected_adapt_dhh_3_2 = 1*(ada_dhh_3_2$qvals <= q)
  } else{
    rejected_adapt_dhh_3_2 = 0
  }
  if (length(dhh_pvals_3_3) >= 1){
    dhh_x_3_3 = data.frame(x = dhh_x_3_3)
    ada_dhh_3_3 <- adapt_gmm(x=dhh_x_3_3, pvals=dhh_pvals_3_3, alphas=alphas,
                             beta_formulas = formulas, model_type = "neural", nclasses= c(5))
    rejected_adapt_dhh_3_3 = 1*(ada_dhh_3_3$qvals <= q)
  } else{
    rejected_adapt_dhh_3_3 = 0
  }
  

 print('r & cPCH-AdaPTGMM & DHH-AdaPT-GMM & adaFilter')
 cat(paste0('2 & ',sum(rejected_adapt_cpch_3_2), ' & ', sum(rejected_adapt_fisher_3_2),' & ',sum(rejected_adapt_dhh_3_2),' & ',sum(rejections_ada_3_2),'\\
3 & ',sum(rejected_adapt_cpch_3_3),' & ',sum(rejected_adapt_fisher_3_3),' & ',sum(rejected_adapt_dhh_3_3),' & ',sum(rejections_ada_3_3),'\\'))
 
 
