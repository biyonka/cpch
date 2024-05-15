library(AdaPTGMM)
library(adaFilter)
library(HDMT)
library(DACT)
library(mediation.test) #method of Miles & Chambaz (2021)

library(ggplot2)
# First read in the arguments listed at the command line
args=(commandArgs(TRUE))# args is now a list of character vectors
# First check to see if arguments are passed.
# Then cycle through each element of the list and evaluate the expressions.
arguments = c('iteration')

if (length (args) == 0) {
  print("No arguments supplied.")
  ## supply default values
  iteration = 1
} else {
  for (i in 1:length(args)) {
    eval (parse (text = paste(arguments[i]," = ",args[[i]],sep="") ))
  }
}

set.seed(iteration)
q = 0.1



lookup_table =read.csv('PATH-TO-CPCH-DIRECTORY/cpch/lookup_table.csv')
print('read in lookup table')
set.seed(iteration)
q = 0.1



#######
#calculate cpch p-value
#######
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


#only for m=r=2, so don't need to specify method
get_adjusted_pval = function(pval, M=2, r=2){
  #read in from lookup table
  #lookup_table[lookup_table$methods == 'Fisher' & lookup_table$m== m & lookup_table$r == r,]
  #get alpha prime that pval is closest to
  #get alpha for that alpha prime
  #if pval closest to 0, interpolate between (0, 0) and (alpha', alpha) for smallest alpha'
  #if pval closest to biggest alpha', interpolate using biggest alpha', alpha, and next biggest alpha', alpha
  #if pval somewhere in middle, interpolate between points
  lookup_subset =  lookup_table[lookup_table$methods == 'Fisher' & lookup_table$m== M & lookup_table$r == r,]
  closest_alpha_prime = lookup_subset$alpha_prime[which.min(abs(pval - lookup_subset$alpha_prime))]
  if(which.min(abs(pval - lookup_subset$alpha_prime)) == 1){
    x1 = 0
    y1 = 0
    x2 = lookup_subset$alpha_prime[1]
    y2 = lookup_subset$alpha[1]
    adjusted_pval =  (y2-y1)/(x2-x1) * (pval - x1) + y1
  }else if(which.min(abs(pval - lookup_subset$alpha_prime)) == length(lookup_subset$alpha_prime)){
    x1 = lookup_subset$alpha_prime[length(lookup_subset$alpha_prime)-1]
    y1 = lookup_subset$alpha[length(lookup_subset$alpha_prime)-1]
    x2 = 1
    y2 = 1
    adjusted_pval =  (y2-y1)/(x2-x1) * (pval - x1) + y1
  } else{
    x1 = closest_alpha_prime
    y1 = lookup_subset$alpha[which.min(abs(pval - lookup_subset$alpha_prime))]
    x2 = lookup_subset$alpha_prime[which.min(abs(pval - lookup_subset$alpha_prime))+1]
    y2 = lookup_subset$alpha[which.min(abs(pval - lookup_subset$alpha_prime))+1]
    adjusted_pval =  (y2-y1)/(x2-x1)* (pval - x1) + y1
  }
  return(adjusted_pval)
}



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

pi1fx = function(x){
  if (x >= 0.95){pi1 = 1}else{pi1 = 0.1}
  return(pi1)
}


run_sim = function(theta){
  # Generate data as in the Multiple PCH Testing with Covariates section of the paper
  m <- 10000
  q = 0.1
  x <- runif(m, min = 0, max = 1)
  pi1 = sapply(x, pi1fx)
  H1 <- 1 * (runif(m) <= pi1)
  H2 <- 1 * (runif(m) <= pi1)
  theta1 <- H1 * theta
  theta2 <- H2 * theta
  z1 <- rnorm(m, mean=theta1)
  z2 <- rnorm(m, mean=theta2)
  # Two sided testing
  pvals1 <- 2 * pnorm(abs(z1), lower.tail = FALSE)
  pvals2 <- 2 * pnorm(abs(z2), lower.tail = FALSE)
  maxp = pmax(pvals1, pvals2)
  z = data.frame(z1, z2)
  unadjusted_cpch = apply(z, 1, calc_cpch_pvals, true_mu = c(0, 0), oracle = F)
  cpch = sapply(unadjusted_cpch, get_adjusted_pval, M=2, r=2)
  mm_pvals = mediation_test(Z, 0.05)$pval
  #true for alternative, false for null
  #in this setting, this is like the case where 
  alts = (theta1 != 0 & theta2 != 0)
  nulls = !alts
  p.matrix = matrix(c(pvals1, pvals2), nrow = m, byrow = FALSE)
  
  #dhh
  dhh_pvals = (10)*maxp[maxp <= 0.1]
  dhh_alts = alts[maxp <= 0.1]
  dhh_nulls = nulls[maxp <= 0.1]
  dhh_x = x[maxp<= 0.1]
  
  #adafilter
  rejected_ada = adaFilter(p.matrix, r=2, alpha = q)$decision
  tpr_adaf = mean(rejected_ada[alts == TRUE])
  fdr_adaf = sum(rejected_ada[nulls])/max(1, sum(rejected_ada))
  # 
  # Run adapt_gmm
  x <- data.frame(x = x)
  formulas = paste("splines::ns(x,df=", c(2,4), ")")
  alphas <- c(q)
  
  #adapt classical
  print('maxp')
    classical <- adapt_gmm(x=x, pvals=maxp, alphas=alphas,
                           beta_formulas = formulas, model_type = "neural", nclasses= c(2,3))
  rejected_adapt_c = 1*(classical$qvals <= q)
  tpr_adapt_c = mean(rejected_adapt_c[alts])
  fdr_adapt_c = sum(rejected_adapt_c[nulls])/max(1, sum(rejected_adapt_c))
  # 
  #adapt cpch
  
  print('cpch')
  ada_cpch <- adapt_gmm(x=x,pvals=cpch, alphas=alphas,
                          beta_formulas = formulas, model_type = "neural", nclasses= c(2,3))
  rejected_adapt_cpch = 1*(ada_cpch$qvals <= q)
  tpr_adapt_cpch = mean(rejected_adapt_cpch[alts])
  fdr_adapt_cpch =  sum(rejected_adapt_cpch[nulls])/max(1, sum(rejected_adapt_cpch))
  
  #adapt dhh
  if (length(dhh_pvals) >= 1){
    dhh_x = data.frame(x = dhh_x)
    ada_dhh <- adapt_gmm(x=dhh_x, pvals=dhh_pvals, alphas=alphas,
                         beta_formulas = formulas, model_type = "neural", nclasses= c(2,3))
    rejected_adapt_dhh = 1*(ada_dhh$qvals <= q)
    tpr_adapt_dhh = sum(rejected_adapt_dhh[dhh_alts])/sum(alts)
    fdr_adapt_dhh =  sum(rejected_adapt_dhh[dhh_nulls])/max(1, sum(rejected_adapt_dhh))
  } else{
    tpr_adapt_dhh = 0
    fdr_adapt_dhh =  0
  }
  
  #dhh BH and Storey
  if (length(dhh_pvals) >= 1){
    dhh_bh = p.adjust(dhh_pvals, method = 'BH', n = m)
    rejected_dhh_bh = 1*(dhh_bh <= q)
    rejected_dhh_s = 1*(storey(dhh_pvals, q)$significant)
    tpr_dhh_bh = sum( rejected_dhh_bh[dhh_alts])/sum(alts)
    fdr_dhh_bh =  sum( rejected_dhh_bh[dhh_nulls])/max(1, sum( rejected_dhh_bh))
    tpr_dhh_s = sum( rejected_dhh_s[dhh_alts])/sum(alts)
    fdr_dhh_s =  sum( rejected_dhh_s[dhh_nulls])/max(1, sum( rejected_dhh_s))
  } else{
    tpr_dhh_bh = 0
    fdr_dhh_bh =  0
    tpr_dhh_s = 0
    fdr_dhh_s =  0
  }
  # 
  # #bhq classical
  bhqvals_c = p.adjust(maxp, method = 'BH', n = m)
  rejected_bh_c = 1*(bhqvals_c <= q)
  tpr_bh_c = mean(rejected_bh_c[alts])
  fdr_bh_c = sum(rejected_bh_c[nulls])/max(1, sum(rejected_bh_c))
  # 
  
  #bhq cpch
  bhqvals_cpch = p.adjust(cpch, method = 'BH', n = m)
  rejected_bh_cpch = 1*(bhqvals_cpch <= q)
  tpr_bh_cpch = mean(rejected_bh_cpch[alts])
  fdr_bh_cpch = sum(rejected_bh_cpch[nulls])/max(1, sum(rejected_bh_cpch))
  
  # #classical storey p-values
  rejected_storey_c = 1*(storey(maxp, q)$significant)
  tpr_storey_c = mean(rejected_storey_c[alts])
  fdr_storey_c = sum(rejected_storey_c[nulls])/max(1, sum(rejected_storey_c))
  # 
  
  #cpch storey p-values
  rejected_storey_cpch = 1*(storey(cpch, q)$significant)
  tpr_storey_cpch = mean(rejected_storey_cpch[alts])
  fdr_storey_cpch =  sum(rejected_storey_cpch[nulls])/max(1, sum(rejected_storey_cpch))
  # 
  #dact BH
  dact_pvals = DACT(p_a=pvals1, p_b=pvals2, correction = 'JC')
  dact_bh = p.adjust(dact_pvals, method ='BH')
  fdr_bh_dact = sum(dact_bh[nulls] <= q, na.rm = TRUE)/max(1, sum(dact_bh <= q, na.rm = TRUE))
  tpr_bh_dact = mean(dact_bh[alts] <= q)
  
  #dact Storey
  dact_storey = 1*(storey(dact_pvals, q)$significant)
  tpr_storey_dact = mean(dact_storey[alts])
  fdr_storey_dact = sum(dact_storey[nulls])/max(1, sum(dact_storey))
  #
  # #dact adapt-GMM
  dact_adapt_gmm <- adapt_gmm(x=x, pvals=dact_pvals, alphas=alphas,
                              beta_formulas = formulas, model_type = "neural", nclasses= c(3,4))
  rejected_adapt_dact = 1*(dact_adapt_gmm$qvals <= q)
  tpr_adapt_dact = mean(rejected_adapt_dact[alts])
  fdr_adapt_dact = sum( rejected_adapt_dact[nulls])/max(1, sum( rejected_adapt_dact))
  #
  # #hdmt p-values
  alphas= null_estimation(p.matrix)
  alpha00 = alphas$alpha00
  alpha01 = alphas$alpha01
  alpha10 = alphas$alpha10
  alpha1 =  alphas$alpha1
  alpha2 =  alphas$alpha2
  hdmt_loc = fdr_est(alpha00, alpha01, alpha10, alpha1, alpha2, p.matrix, exact = 1)
  fdr_hdmt = sum(hdmt_loc[nulls] <= q)/max(1, sum(hdmt_loc <= q))
  tpr_hdmt = mean(hdmt_loc[alts] <= q)
  
  
  
  #Miles & Chambaz Mediation Test Methods
  adapt_mm <- adapt_gmm(x=x, pvals=mm_pvals, alphas=alphas,
                          beta_formulas = formulas, model_type = "neural", nclasses= c(2,3))
  rejected_adapt_mm = 1*(adapt_mm$qvals <= q)
  tpr_adapt_mm = mean(rejected_adapt_mm[alts])
  fdr_adapt_mm = sum(rejected_adapt_mm[nulls])/max(1, sum(rejected_adapt_mm))

  # #bhq minimax test
  bhqvals_mm = p.adjust(mm_pvals , method = 'BH', n = m)
  rejected_bh_mm = 1*(bhqvals_mm <= q)
  tpr_bh_mm = mean(rejected_bh_mm[alts])
  fdr_bh_mm = sum(rejected_bh_mm[nulls])/max(1, sum(rejected_bh_mm))
  
  # # storey minimax test p-values
  rejected_storey_mm = 1*(storey(mm_pvals, q)$significant)
  tpr_storey_mm = mean(rejected_storey_mm[alts])
  fdr_storey_mm = sum(rejected_storey_mm[nulls])/max(1, sum(rejected_storey_mm))
  
  
  return(c(tpr_adaf, tpr_hdmt, 
           tpr_adapt_c, tpr_adapt_cpch, tpr_adapt_dhh, tpr_adapt_dact,
           tpr_bh_c, tpr_bh_cpch, tpr_bh_dact, tpr_dhh_bh, 
           tpr_storey_c, tpr_storey_cpch, tpr_storey_dact, tpr_dhh_s, tpr_adapt_mm, tpr_bh_mm, tpr_storey_mm,
           fdr_adaf, fdr_hdmt,
           fdr_adapt_c, fdr_adapt_cpch, fdr_adapt_dhh, fdr_adapt_dact,
           fdr_bh_c, fdr_bh_cpch, fdr_bh_dact, fdr_dhh_bh,
           fdr_storey_c, fdr_storey_cpch, fdr_storey_dact, fdr_dhh_s, fdr_adapt_mm, fdr_bh_mm, fdr_storey_mm))
}



thetas = seq(1, 4, 1)
output = sapply(thetas, run_sim)
df = data.frame('it' = iteration, 'theta' =rep(thetas, each = length(output[,1])), 
                'comb_method' = rep(c('AdaFilter', 'HDMT', 'Max-P', 'cPCH', 'DHH', 'DACT',
                                      rep(c('Max-P', 'cPCH', 'DACT', 'DHH'), 2), c('MM', 3)), 2*length(thetas)),
                'mt_method' = rep(c(c('AdaFilter', 'HDMT'), rep('AdaPT-GMM', 4), rep('BH', 4), rep('Storey', 4), c('AdaPT-GMM', 'BH', 'Storey')),  2*length(thetas)),
                'metric' = rep(c(rep('TPR', 17), rep('FDR', 17)), 4),
                'value' = c(output[,1], output[,2], output[,3], output[,4])
)


#change to where you want data to be saved
write.csv(df,file=paste0("Data/iter_", iteration, ".csv"))

