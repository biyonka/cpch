library(ggplot2)
library(dplyr)
library(readr)
setwd("~/Documents/Research/minmax/cpch/validity")
#setwd("~/Documents/Research/cpch/validity") #set to location of validity data


#data used to reproduce the plots in the Approximate Validity section of the paper
combined = read.csv('Data/combined_data_m_3_r_3.csv')

#make side-by-side facet plots
get_nonzero_mu = function(mu_vec){
  substring = tail(strsplit(mu_vec, split = " ")[[1]], 1)
  substr(substring, 1, nchar(substring)-1)
}

nonzero_mus = sapply(combined$mu, get_nonzero_mu)
combined$nonzero_mus = nonzero_mus


#qqplot m=3, r=2, m=r=2
qq_3_2 = combined[combined$m == 3 & combined$r == 2, ]
labelling = function(string){paste0('m = ', string, ', r = 2')}
qq_3_2$n_lab = sapply(qq_3_2 $m, labelling)


qq_2_2 = combined[combined$m == 2 & combined$r == 2, ]
qq_2_2$n_lab = 'r = m = 2'

test = data.frame('nonzero_mus'= c(qq_2_2$nonzero_mus, rep(qq_3_2$nonzero_mus, 3)), 
                  'n_lab'= c(rep('r = m = 2', nrow(qq_2_2)), rep(qq_3_2$n_lab, 3)), 
                  'Combination'= c(rep('Bonferroni', nrow(qq_2_2)), rep('Fisher', nrow(qq_3_2)), rep('Simes', nrow(qq_3_2)), rep('Bonferroni', nrow(qq_3_2))),
                  'Quantiles'= c(qq_2_2$cpch_f, qq_3_2$cpch_f, qq_3_2$cpch_s, qq_3_2$cpch_b)
)

test$n_lab = factor(test$n_lab, levels = c('r = m = 2', 'm = 3, r = 2'))
test$Combination = factor(test$Combination, levels = c('Bonferroni', 'Simes', 'Fisher'))
alltogether = ggplot(test, aes(sample = Quantiles, color = nonzero_mus)) + 
  geom_qq(aes(color = nonzero_mus), size = 0.3, distribution = stats::qunif) +
  facet_wrap(~ n_lab + Combination, nrow = 1) +
  geom_abline(slope = 1, intercept = 0) +
  theme_light() + 
  theme(strip.text = element_text(colour = 'black')) +
  theme(strip.background =element_rect(fill="lightgray")) + 
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  labs('color' = expression(theta[(m)])) + xlab('Theoretical Unif[0, 1]') + 
  ylab('Empirical cPCH')+ scale_x_continuous(breaks=c(0,0.5,1)) +
  theme(aspect.ratio = 1,strip.text.x = element_text(size = 35),
        axis.title =element_text(size= 37), axis.text = element_text(size = 35),
        legend.title=element_text(size = 40), legend.text=element_text(size=30), 
        panel.border = element_blank(), panel.spacing.x = unit(6, "mm"))

#Set directory to where you want plots to be saved
ggsave(filename = "cpch_qqplot_all.png",
       plot = alltogether, path = 'Plots',
       bg = 'white', height = 7, width = 24)


 #KS plots
ks_unifs = sapply(seq(1, 10000, 1), function(i){
  1-ks.test(runif(10000), runif)[[1]]
}
)

ks_3_2 = read.csv('Data/ks_df_m_3_r_2.csv')
ks_2_2 = read.csv('Data/ks_df_m_2_r_2.csv')

methods22 = c(rep('Fisher', nrow(ks_2_2)), rep('Simes', nrow(ks_2_2)), rep('Bonferroni', nrow(ks_2_2)))
methods32 = c(rep('Fisher', nrow(ks_3_2)), rep('Simes', nrow(ks_3_2)), rep('Bonferroni', nrow(ks_3_2)))
ks22 = c(ks_2_2$ks_f,  ks_2_2$ks_s, ks_2_2$ks_b)
ks32 = c(ks_3_2$ks_f,  ks_3_2$ks_s, ks_3_2$ks_b)


ks_df = data.frame(config = c(rep('m = 2, r = 2', length(ks22)), rep('m = 3, r = 2', length(ks32))),
                   Method = c(methods22, methods32),
                   nonzero_mu = rep(seq(0, 4, 0.2), 2),
                   ks = c(ks22, ks32)
)


ks_plot = ggplot(ks_df, aes(x = nonzero_mu, y = ks, color = Method)) +
  geom_line(size = 1) + ylim(c(0, max(0.2))) + 
  theme_light() +
  theme(strip.text = element_text(colour = 'black')) +
  theme(strip.background =element_rect(fill="lightgray")) + 
  geom_hline(aes(yintercept = mean(ks_unifs)), linetype = 2, color = 'red') +
  # geom_hline(aes(yintercept = mean(ks_unifs), linetype = 'Avg K-S Distance \n between 10000 Independent Samples \n from Unif[0, 1]'), color = 'red') +
  theme(legend.key.height=unit(2, "cm")) + 
  facet_wrap(~ config, ncol = 2) +
  scale_linetype_manual(name = "", values = c(2, 2), guide = guide_legend(label = TRUE)) +
  xlab(expression(theta[(m)])) + ylab('K-S Distance') + labs(color = expression(paste('Combining \n Function'))) +
  theme(aspect.ratio = 1,strip.text.x = element_text(size = 22),
        axis.title =element_text(size= 23), axis.text = element_text(size = 21),
        legend.title=element_text(size = 21), legend.text=element_text(size=20), 
        panel.border = element_blank())

#Set directory to where you want plots to be saved
ggsave(filename = "cpch_ks_plot_all.png",
       plot = ks_plot, path = 'Plots', bg = 'white', 
       height = 7, width = 14)



#qqplot m=r=3
qq_3_3 = combined[combined$m == 3 & combined$r == 3, ]
qq_3_3$n_lab = 'r=m=3'

#make side-by-side facet plots
get_nonzero_mu = function(mu_vec){
  substring = tail(strsplit(mu_vec, split = " ")[[1]], 2)
  x = str(substr(substring, 1, nchar(substring)-1))
  final = paste0('(', x[1], ',', x[2], ')')
  return(final)
}

nonzero_mu = sapply(qq_3_3$mu, get_nonzero_mu)
qq_3_3$nonzero_mus = nonzero_mu

qq_33 = ggplot(qq_3_3, aes(sample = cpch_f, color = mu)) + 
  geom_qq(aes(color = mu), size = 0.3, distribution = stats::qunif) +
  theme_light() +
  theme(strip.text = element_text(colour = 'black')) +
  theme(strip.background =element_rect(fill="lightgray")) + 
  geom_abline(slope = 1, intercept = 0) +
  # guides(colour = guide_legend(override.aes = list(size=2))) + 
  labs('color' = expression(theta))  + xlab('Theoretical Unif[0, 1]') + 
  ylab('Empirical cPCH') 

#Set directory to where you want plots to be saved
ggsave(filename = "qq_3_3.png",
       plot = qq_33, path = 'Plots', bg = 'white', 
       height = 7, width = 20)
