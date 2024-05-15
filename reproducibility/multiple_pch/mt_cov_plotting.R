library(ggplot2)
library(dplyr)
library(readr)
#set directory to where data from mt_covariates_sim.R is saved
setwd("PATH")

#set to name of data file saved from running compile_csvs_mt_cov.R tocompile csv output from mt_covariates_sim.R 
final_df = read.csv('mt_cov_data_10000_100_reps.csv')
#Rename TPR to Power
final_df$metric = rep(c('FDR', 'Power'), nrow(final_df)/2)


ggplot(final_df, aes(x = theta, y = value, linetype = mt_method, color= comb_method)) + geom_line() + 
  theme_minimal() + facet_wrap(~metric) + labs(color = "Combination\nMethod",  linetype = 'Multiple Testing\nMethod')  + 
  scale_linetype_manual(values=c(6, 1, 4, 2, 3)) + geom_hline(yintercept = 0.1, linetype = 'dotted') +
  ylab('Power') + xlab(expression(paste(theta)))+ 
  geom_errorbar(aes(ymin=value-2*se, ymax=value+2*se), width=.05, linetype = 1) +
  theme(aspect.ratio = 1, strip.text.x = element_text(size =14),
        axis.title =element_text(size=14), axis.text = element_text(size = 12),
        legend.title=element_text(size=14), legend.text=element_text(size=12), title = element_text(size = 12))
