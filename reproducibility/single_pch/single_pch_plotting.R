library(dplyr)
library(ggplot2)
#change directory to where data output from robustness_sim.py and single_pch_sim.py is saved
setwd("~/Documents/Research/cpch/single_pch/Data/") 

m = 4
files = list.files(path=paste0('m_',m), full.names = TRUE) 
df <- files %>% 
  lapply(read.csv) %>% 
  bind_rows


correction_method = rep(c(rep('Standard', 3),rep('cPCH', 3), rep('cPCH Oracle', 3)), nrow(df)/9)
comb_method = rep(c('Fisher', 'Simes', 'Bonferroni'), nrow(df)/3)

df['correction_method'] = correction_method
df['correction_method'] = factor(correction_method, levels = c('Standard','cPCH','cPCH Oracle'))
df['comb_method'] = comb_method
#sort df columns so cpch methods are
df = df[order(match(df$correction_method, c('Standard','cPCH', 'cPCH Oracle'))),]

df_null = df[df$k < df$r,]
df_alt = df[df$k >= df$r,]

labelling = function(string){paste0("r =",  string)}

df_null$r = sapply(df_null$r, labelling)
df_alt$r = sapply(df_alt$r, labelling)
df$r = sapply(df$r, labelling)
labelling_1 = function(string){paste0('r* =', string)}
df_null$k = sapply(df_null$k, labelling_1)
df_alt$k = sapply(df_alt$k, labelling_1)
df$k = sapply(df$k, labelling_1)

single_pc_alt = ggplot(data=df_alt, aes(x=sig_strength, y=p_reject, color = correction_method, linetype = comb_method)) +
  geom_line(size = 0.75)+ 
  theme_light() +
  theme(strip.text = element_text(colour = 'black')) +
  geom_hline(aes(yintercept = 0.05), linetype = 4, color = 'gray42') +
  theme(strip.background =element_rect(fill="lightgray")) + 
  facet_wrap(~r + k, nrow = 2, dir = 'h', labeller = labeller(.multi_line = FALSE)) +
  theme(panel.border = element_blank()) + 
  theme(aspect.ratio = 1, strip.text.x = element_text(size = 15),
        axis.title =element_text(size=15), axis.text = element_text(size = 12),
      legend.title=element_text(size=15), legend.text=element_text(size=12))+
  labs(linetype = 'Combination Method', color = "Correction Method") +
  ylab('Power') +xlab(expression(paste(theta))) +
 geom_errorbar(aes(ymin=p_reject-2*se, ymax=p_reject+2*se), show.legend = FALSE, width=.1)


null_no_orac = df_null[(df_null$correction_method != 'cPCH Oracle'),]
single_pc_null = ggplot(data=null_no_orac, aes(x=sig_strength, y=p_reject, color = correction_method, linetype = comb_method)) +
  geom_errorbar(aes(ymin=p_reject-2*se, ymax=p_reject+2*se), width=.1, show.legend = FALSE) + geom_line(size = 0.75)+ #geom_point(size = 1.5) +
  theme_light() +
  theme(strip.text = element_text(colour = 'black')) +
  geom_hline(aes(yintercept = 0.05), linetype = 4, color = 'gray42') +
  theme(strip.background =element_rect(fill="lightgray")) + 
  facet_wrap(~r + k, nrow = 2, dir = 'h', labeller = labeller(.multi_line = FALSE)) +
  theme(panel.border = element_blank()) + scale_shape_manual(values=c(3, 16, 17)) + 
  theme(aspect.ratio = 1, strip.text.x = element_text(size = 15),
        axis.title =element_text(size=15), axis.text = element_text(size = 12),
        legend.title=element_text(size=15), legend.text=element_text(size=12))+
  labs(linetype = 'Combination Method', color = "Correction Method") +
  ylab('Type I Error') +xlab(expression(paste(theta))) +
  ylim(0, 0.1)

single_pc = ggplot(data=df, aes(x=sig_strength, y=p_reject, color = correction_method, linetype = comb_method)) +
  geom_line(size = 0.75)+ geom_point(aes( fill = correction_method), size = 1.5) +
  theme_light() +
  theme(strip.text = element_text(colour = 'black')) +
  geom_hline(aes(yintercept = 0.05), linetype = 4, color = 'gray42') +
  theme(strip.background =element_rect(fill="lightgray")) + 
  facet_wrap(~r + k, nrow = 2, dir = 'h', labeller = labeller(.multi_line = FALSE)) +
  theme(panel.border = element_blank()) + scale_shape_manual(values=c(3, 16, 17)) + 
  scale_linetype_manual(values=c(5, 1, 3)) + 
  theme(aspect.ratio = 1, strip.text.x = element_text(size = 15),
        axis.title =element_text(size=15), axis.text = element_text(size = 12),
        legend.title=element_text(size=12), legend.text=element_text(size=12))+
  labs(linetype = 'Combination Method', color = "Correction Method") +
  ylab('TPR') + xlab(expression(paste(theta)))
  
#change directory to where you want plot to be saved
ggsave(filename = paste0("single_pc_plot_m_", m,".png"),
       plot = single_pc, path = 'Plots', bg = 'white',
       height = 7, width = 12)
      
#change directory to where you want plot to be saved
ggsave(filename = paste0("single_pc_alt_plot_m_", m,".png"),
       plot = single_pc_alt, path = 'Plots', bg = 'white',
       height = 7, width = 12)
#change directory to where you want plot to be saved
ggsave(filename = paste0("single_pc_null_plot_m_", m,".png"),
       plot = single_pc_null, path = '~/Documents/Research/cpch/Plots/', bg = 'white',
       height = 7, width = 14)


###########################   
#### ROBUSTNESS PLOTS #####
########################### 
m=3
#change directory to where data is saved
setwd("~/Documents/Research/cpch/single_pch/Data/robustness")
rob_files = list.files(path=paste0('m_',m), full.names = TRUE) 
r_df <- rob_files %>% 
  lapply(read.csv) %>% 
  bind_rows 

methods = c(rep('CPCH-Fisher-t',5), rep('CPCH-Simes-t', 5), rep('CPCH-Bonferroni-t', 5), rep('CPCH-Fisher-norm', 5), rep('CPCH-Simes-norm', 5), rep('CPCH-Bonferroni-norm', 5)) 
r_df['methods'] = rep(methods, nrow(r_df)/30)
distributional_assumption = factor(rep(c(rep('t', 15), rep('norm', 15)), nrow(r_df)/30), levels = c('t', 'norm'))
comb_method = factor(rep(c(rep('Fisher', 5), rep( 'Simes', 5), rep('Bonferroni', 5)), nrow(r_df)/15), levels = c('Bonferroni', 'Simes', 'Fisher'))
r_df['distribution'] = distributional_assumption 
r_df['comb_method'] = comb_method

r_df_3 = r_df[r_df$sig_strength == 4,]
r_df_3 = r_df_3[order(r_df_3$comb_method),]

r_df_null_3 = r_df_3[r_df_3$k < r_df_3$r,]
r_df_alt_3 = r_df_3[r_df_3$k >= r_df_3$r,]

labelling = function(string){paste('r=', string)}
r_df_null_3$r = sapply(r_df_null_3$r, labelling)
r_df_alt_3$r = sapply(r_df_alt_3$r, labelling)
labelling_1 = function(string){paste0('r* =', string)}
r_df_null_3$k = sapply(r_df_null_3$k, labelling_1)
r_df_alt_3$k = sapply(r_df_alt_3$k, labelling_1)



null = ggplot(data=r_df_null_3, aes(x=dof, y=p_reject, linetype = comb_method,  color = distribution)) +
  geom_line(size = 0.75)+ 
  theme_light() +
  theme(strip.text = element_text(colour = 'black')) +
  geom_hline(aes(yintercept = 0.05), linetype = 4, color = 'gray42') +
  theme(strip.background =element_rect(fill="lightgray")) + 
  facet_wrap(~r + k, ncol = 5, dir = 'h', labeller = labeller(.multi_line = FALSE)) +
  theme(panel.border = element_blank()) + scale_shape_manual(values=c(3, 16, 17)) + 
  scale_linetype_manual(values=c(3, 5, 1)) + 
  theme(aspect.ratio = 1, strip.text.x = element_text(size = 15),
        axis.title =element_text(size=15), axis.text = element_text(size = 12),
        legend.title=element_text(size=12), legend.text=element_text(size=11))+
  labs(linetype = "Combination Method", color = 'Distributional Assumption') +
  ylab('Type I Error') +scale_x_continuous(name = 'DOF', breaks = c(1, 3, 5, 7, 9))+  ggtitle('Null Configurations') +
  geom_errorbar(aes(ymin=p_reject-2*se, ymax=p_reject+2*se), show.legend = FALSE, width=.1, linetype = 1)  



alt = ggplot(data=r_df_alt_3, aes(x=dof, y=p_reject,  linetype = comb_method,  color = distribution)) +
  geom_line(size = 0.75)+
  theme_light() +
  theme(strip.text = element_text(colour = 'black')) +
  geom_hline(aes(yintercept = 0.05), linetype = 4, color = 'gray42') +
  theme(strip.background =element_rect(fill="lightgray")) + 
  facet_wrap(~r + k, ncol = 5, dir = 'h', labeller = labeller(.multi_line = FALSE)) +
  theme(panel.border = element_blank()) + scale_shape_manual(values=c(3, 16, 17)) + 
  scale_linetype_manual(values=c(3, 5, 1)) + 
  theme(aspect.ratio = 1, strip.text.x = element_text(size = 15),
        axis.title =element_text(size=15), axis.text = element_text(size = 12),
        legend.title=element_text(size=12), legend.text=element_text(size=11))+
  labs(linetype = "Combination Method", color = 'Distributional Assumption') + ggtitle('Alternative Configurations') +
  ylab('Power') +scale_x_continuous(name = '', breaks = c(1, 3, 5, 7, 9))+ 
  geom_errorbar(aes(ymin=p_reject-2*se, ymax=p_reject+2*se), width=.1, show.legend = FALSE) 


#change directory to where you want plot to be saved
ggsave(filename = "robustness_plot_alt.png",
       plot =alt, path = 'Plots', bg = 'white',
       height = 7, width = 15, dpi = 500)


#change directory to where you want plot to be saved
ggsave(filename = "robustness_plot_null.png",
       plot = null, path = 'Plots', bg = 'white',
       height = 7, width = 15, dpi = 500)

library(ggpubr)
final_robust_plot = ggarrange(alt, null, ncol=1, nrow=2, common.legend = TRUE, legend="right", heights = c(1, 1), widths = c(1, 1))

#change directory to where you want plot to be saved
ggsave(filename = "robustness_plot.png",
       plot = final_robust_plot, path = 'Plots', bg = 'white',
       height = 8.71, width = 12.31, dpi = 500)


