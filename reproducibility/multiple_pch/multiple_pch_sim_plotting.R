library(dplyr)
library(ggplot2)
#set to where data output from running multiple_pch_sim.py is saved
setwd("~/Documents/Research/cpch/multiple_pch/Data/")

files = list.files(path='M_2000_m_4/final_data', full.names = TRUE) 
df <- files %>% 
    lapply(read.csv) %>% 
    bind_rows 

correction_method = c(rep('cPCH',6), rep('Standard', 6), rep('DHH', 6), 'AdaFilter')
comb_method = c(rep(c('Fisher', 'Simes', 'Bonferroni'), 6),'AdaFilter')
mt_method = c(rep(c(rep('BH', 3),rep('Storey', 3)), 3), 'AdaFilter')
df['correction_method'] = factor(rep(correction_method, nrow(df)/19), levels = c( 'AdaFilter', 'Standard','DHH','cPCH'))
df['comb_method'] = rep(comb_method, nrow(df)/19)
df['mt_method'] = factor(rep(mt_method, nrow(df)/19), levels = c( 'AdaFilter', 'BH','Storey'))
labelling = function(string){paste('pi[1]==', string)}
df$pi1 = sapply(df$pi1, labelling)


labelling_1 = function(string){paste('theta ==', string)}
df$ss = sapply(df$ss, labelling_1)
df = df[order(df$correction_method),]
labelling2 = function(string){paste('p ==', string)}
df$p = sapply(df$p, labelling2)

power_df = df[df$type == 'Power',]

ggplot(data=power_df, aes(x=r, y=p_reject, color = correction_method, linetype = comb_method, shape = mt_method)) +
  geom_line(size = 1)+ geom_point(size = 3) + 
  scale_x_continuous(breaks = seq(2, 4, 1)) +
  theme_light() +
  theme(strip.text = element_text(colour = 'black')) +
  theme(strip.background =element_rect(fill="lightgray")) + 
  facet_wrap(~pi1 + ss + p, ncol = 6, dir = 'h', labeller = labeller(.cols = label_parsed, .multi_line = FALSE)) +
  theme(panel.border = element_blank()) + scale_shape_manual(values=c(3, 16, 17)) + 
 scale_linetype_manual(values=c(4, 2, 1, 3)) + 
  theme(aspect.ratio = 1, strip.text.x = element_text(size = 12),
  axis.title =element_text(size=15), axis.text = element_text(size = 12),
  legend.title=element_text(size=15), legend.text=element_text(size=12))+
  labs(shape = "Multiple Testing\n Method", linetype = 'Combination Method', color = "Correction Method") +
  ylab('Power')  
  
final_plot =  ggplot(data=power_df, aes(x=r, y=p_reject, color = correction_method, linetype = comb_method, shape = mt_method)) +
  geom_line(size = 1)+ geom_point(size = 2) +
  scale_x_continuous(breaks = seq(2, 4, 1)) +
  theme_light() +
  theme(strip.text = element_text(colour = 'black')) +
  theme(strip.background =element_rect(fill="lightgray")) + 
  facet_wrap(~pi1 + ss + p, ncol = 6, dir = 'h', labeller = labeller(.cols = label_parsed, .multi_line = FALSE)) +
  theme(panel.border = element_blank()) + scale_shape_manual(values=c(3, 16, 17)) + 
  scale_linetype_manual(values=c(4, 2, 1, 3)) + 
 theme(aspect.ratio = 1, #strip.text.x = element_text(size = 12),
       axis.title =element_text(size=15), axis.text = element_text(size = 13),
        legend.title=element_text(size=15), legend.text=element_text(size=13))+
  labs(shape = "Multiple Testing\n Method", linetype = 'Combination Method', color = "Correction Method") +
  ylab('Power') 
  
#set directory to where you want plots to be saved
ggsave(filename = "mt_power_plot_m_4.png",
       plot = final_plot, path = '~/Documents/Research/cpch/Plots/', bg = 'white',
       height = 8, width = 12)


fdr_df = df[df$type == 'FDR',]

final_plot_fdr = ggplot(data=fdr_df, aes(x=r, y=p_reject, color = correction_method, linetype = comb_method, shape = mt_method)) +
  geom_line(size = 1)+ geom_point(size = 2) +
  theme_light() +
  scale_x_continuous(breaks = seq(2, 4, 1)) +
  theme(strip.text = element_text(colour = 'black')) +
  theme(strip.background =element_rect(fill="lightgray")) + 
  facet_wrap(~pi1 + ss + p, ncol = 8, labeller = labeller(.cols = label_parsed, .multi_line = FALSE)) +
  theme(panel.border = element_blank()) + scale_shape_manual(values=c(3, 16, 17)) + 
  scale_linetype_manual(values=c(10, 5, 1, 3)) + 
  theme(aspect.ratio = 1,# strip.text.x = element_text(size = 11),
        axis.title =element_text(size=15), axis.text = element_text(size = 12),
        legend.title=element_text(size=15), legend.text=element_text(size=12)) +
  labs(shape = "Multiple Testing\nMethod", linetype = 'Combination Method', color = "Correction Method") +
  geom_hline(aes(yintercept = 0.1), linetype = 2, color = 'black') +
  ylab('FDR')

#set directory to where you want plots to be saved
ggsave(filename = "mt_fdr_plot_m_4.png",
       plot = final_plot_fdr, path = 'Plots', bg = 'white',
       height = 8, width = 12)


 