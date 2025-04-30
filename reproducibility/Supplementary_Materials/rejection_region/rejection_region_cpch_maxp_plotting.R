
#import things
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
#calculate cpch p-values
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


#make grid
xs = seq(0, 5, 0.01)
ys = xs #+ 0.0001

df = expand.grid(xs, ys)
#df = data.frame(t1=xs, t2=ys)

alpha = 0.05

#apply function to rows of data frame
test = apply(df, 1, function(row){
  signal_cpch = rep(NA, nrow(df))#matrix(NA,length(xs),length(ys))
  signal_maxp = rep(NA, nrow(df))#matrix(NA,length(xs),length(ys))
  
  pval_maxp =  2 * max(pnorm(abs(row[1]), 0, 1, lower.tail = FALSE), pnorm(abs(row[2]), 0,1, lower.tail = FALSE))

  #minmax
  pval_cpch = calc_cpch_pvals(row, true_mu=c(0, 0))
  
return(c(pval_maxp <= alpha, 
         pval_cpch <= 0.0425,
         pval_cpch <= alpha))
 # signal_cpch[j,k] = pval_cpch <= alpha
  #signal_maxp[j,k] = pval_maxp <= alpha
})



final = data.frame('T1' = rep(df[,1], 3), 'T2' = rep(df[,2], 3),
'method' = c(rep('Max-p', length(test[1,])), rep('cPCH', length(test[2,])), rep('Unadjusted cPCH', length(test[3,]))), 
'Rejection'=c(test[1,], test[2,], test[3,])
)

final$Rejection[is.na(final$Rejection)] == TRUE

# The palette with black:
cbbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

write.csv(final, '~/Documents/Research/cpch/Data/rejection_regions_data.csv')

#final = read.csv( '~/Documents/Research/cpch/Data/rejection_regions_data.csv')
# Plotting rejection regions with transparency to show overlapping regions
rejection_region_plot = ggplot() +
  # First layer for Max-p
  geom_tile(data = subset(final, method == "Max-p" & Rejection == 1), 
            aes(x = T1, y = T2), 
            fill ="#0072B2", alpha = 1
            ) +
  # # Second layer for cPCH with transparency so Max-p is visible underneath
  # geom_tile(data = subset(final, method == "cPCH" & Rejection == 1),
  #           aes(x = T1, y = T2),
  #           fill = "#CC79A7", alpha = 0.5) +
  # Second layer for cPCH with transparency so Max-p is visible underneath
  geom_tile(data = subset(final, method == "Unadjusted cPCH" & Rejection == 1),
            aes(x = T1, y = T2),
            fill =  "#D55E00", alpha = 0.6) +
  labs(title = "",
       x = expression(T[1]),
       y = expression(T[2])) +
  theme_minimal() + scale_x_continuous(expand=c(0, 0)) +scale_y_continuous(expand=c(0, 0)) +
  theme(legend.position = "none", aspect.ratio = 1,
        axis.title =element_text(size=17), axis.text = element_text(size = 15)) + coord_cartesian() 


#a tiny strip where T_1 = 1.96, T_2 >= 2ish, and vice versa, which is in  Max-P's rejection region but not unadjusted cPCH
rejection_region_plot


#change directory to where you want plot to be saved
ggsave(filename = paste0("rejection_region_plot_unadjusted_cpch.eps"),
       plot = rejection_region_plot, path = '~/Documents/Research/cpch/Plots/', bg = 'white',
       height = 7, width = 7)






# 
# ggplot(pval_ses) +  geom_point(aes(x=mean, y=se, color=method)) +
#   scale_x_log10(
#   breaks = scales::trans_breaks("log10", function(x) 10^x),
#   labels = scales::trans_format("log10", scales::math_format(10^.x))
# ) +
#   scale_y_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   )

# 
# + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")#ylim(c(0,0.1e-01)) #+  xlim(c(0,0.1))

# 
# 
# 
# ggplot(final, aes(x = T1, y = T2)) +
#   geom_tile(aes(fill = interaction(method, Rejection)), alpha = 0.3) + scale_fill_brewer() +  
# #  labs(title = "Rejection regions") +
#   #geom_text(aes(label = Power), size = 2) +
#   #facet_wrap(~ method, nrow = 2)+
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1), legend.position='none')# +
#  # scale_x_discrete(name = 'x', breaks=c('-5','-4', '-3', '-2', '-1', "0","1","2", '3', '4', '5')) +
#  # scale_y_discrete(name = 'y', breaks=c('-5','-4', '-3', '-2', '-1', "0","1","2", '3', '4', '5'))
# 








# 
# # Load necessary libraries
# library(ggplot2)
# library(dplyr)
# 
# data = final %>% filter(method == "cPCH")
# # Define a function to find boundary points
# find_boundaries <- function(data) {
#   # Create shifted versions of the data to check neighbors
#   right <- data %>% mutate(T1 = T1 + 0.05, direction = "right")
#   left <- data %>% mutate(T1 = T1 - 0.05, direction = "left")
#   up <- data %>% mutate(T2 = T2 + 0.05, direction = "up")
#   down <- data %>% mutate(T2 = T2 - 0.05, direction = "down")
#   
#   # Combine original data with shifted data to check for boundary conditions
#   neighbors <- bind_rows(right, left, up, down) %>%
#     select(T1, T2, method, direction) %>%
#     left_join(data, by = c("T1", "T2", "method")) %>%
#     mutate(is_boundary = if_else(is.na(Rejection) | Rejection == 0, TRUE, FALSE))
#   
#   test = neighbors[!is.na(neighbors$Rejection),]
#   return(
#          
#          test[test$is_boundary == TRUE,]
#          
#         # $& neighbors$Rejection == TRUE,
#          )
#   # Filter the original points where at least one neighbor is not in the rejection region
#   boundaries <- data %>%
#     filter(Rejection == TRUE) %>%
#     left_join(neighbors %>% filter(is_boundary == TRUE), by = c("T1", "T2", "method")) %>%
#     filter(!is.na(direction)) %>%
#     select(T1, T2, method) %>%
#     distinct()
#   
#   return(boundaries)
# }
# 
# # Find boundaries for each method
# cPCH_boundary <- find_boundaries(final %>% filter(method == "cPCH"))
# Maxp_boundary <- find_boundaries(final %>% filter(method == "Max-p"))
# 
# # Plotting the rejection region boundaries
# ggplot() +
#   # Boundary for Max-p rejection region
#   geom_path(data = Maxp_boundary, aes(x = T1, y = T2), color = "blue", size = 1, linetype = "dashed") +
#   # Boundary for cPCH rejection region
#   geom_path(data = cPCH_boundary, aes(x = T1, y = T2), color = "red", size = 1, linetype = "solid") +
#   labs(title = "Outline of Rejection Regions for 'cPCH' and 'Max-p' Methods",
#        x = "T1",
#        y = "T2") +
#   theme_minimal() +
#   theme(legend.position = "none")
# 




