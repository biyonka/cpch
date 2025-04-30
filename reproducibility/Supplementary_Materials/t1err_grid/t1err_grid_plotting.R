#plotting for t1err grid for m=3, r=3, r^*=2 

#change directory to where data output from robustness_sim.py and single_pch_sim.py is saved
#setwd("~/Documents/Research/cpch/single_pc_power/Data_t1err_grid") 

m = 3
files = list.files(path="~/Documents/Research/cpch/single_pc_power/Data_t1err_grid", full.names = TRUE) 
t1err_grid_df <- files %>% 
  lapply(read.csv) %>% 
  bind_rows




files_1 = list.files(path="~/Documents/Research/cpch/single_pc_power/Data_t1err_grid_1", full.names = TRUE) 
t1err_grid_df_1 <- files_1 %>% 
  lapply(read.csv) %>% 
  bind_rows

t1err_grid_df$total_reject = t1err_grid_df_1$total_reject + t1err_grid_df$total_reject




prob_rej = sapply(seq(1, nrow(t1err_grid_df), 2), function(i){
  t1 = t1err_grid_df$theta1[i]
  t2 = t1err_grid_df$theta2[i]
  if(t1==t2){
    return(t1err_grid_df[t1err_grid_df$theta1==t1 & t1err_grid_df$theta2==t2,]$total_reject/(4*10000))
  }else{
    total_reject = t1err_grid_df[t1err_grid_df$theta2==t1 & t1err_grid_df$theta1==t2,]$total_reject +t1err_grid_df[t1err_grid_df$theta2==t2 & t1err_grid_df$theta1==t1,]$total_reject 
    return(total_reject/(4*10000))
  }
})

# t1=0.2
# t2=0.4
# t1err_grid_df[t1err_grid_df$theta2==t1 & t1err_grid_df$theta1==t2,]

t1err_grid_df['p_reject'] =as.vector(prob_rej)


cbbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

min_val = min(t1err_grid_df['p_reject'])
max_val = max(t1err_grid_df['p_reject'])


b=c(min_val, 0.05, 0.055, 0.06, 0.0625, max_val)

heatmap = ggplot() +
   theme_minimal() +
  # First layer for Max-p
  geom_tile(data = subset(t1err_grid_df, method ="cPCH-Fisher"), 
            aes(x = theta1, y = theta2, fill=p_reject))+ 
  scale_fill_gradientn(limits = c(min_val, max_val),
                       colours=c("white","lightgray",  "gray70", "gray55", 'black'),
                      # colours=c("#0072B2","#56B4E9",  "#009E73", "#CC79A7", 'black'),
                      # colours=c("navyblue",'darkblue', "darkmagenta", "red", 'black'),
                       breaks=b, labels=format(b)) +
  labs(
    x = expression(theta[1]),
    y = expression(theta[2])
  ) + 
   theme(legend.position = "none", aspect.ratio = 1,
            axis.title =element_text(size=17), axis.text = element_text(size = 15)) 
                  
heatmap



t1err_fisher = t1err_grid_df[t1err_grid_df$method == 'cPCH-Fisher',]


# Function to aggregate heatmap data by reducing resolution
reduce_heatmap_resolution <- function(data, x_col, y_col, value_col, bin_size) {
  # Bin the x and y columns
  data <- data %>%
    mutate(
      x_bin = floor(!!sym(x_col) / bin_size) * bin_size,# + bin_size / 2,
      y_bin = floor(!!sym(y_col) / bin_size) * bin_size# + bin_size / 2
    )
  
  # Aggregate values by x_bin and y_bin
  aggregated_data <- data %>%
    group_by(x_bin, y_bin) %>%
    summarize(p_reject = mean(!!sym(value_col)), .groups = "drop")
  
  return(aggregated_data)
}


# testing = t1err_fisher %>%
#   mutate(
#     x_bin = floor(!!sym('theta1') / 0.2) * 0.2,# + 0.2 / 2,
#     y_bin = floor(!!sym('theta2') / 0.2) * 0.2# + 0.2 / 2
#   )


smaller_fisher =reduce_heatmap_resolution(t1err_fisher, "theta1", "theta2", "p_reject", 0.2)

# averages <- sapply(seq(1, length(t1err_fisher$p_reject) - 1, by = 2), function(i) {
#   mean(t1err_fisher$p_reject[i:(i + 1)])
# })

min_val = min(smaller_fisher['p_reject'])
max_val = max(smaller_fisher['p_reject'])
b=c(0, 0.05, 0.055, 0.06, 0.061, 0.06325, max_val, max_val)


smaller_fisher_1 = smaller_fisher[smaller_fisher$x_bin <= 4.6 &smaller_fisher$y_bin <= 4.6, ]
heatmap_smaller_res = ggplot() +
  theme_minimal() +
  # First layer for Max-p
  geom_tile(data = smaller_fisher_1, 
            aes(x = x_bin, y = y_bin, fill=p_reject))+ 
  scale_fill_gradientn(limits = c(min_val, max_val),
                       colours=c("white","lightgray", "gray70", "gray60",'gray50','gray45',  'gray12','black'),
                       # colours=c("#0072B2","#56B4E9",  "#009E73", "#CC79A7", 'black'),
                       # colours=c("navyblue",'darkblue', "darkmagenta", "red", 'black'),
                       breaks=b, labels=format(b)) +
  labs(
    x = expression(theta[1]),
    y = expression(theta[2])
  ) + 
  theme(legend.position = "none", aspect.ratio = 1,
        axis.title =element_text(size=17), axis.text = element_text(size = 15)) 

heatmap_smaller_res 

#scale_fill_gradient(low = "white", high = "black")# + 
# scale_fill_distiller(palette = "YlGn", direction=2, 
#   values = seq(0, 0.01, 0.1))


#change directory to where you want plot to be saved
ggsave(filename = paste0("heatmap_smaller_res.eps"),
       plot = heatmap_smaller_res, path = '~/Documents/Research/cpch/Plots/', bg = 'white',
       height = 7, width = 7)

