library(ggplot2)
library(dplyr)
library(readr)

# First read in the arguments listed at the command line
args=(commandArgs(TRUE))# args is now a list of character vectors
# First check to see if arguments are passed.
# Then cycle through each element of the list and evaluate the expressions.
arguments = c('iter')

if (length (args) == 0) {
  print("No arguments supplied.")
  ## supply default values
  iter = 1
} else {
  for (i in 1:length(args)) {
    eval (parse (text = paste(arguments[i]," = ",args[[i]],sep="") ))
  }
}

total_iterations = 100 #change this to your specifications 

mypath = paste0('Data', sep ="")

files = list.files(path=mypath, full.names = TRUE) 
if (length(files) == total_iterations){
  df <- files %>% 
    lapply(read_csv) %>% 
    bind_rows 
}

df = df[,-1 ]

#get mean of value grouped by (theta mt_method comb_method metric)
means = df %>% group_by(theta, mt_method, comb_method, metric) %>%  summarise(across(-it, mean, na.rm = TRUE))
#get sd of value grouped by (theta mt_method comb_method metric)
ses = df %>% group_by(theta, mt_method, comb_method, metric) %>%  summarise(across(-it, function(x){sd(x, na.rm = TRUE)/(total_iterations ** 0.5)}))
#get sd of value grouped by (theta mt_method comb_method metric)

means['se'] = ses['value']


#replace path with where you want data to be saved
write.csv(means, paste0('Data/mt_cov_data_10000_',total_iterations,'_reps.csv'))

