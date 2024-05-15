library(ggplot2)
library(dplyr)
library(readr)

# First read in the arguments listed at the command line
args=(commandArgs(TRUE))# args is now a list of character vectors
# First check to see if arguments are passed.
# Then cycle through each element of the list and evaluate the expressions.
arguments = c('m', 'M', 'r', 'pi1', 'p', 'sig_strength')

if (length (args) == 0) {
  print("No arguments supplied.")
  ## supply default values
  m = 2
  M = 2000
  r = 2
  p = 0.0
  pi1 = 0.2
  sig_strength = 2
} else {
  for (i in 1:length(args)) {
    eval (parse (text = paste(arguments[i]," = ",args[[i]],sep="") ))
  }
}

power_iters = 20
total_iterations = 50

#Change mypath to where corresponding data from multiple_pch_sim.py is saved
#because R will paste float of form 1.0 as just 1 and the data generated from the cluster is 1.0
if (p == 1){
  mypath = paste('Data/M_', M, '_m_', m, '/r_', r, '_pi1_', pi1, '_p_1.0_ss_', sig_strength, sep ="")
} else if (p == 0){
  mypath = paste('Data/M_', M, '_m_', m, '/r_', r, '_pi1_', pi1, '_p_0.0_ss_', sig_strength, sep ="")
} else{
  mypath = paste('Data/M_', M, '_m_', m, '/r_', r, '_pi1_', pi1, '_p_', p, '_ss_', sig_strength, sep ="")
}


files = list.files(path=mypath, full.names = TRUE) 
if (length(files) == total_iterations){
  df <- files %>% 
    lapply(read_csv) %>% 
    bind_rows 
  fdr_power = unlist(apply(df, 2, mean, na.rm = TRUE))
  total = length(fdr_power)
  power_iters = total_iterations * power_iters
  errors = unlist(apply(df, 2, sd, na.rm = TRUE)/(power_iters ** 0.5))
  methods =  c('cPCH-Fisher-BH', 'cPCH-Simes-BH', 'cPCH-Bonferroni-BH',
               'cPCH-Fisher-Storey', 'cPCH-Simes-Storey', 'cPCH-Bonferroni-Storey',
               'Fisher-BH', 'Simes-BH', 'Bonferroni-BH',
               'Fisher-Storey', 'Simes-Storey', 'Bonferroni-Storey',
               'Fisher-DHH', 'Simes-DHH', 'Bonferroni-DHH',
               'Fisher-sDHH', 'Simes-sDHH', 'Bonferroni-sDHH', 'adaFilter')
  fdr_power_df = data.frame(M=M, m=m, r=r, pi1 = pi1, p=p, ss = sig_strength,
                            method = rep(methods, 2), 
                            type = c(rep('FDR', length(methods)), rep('Power', length(methods))), 
                            p_reject = fdr_power[7:total], 
                            se = errors[7:total])
  
  #change data_path to where you want data to be saved
  #R pastes float of form 1.0 as just 1 and the data generated from the cluster uses 1.0
  if (p == 1){
    data_path = paste('Data/M_', M, '_m_', m, '/final_data/r_', r, '_pi1_',pi1, '_p_1.0_ss_', sig_strength, sep ="")
  } else if (p == 0){
    data_path = paste('Data/M_', M, '_m_', m, '/final_data/r_', r, '_pi1_',pi1, '_p_0.0_ss_', sig_strength, sep ="")
  } else{
    data_path = paste('Data/M_', M, '_m_', m, '/final_data/r_', r, '_pi1_',pi1, '_p_', p, '_ss_', sig_strength, sep ="")
  }
  
  
  write.csv(fdr_power_df, file = paste0(data_path, '.csv'))
  
} else {
  print('waiting for full iterations')
}



