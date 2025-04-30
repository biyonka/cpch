#!/bin/bash
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --time=12:00:00
#SBATCH --partition=sapphire,shared,seas_compute,janson,janson_cascade,janson_bigmem,serial_requeue  # Partition to submit to
#SBATCH --mem=16GB # Change to memory usage for your cluster


module load python/3.10.9-fasrc01
module load R/4.2.2-fasrc01
source activate gurobi_env
export R_LIBS_USER=$HOME/apps/R_4.0.5:$R_LIBS_USER

#export R_LIBS_USER=$HOME/apps/R_4.0.5:$R_LIBS_USER #Change for your cluster
python single_pch_sim_dependence.py ${m} ${r} ${k} ${sig_strength} ${alpha} ${rho} ${replicates} ${N}
