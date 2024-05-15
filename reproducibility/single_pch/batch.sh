#!/bin/bash
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --time=12:00:00
#SBATCH --partition=     # Partition to submit to
#SBATCH --mem=   # Change to memory usage for your cluster
#SBATCH -o output/hostname_%j.out # Change to desired folder and file to which STDOUT will be written
#SBATCH -e err/hostname_%j.err # Change to desired folder and file to which STDERR will be written

module load Anaconda3/2020.11 #change to load Python and Anaconda3 python packages for your cluster
module load R/4.0.5-fasrc01 #change to load R for your cluster
export R_LIBS_USER=$HOME/apps/R_4.0.5:$R_LIBS_USER #Change for your cluster
python single_pc_sim.py ${m} ${r} ${k} ${sig_strength} ${alpha} ${replicates} ${N}
