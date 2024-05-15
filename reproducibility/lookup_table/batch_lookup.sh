#!/bin/bash
#SBATCH -n 32 # Number of cores set to 1 since we aren't parallelizing over the 32 cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --time=12:00:00
#SBATCH --partition=janson,janson_bigmem,serial_requeue # Partition to submit to
#SBATCH --mem=20000 # Memory pool for all cores (see also --mem-per-cpu) can lower memory
#SBATCH -o output/hostname_%j.out # File to which STDOUT will be written
#SBATCH -e err/hostname_%j.err # File to which STDERR will be written

module load python/3.10.9-fasrc01
module load R/4.2.2-fasrc01
source activate gurobi_env
python generate_lookup.py ${m} ${r} ${alpha} ${N}
