#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH -p murphy
#SBATCH --mem-per-cpu=4000
#SBATCH -t 0-6:00 # Running time of 2 hours

# --ntasks-per-node will be used in doParallel.R to specify the number of 
# cores to use on the machine. 
# The murphy lab partition has 32 cores per node.
# The general partition has 64 cores per node.

cd $SLURM_SUBMIT_DIR

echo $SLURM_SUBMIT_DIR

# Source the modules
source new-modules.sh

# Load required software modules

# Simple R setup
module load R/3.3.3-fasrc01

# Run program
Rscript estimation.R $BARBETA

