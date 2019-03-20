#!/bin/bash
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH -p murphy
#SBATCH --mem=120000
#SBATCH -t 1-00:00              # Runtime in D-HH:MM
# -c will be used in doParallel.R to specify the number of
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
srun -c $SLURM_CPUS_PER_TASK Rscript ti_estimation.R $BARBETA

