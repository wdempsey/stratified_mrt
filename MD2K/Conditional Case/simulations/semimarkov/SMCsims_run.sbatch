#!/bin/bash
#SBATCH -p murphy
#SBATCH -n 20
#SBATCH -N 1                       # All cores on one machine
#SBATCH --mem-per-cpu=4000
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH --continguous # Ensure that all of the cores are on the same Infiniband network
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=dempsey.walter@gmail.com # Email to which notifications will be sent

cd $SLURM_SUBMIT_DIR

# Source the modules
source new-modules.sh

# Load required software modules

# Simle R setup
module load R/3.2.0-fasrc01

# Rmpi setup
# module load intel/15.0.0-fasrc01
# module load mvapich2/2.0-fasrc03
# module load R/3.2.2-fasrc02

# Run program
# R CMD BATCH --no-save --no-restore '--args barbeta=b day=d' semi_treatment_args.R
srun Rscript semi_estimation.R ${BARBETA} > output_${BARBETA}.Rout
