#!/bin/bash
# SMC_run.sbatch
#
#SBATCH -p murphy # partition
#SBATCH -n 1 # one core
#SBATCH -N 1 # on one machine
#SBATCH -t 2-00:00              # Runtime in D-HH:MM
#SBATCH --mem 1000 # Memory request
#

#SBATCH --mail-user=dempsey.walter@gmail.com
#SBATCH --mail-type=ALL

##   Didn't know if SLURM scripts execute in the home directory

cd $SLURM_SUBMIT_DIR

# Source the modules
source new-modules.sh

# Load required software modules

# Rmpi setup                                                                                                                                                               
module load intel/15.0.0-fasrc01                                                                                                                                          
module load mvapich2/2.0-fasrc03                                                                                                                                          
module load R/3.2.2-fasrc02 

# module load R/3.2.0-fasrc01

# Run program
# R CMD BATCH --no-save --no-restore '--args barbeta=b day=d' semi_treatment_args.R
Rscript semi_treatment.R ${BARBETA} ${DAY} > output_p${BARBETA}_s${DAY}.Rout
