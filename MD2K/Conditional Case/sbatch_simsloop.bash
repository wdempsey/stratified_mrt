#!/bin/bash
for BARBETA in 0.020 0.025 0.030; do
    #
    echo "${BARBETA}"
    export BARBETA
    #
    sbatch -o MCsim_${BARBETA}.stdout.txt \
	   -e MCsim_${BARBETA}.stdout.txt \
           --job-name=my_analysis_${BARBETA} \
	   simsloop.sbatch
    #
    sleep 1 # pause to be kind to the scheduler
done



