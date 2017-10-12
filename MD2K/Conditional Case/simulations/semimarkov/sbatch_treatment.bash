#!/bin/bash
for BARBETA in 0.020 0.025 0.030; do
    for DAY in $(seq 2 10); do
	#
	echo "${BARBETA}, ${DAY}"
	export BARBETA DAY
	#
	sbatch -o out_p${BARBETA}_s${DAY}.stdout.txt \
	       -e out_p${BARBETA}_s${DAY}.stdout.txt \
	       --job-name=my_analysis_${BARBETA}_${DAY} \
	       SMC_treatment.sbatch
	#
	sleep 1 # pause to be kind to the scheduler
    done
done



