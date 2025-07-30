#!/bin/bash

ages=(9 13 17 23)

for mc_seed in {1..100}
do
	for age in "${ages[@]}"
	do
		    echo "Submitting job for age=$age and mc_seed=$mc_seed"
		    sbatch --export=ALL,mc_seed=$mc_seed,age=$age 02_fit_i_R_severity_MC.slurm
	done
done
