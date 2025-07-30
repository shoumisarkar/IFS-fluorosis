#!/bin/bash

ages=(23) #(9 13 17 23) #(9 13 17 23) #(9 13 17 23)

#missing_b=(12 66 77 78 81 82 85 86)

for mc_seed in {269..500} #for b in "${missing_b[@]}"
do
	for age in "${ages[@]}"
	do
		    echo "Submitting job for age=$age and mc_seed=$mc_seed"
		    sbatch --export=ALL,mc_seed=$mc_seed,age=$age 02_fit_i_R_severity_MC.slurm
	done
done
