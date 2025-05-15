#!/bin/bash

ages=(9 13 17 23) #(9 13 17 23)

#missing_b=(12 66 77 78 81 82 85 86)

mc_seed=(2)

for chunk in {2..9} #for b in "${missing_b[@]}"
do
	for age in "${ages[@]}"
	do
		    echo "Submitting job for age=$age and mc_seed=$mc_seed , chunk=$chunk "
		    sbatch --export=ALL,mc_seed=$mc_seed,age=$age,chunk=$chunk 02_fit_i_R_combined_MC.slurm
	done
done
