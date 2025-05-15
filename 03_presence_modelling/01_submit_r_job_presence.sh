#!/bin/bash

ages=(9 13 17 23)
corstrs=("jackknifed") #("independence" "exchangeable" "ar1" "jackknifed")

for age in "${ages[@]}"
do
	for corstr in "${corstrs[@]}"
	do

	    echo "Submitting job for age=$age and corstr=$corstr"
	    sbatch --export=ALL,age=$age,corstr=$corstr 02_fit_i_R_presence.slurm
	done
done