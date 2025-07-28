#!/bin/bash

ages=(9 13 17 23) 
corstrs=("ar1" "independence" "exchangeable" "jackknifed")

for b in {1..500} 
do
	for age in "${ages[@]}"
	do
		for corstr in "${corstrs[@]}"
		do
		    echo "Submitting job for age=$age and corstr=$corstr, b=$b"
		    sbatch --export=ALL,b=$b,age=$age,corstr=$corstr 02_fit_i_R_presence_BS.slurm
		done
	done
done
