#!/bin/bash

ages=(9 13 17 23) #(9 13 17 23)
corstrs=("independence") #("ar1" "independence" "exchangeable" "jackknifed")
#missing_b=(12 66 77 78 81 82 85 86)

for b in {301..500} #for b in "${missing_b[@]}"
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

