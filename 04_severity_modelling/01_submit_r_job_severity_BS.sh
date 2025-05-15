#!/bin/bash

ages=(9 13 17 23) #(9 13 17 23)
corstrs=("independence") #("jackknifed") #("ar1" "independence" "exchangeable")
#missing_b=(3 4 5 7 18 22 23 24 29 35 36 38 40 41 43 44 47 53 61 66 68 69 70 76 78 81 83 87 96)

for b in {121..200} #for b in "${missing_b[@]}"
do
	for age in "${ages[@]}"
	do
		for corstr in "${corstrs[@]}"
		do
		    echo "Submitting job for age=$age and corstr=$corstr, b=$b"
		    sbatch --export=ALL,b=$b,age=$age,corstr=$corstr 02_fit_i_R_severity_BS.slurm
		done
	done
done