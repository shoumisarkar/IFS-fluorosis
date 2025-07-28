#!/bin/bash

ages=(9 13 17 23)
corstrs_pres=("jackknifed") #("ar1" "independence" "exchangeable" "jackknifed")
corstrs_sev=("jackknifed")

for b in {1..500} 
do
	for age in "${ages[@]}"
	do
		for corstr_pres in "${corstrs_pres[@]}"
		do
			for corstr_sev in "${corstrs_sev[@]}"
			do		
 				echo "Submitting job for age=$age, corstr_pres=$corstr_pres, corstr_sev=$corstr_sev, b=$b"
		    		sbatch --export=ALL,b=$b,age=$age,corstr_pres=$corstr_pres,corstr_sev=$corstr_sev 02_fit_i_R_combined_BS.slurm
			done
		done
	done
done
