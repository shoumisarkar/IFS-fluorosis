#!/bin/bash

ages=(9 13 17 23)
corstrs_pres=("independence") #("independence" "exchangeable" "ar1" "jackknifed")
corstrs_sev=("independence") #("independence" "exchangeable" "ar1" "jackknifed")

for age in "${ages[@]}"
do
	for corstr_pres in "${corstrs_pres[@]}"
	do
		for corstr_sev in "${corstrs_sev[@]}"
		do
			echo "Submitting job for age=$age, corstr_pres=$corstr_pres, corstr_sev=$corstr_sev"
	    		sbatch --export=ALL,age=$age,corstr_pres=$corstr_pres,corstr_sev=$corstr_sev 02_fit_i_R_combined.slurm
		done
	done
done
