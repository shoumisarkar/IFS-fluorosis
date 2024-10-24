#!/bin/bash

ages=(9) #(9 13 17 23)
corstrs_pres=("jackknifed") #("ar1" "independence" "exchangeable")
corstrs_sev=("jackknifed")
#missing_b=(4 5 7 15 18 29 34 35 36 40 41 44 66 68 69 70 81 82 83 84 87 96)

for b in {81..120} #"${missing_b[@]}" #{101..120}
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