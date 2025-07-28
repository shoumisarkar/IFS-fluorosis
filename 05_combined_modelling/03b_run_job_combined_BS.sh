#!/bin/sh

sleepsecs=$[ ( $RANDOM % 100 ) + 100 ]s
sleep $sleepsecs

cd /path/to/Fluorosis/Codes/05_combined_modelling

pwd; hostname; date

echo "b=$b, age=$age, corstr_pres=$corstr_pres, corstr_sev=$corstr_sev"  # Debugging line to check the value of age

echo " This script runs at: \n"
date
echo "\n\n"

module load R/4

Rscript fit_combined_BS.R $b $age $corstr_pres $corstr_sev
