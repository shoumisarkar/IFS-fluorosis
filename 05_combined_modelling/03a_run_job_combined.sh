#!/bin/sh

sleepsecs=$[ ( $RANDOM % 100 ) + 100 ]s
sleep $sleepsecs

cd /path/to/Fluorosis/Codes/05_combined_modelling

pwd; hostname; date

echo "The value of age is $age, corstr_pres is $corstr_pres, corstr_sev is $corstr_sev"  # Debugging line to check the value of age

echo " This script runs at: \n"
date
echo "\n\n"

module load R/4

Rscript fit_combined_whole_data_based.R $age $corstr_pres $corstr_sev
