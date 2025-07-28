#!/bin/sh

sleepsecs=$[ ( $RANDOM % 100 ) + 100 ]s
sleep $sleepsecs

cd /path/to/Fluorosis/Codes/04_severity_modelling

pwd; hostname; date

echo "The value of age is $age, corstr is $corstr"  # Debugging line to check the value of age

echo " This script runs at: \n"
date
echo "\n\n"

module load R/4

Rscript fit_severity_whole_data_based.R $age $corstr
