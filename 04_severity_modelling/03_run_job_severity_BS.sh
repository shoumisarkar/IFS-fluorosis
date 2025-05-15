#!/bin/sh

sleepsecs=$[ ( $RANDOM % 100 ) + 100 ]s
sleep $sleepsecs

cd /blue/somnath.datta/shoumisarkar/Fluorosis/Codes/04_severity_modelling

pwd; hostname; date

echo "b=$b, age=$age, corstr=$corstr"  # Debugging line to check the value of age

echo " This script runs at: \n"
date
echo "\n\n"

module load R/4

Rscript fit_severity_BS.R $b $age $corstr
