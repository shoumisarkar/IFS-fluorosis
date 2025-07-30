#!/bin/sh

sleepsecs=$[ ( $RANDOM % 100 ) + 100 ]s
sleep $sleepsecs

cd /path/to/Fluorosis/Codes/06_simulation/N=50

pwd; hostname; date

echo "mc_seed=$mc_seed, age=$age"  # Debugging line to check the value of age

echo " This script runs at: \n"
date
echo "\n\n"

module load R/4

Rscript 01c_sim_dat_MC_N_50_combined.R $mc_seed $age
