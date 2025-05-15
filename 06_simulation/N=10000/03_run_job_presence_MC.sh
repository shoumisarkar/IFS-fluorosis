#!/bin/sh

sleepsecs=$[ ( $RANDOM % 100 ) + 100 ]s
sleep $sleepsecs

cd /blue/somnath.datta/shoumisarkar/Fluorosis/Codes/06_simulation/N=10000

pwd; hostname; date

echo "mc_seed=$mc_seed, age=$age"  # Debugging line to check the value of age

echo " This script runs at: \n"
date
echo "\n\n"

module load R/4

Rscript 01a_sim_dat_MC_N_10000_presence.R $mc_seed $age