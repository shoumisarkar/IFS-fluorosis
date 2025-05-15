#!/bin/sh

sleepsecs=$[ ( $RANDOM % 100 ) + 100 ]s
sleep $sleepsecs

cd /blue/somnath.datta/shoumisarkar/Fluorosis/Codes/06_simulation/N=10000_chunks

pwd; hostname; date

echo "mc_seed=$mc_seed, age=$age, chunk=$chunk "  # Debugging line to check the value of age

echo " This script runs at: \n"
date
echo "\n\n"

module load R/4

Rscript 01c_sim_dat_MC_N_10000_combined_chunks.R $mc_seed $age $chunk