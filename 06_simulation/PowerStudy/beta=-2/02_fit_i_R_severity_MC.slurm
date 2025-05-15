#!/bin/bash

#SBATCH --job-name=s30_m2_MC.slurm
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=shoumisarkar@ufl.edu     # Where to send mail    
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16	       # Allocate 16 cores to the task
#SBATCH --mem=30GB                     # Job memory request
#SBATCH --time=12:00:00               # Time limit hrs:min:sec
#SBATCH --account=somnath.datta
#SBATCH --qos=somnath.datta-b
#SBATCH --output=s30_m2_MC_R_%j_%A_%a.log   # Standard output and error log

echo "Running job for mc_seed=$mc_seed, with age=$age"
bash 03_run_job_severity_MC.sh
