#!/bin/bash

#SBATCH --job-name=pjj50_MC.slurm
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=name@email.com     # Where to send mail    
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16	       # Allocate 16 cores to the task
#SBATCH --mem=30GB                     # Job memory request
#SBATCH --time=00:30:00               # Time limit hrs:min:sec
#SBATCH --account=qos_name
#SBATCH --qos=qos_name
#SBATCH --output=pjj50_MC_R_%j_%A_%a.log   # Standard output and error log

echo "Running job for mc_seed=$mc_seed, with age=$age"
bash 03_run_job_presence_MC.sh
