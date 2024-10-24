#!/bin/bash

#SBATCH --job-name=pres.slurm
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=shoumisarkar@ufl.edu     # Where to send mail    
#SBATCH --nodes=5
#SBATCH --ntasks=20                    # Run on a single CPU
#SBATCH --mem=30GB                     # Job memory request
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --account=biostat-dept
#SBATCH --qos=biostat-dept-b
#SBATCH --output=pres_R_%j_%A_%a.log   # Standard output and error log

echo "Running job with age=$age and corstr=$corstr"
bash 03_run_job_presence.sh
