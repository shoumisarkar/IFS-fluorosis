#!/bin/bash

#SBATCH --job-name=sev.slurm
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=name@email.com     # Where to send mail    
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16	       # Allocate 16 cores to the task
#SBATCH --mem=30GB                     # Job memory request
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --account=qos_name
#SBATCH --qos=qos_name
#SBATCH --output=sev_R_%j_%A_%a.log   # Standard output and error log

echo "Running job with age=$age and corstr=$corstr"
bash 03_run_job_severity.sh
