#!/bin/bash

#SBATCH --job-name=pi_ext_BS.slurm
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=shoumisarkar@ufl.edu     # Where to send mail    
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16	       # Allocate 16 cores to the task
#SBATCH --mem=30GB                     # Job memory request
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --account=qiu
#SBATCH --qos=qiu-b
#SBATCH --output=pi_ext_BS_R_%j_%A_%a.log   # Standard output and error log

echo "Running job for b=$b, with age=$age and corstr=$corstr"
bash 03_run_job_presence_BS.sh
