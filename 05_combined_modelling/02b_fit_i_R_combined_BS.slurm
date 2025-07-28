#!/bin/bash

#SBATCH --job-name=cjj_BS.slurm
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=name@email.com     # Where to send mail    
#SBATCH --nodes=1		       #Use only 1 node
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16	       # Allocate 16 cores to the task
#SBATCH --mem=30GB                     # Job memory request
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --account=qos_name
#SBATCH --qos=qos_name
#SBATCH --output=cjj_ext_BS_R_%j_%A_%a.log   # Standard output and error log

echo "Running job for b=$b, with age=$age, corstr_pres=$corstr_pres, corstr_sev=$corstr_sev"
bash 03_run_job_combined_BS.sh
