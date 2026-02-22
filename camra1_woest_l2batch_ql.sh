#!/bin/bash 
#SBATCH --job-name="camra0_ql_woest"
#SBATCH --time=2:00:00
#SBATCH --mem=128G
#SBATCH --account=woest
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --array=15,18,24,25


source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate cao_3_11

time /home/users/cjwalden/git/camra-radar-utils/make_woest_iop_quicklooks.py -d 202308${SLURM_ARRAY_TASK_ID}



