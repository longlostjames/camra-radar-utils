#!/bin/bash 
#SBATCH --partition=short-serial
#SBATCH --job-name=camra_array
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=15:00:00
#SBATCH --array=11,12,17,25
#SBATCH --mem=24G

source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate cao_3_11

time /home/users/cjwalden/git/camra-radar-utils/proc_camra2ncas_woest_iop.py -d 202307${SLURM_ARRAY_TASK_ID}

