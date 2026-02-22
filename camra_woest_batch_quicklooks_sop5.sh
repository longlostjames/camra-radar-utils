#!/bin/bash 
#SBATCH --partition=short-serial 
#SBATCH --job-name=camra_quicklooks
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=20:00:00
#SBATCH --array=13
#SBATCH --mem=56G

source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate cao_3_11

time /home/users/cjwalden/git/camra-radar-utils/make_woest_sop_quicklooks.py -d 202307${SLURM_ARRAY_TASK_ID}



