#!/bin/bash 
#SBATCH --account=woest
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --job-name=camra_array
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --array=22
#SBATCH --mem=64G

source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate cao_3_11

time /home/users/cjwalden/git/camra-radar-utils/proc_camra2ncas_woest_other.py -d 202306${SLURM_ARRAY_TASK_ID}





