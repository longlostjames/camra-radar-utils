#!/bin/bash 
#!/bin/bash 
#SBATCH --account=woest
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --job-name=ql_sop0
#SBATCH -o slurm_logs/%A_%a.out
#SBATCH -e slurm_logs/%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --array=11-13,17,18,25
#SBATCH --mem=128G

source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate cao_3_11

time /home/users/cjwalden/git/camra-radar-utils/make_woest_sop_quicklooks.py -d 202307${SLURM_ARRAY_TASK_ID}





