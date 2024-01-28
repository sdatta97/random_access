#!/bin/bash

#SBATCH --mem=128GB
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=./logs/%a_%A.out
#SBATCH --error=./logs/%a_%A.err


#SBATCH --job-name=randomaccess
#SBATCH --mail-type=END
#SBATCH --mail-user=sd3927@nyu.edu

echo "$SLURM_ARRAY_TASK_ID"

module purge
module load python/3.8

python3 modified_mosely.py
exit
