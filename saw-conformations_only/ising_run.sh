#!/bin/bash
#SBATCH --job-name=ising
#SBATCH --output=ising-%a.slurm.log
#SBATCH --array=0-49

module purge
module load Python/Anaconda_v10.2019
source activate main_env

idx=$SLURM_ARRAY_TASK_ID
python cy_ising_for_cluster.py -i struct_conf_$idx.dat -o Data/conf_data_$idx.npz