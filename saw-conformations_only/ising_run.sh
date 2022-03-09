#!/bin/bash
#SBATCH --job-name=ising
#SBATCH --output=/home/rbmoskalenko/saw-conformations_only/slurm_logs/ising-%a.slurm.log
#SBATCH --array=0-199

module purge
module load Python/Anaconda_v10.2019
source activate main_env



idx=$SLURM_ARRAY_TASK_ID
for i in $(seq $(($idx * 50)) 1 $((($idx+1)*50-1))); do
	python cy_ising_for_cluster.py -i parameters/struct_conf_$i.dat -o Data/conf_data_$i.npz
done