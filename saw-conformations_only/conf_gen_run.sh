#!/bin/bash
#SBATCH --job-name=confgen
#SBATCH --output=/home/rbmoskalenko/saw-conformations_only/slurm_logs/confgen-%a.slurm.log
#SBATCH --array=0-99

cd parameters

idx=$SLURM_ARRAY_TASK_ID
for i in $(seq $(($idx * 100)) 1 $((($idx+1)*100-1))); do
	../saw_cnf _conf_$i
done