#!/bin/bash
#SBATCH --job-name=congen
#SBATCH --output=congen-%a.slurm.log
#SBATCH --array=0-49

idx=$SLURM_ARRAY_TASK_ID
./saw_cnf _conf_$idx