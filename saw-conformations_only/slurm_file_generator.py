# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 19:04:01 2021

@author: Роман
"""

fname = 'script.sh'
par_fname = r'_conf_'


num_file = open('parameters/num_of_files.txt', 'tr')
N = int(num_file.readline())
num_file.close()

file = open(fname, 'tw', newline='\n')

file.write(r'#!/bin/bash'+'\n')
file.write(r'#SBATCH --job-name=conf_ising'+'\n')
file.write(r'#SBATCH --output=gen_test.slurm.log'+'\n')
file.write(r'#SBATCH --ntasks='+str(N)+'\n')
file.write(r'#SBATCH --cpus-per-task=1'+'\n')
file.write('\n')
file.write(r'module load Python/Anaconda_v10.2019'+'\n')
file.write(r'source activate main_env'+'\n')
file.write('\n')


for i in range(N):
    file.write('srun -n 1 '+'./saw_cnf '+par_fname+str(i)+' &\n')
    
file.write('\nwait\n\n')    

for i in range(N):
    file.write('srun -n 1 '+'python cy_ising_for_cluster.py -i struct_conf_'+str(i)+r'.dat -o Data/conf_data_'+str(i)+'.npz &\n')

file.write('\nwait\n\n') 

file.close()

