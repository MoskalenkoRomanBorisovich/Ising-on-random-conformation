#!/bin/bash
#SBATCH --job-name=conf_ising
#SBATCH --output=gen_test.slurm.log
#SBATCH --ntasks=100
#SBATCH --cpus-per-task=1

module load Python/Anaconda_v10.2019
source activate main_env

srun -n 1 ./saw_cnf _conf_0 &
srun -n 1 ./saw_cnf _conf_1 &
srun -n 1 ./saw_cnf _conf_2 &
srun -n 1 ./saw_cnf _conf_3 &
srun -n 1 ./saw_cnf _conf_4 &
srun -n 1 ./saw_cnf _conf_5 &
srun -n 1 ./saw_cnf _conf_6 &
srun -n 1 ./saw_cnf _conf_7 &
srun -n 1 ./saw_cnf _conf_8 &
srun -n 1 ./saw_cnf _conf_9 &
srun -n 1 ./saw_cnf _conf_10 &
srun -n 1 ./saw_cnf _conf_11 &
srun -n 1 ./saw_cnf _conf_12 &
srun -n 1 ./saw_cnf _conf_13 &
srun -n 1 ./saw_cnf _conf_14 &
srun -n 1 ./saw_cnf _conf_15 &
srun -n 1 ./saw_cnf _conf_16 &
srun -n 1 ./saw_cnf _conf_17 &
srun -n 1 ./saw_cnf _conf_18 &
srun -n 1 ./saw_cnf _conf_19 &
srun -n 1 ./saw_cnf _conf_20 &
srun -n 1 ./saw_cnf _conf_21 &
srun -n 1 ./saw_cnf _conf_22 &
srun -n 1 ./saw_cnf _conf_23 &
srun -n 1 ./saw_cnf _conf_24 &
srun -n 1 ./saw_cnf _conf_25 &
srun -n 1 ./saw_cnf _conf_26 &
srun -n 1 ./saw_cnf _conf_27 &
srun -n 1 ./saw_cnf _conf_28 &
srun -n 1 ./saw_cnf _conf_29 &
srun -n 1 ./saw_cnf _conf_30 &
srun -n 1 ./saw_cnf _conf_31 &
srun -n 1 ./saw_cnf _conf_32 &
srun -n 1 ./saw_cnf _conf_33 &
srun -n 1 ./saw_cnf _conf_34 &
srun -n 1 ./saw_cnf _conf_35 &
srun -n 1 ./saw_cnf _conf_36 &
srun -n 1 ./saw_cnf _conf_37 &
srun -n 1 ./saw_cnf _conf_38 &
srun -n 1 ./saw_cnf _conf_39 &
srun -n 1 ./saw_cnf _conf_40 &
srun -n 1 ./saw_cnf _conf_41 &
srun -n 1 ./saw_cnf _conf_42 &
srun -n 1 ./saw_cnf _conf_43 &
srun -n 1 ./saw_cnf _conf_44 &
srun -n 1 ./saw_cnf _conf_45 &
srun -n 1 ./saw_cnf _conf_46 &
srun -n 1 ./saw_cnf _conf_47 &
srun -n 1 ./saw_cnf _conf_48 &
srun -n 1 ./saw_cnf _conf_49 &
srun -n 1 ./saw_cnf _conf_50 &
srun -n 1 ./saw_cnf _conf_51 &
srun -n 1 ./saw_cnf _conf_52 &
srun -n 1 ./saw_cnf _conf_53 &
srun -n 1 ./saw_cnf _conf_54 &
srun -n 1 ./saw_cnf _conf_55 &
srun -n 1 ./saw_cnf _conf_56 &
srun -n 1 ./saw_cnf _conf_57 &
srun -n 1 ./saw_cnf _conf_58 &
srun -n 1 ./saw_cnf _conf_59 &
srun -n 1 ./saw_cnf _conf_60 &
srun -n 1 ./saw_cnf _conf_61 &
srun -n 1 ./saw_cnf _conf_62 &
srun -n 1 ./saw_cnf _conf_63 &
srun -n 1 ./saw_cnf _conf_64 &
srun -n 1 ./saw_cnf _conf_65 &
srun -n 1 ./saw_cnf _conf_66 &
srun -n 1 ./saw_cnf _conf_67 &
srun -n 1 ./saw_cnf _conf_68 &
srun -n 1 ./saw_cnf _conf_69 &
srun -n 1 ./saw_cnf _conf_70 &
srun -n 1 ./saw_cnf _conf_71 &
srun -n 1 ./saw_cnf _conf_72 &
srun -n 1 ./saw_cnf _conf_73 &
srun -n 1 ./saw_cnf _conf_74 &
srun -n 1 ./saw_cnf _conf_75 &
srun -n 1 ./saw_cnf _conf_76 &
srun -n 1 ./saw_cnf _conf_77 &
srun -n 1 ./saw_cnf _conf_78 &
srun -n 1 ./saw_cnf _conf_79 &
srun -n 1 ./saw_cnf _conf_80 &
srun -n 1 ./saw_cnf _conf_81 &
srun -n 1 ./saw_cnf _conf_82 &
srun -n 1 ./saw_cnf _conf_83 &
srun -n 1 ./saw_cnf _conf_84 &
srun -n 1 ./saw_cnf _conf_85 &
srun -n 1 ./saw_cnf _conf_86 &
srun -n 1 ./saw_cnf _conf_87 &
srun -n 1 ./saw_cnf _conf_88 &
srun -n 1 ./saw_cnf _conf_89 &
srun -n 1 ./saw_cnf _conf_90 &
srun -n 1 ./saw_cnf _conf_91 &
srun -n 1 ./saw_cnf _conf_92 &
srun -n 1 ./saw_cnf _conf_93 &
srun -n 1 ./saw_cnf _conf_94 &
srun -n 1 ./saw_cnf _conf_95 &
srun -n 1 ./saw_cnf _conf_96 &
srun -n 1 ./saw_cnf _conf_97 &
srun -n 1 ./saw_cnf _conf_98 &
srun -n 1 ./saw_cnf _conf_99 &

wait

srun -n 1 python cy_ising_for_cluster.py -i struct_conf_0.dat -o Data/conf_data_0.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_1.dat -o Data/conf_data_1.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_2.dat -o Data/conf_data_2.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_3.dat -o Data/conf_data_3.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_4.dat -o Data/conf_data_4.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_5.dat -o Data/conf_data_5.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_6.dat -o Data/conf_data_6.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_7.dat -o Data/conf_data_7.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_8.dat -o Data/conf_data_8.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_9.dat -o Data/conf_data_9.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_10.dat -o Data/conf_data_10.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_11.dat -o Data/conf_data_11.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_12.dat -o Data/conf_data_12.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_13.dat -o Data/conf_data_13.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_14.dat -o Data/conf_data_14.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_15.dat -o Data/conf_data_15.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_16.dat -o Data/conf_data_16.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_17.dat -o Data/conf_data_17.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_18.dat -o Data/conf_data_18.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_19.dat -o Data/conf_data_19.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_20.dat -o Data/conf_data_20.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_21.dat -o Data/conf_data_21.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_22.dat -o Data/conf_data_22.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_23.dat -o Data/conf_data_23.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_24.dat -o Data/conf_data_24.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_25.dat -o Data/conf_data_25.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_26.dat -o Data/conf_data_26.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_27.dat -o Data/conf_data_27.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_28.dat -o Data/conf_data_28.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_29.dat -o Data/conf_data_29.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_30.dat -o Data/conf_data_30.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_31.dat -o Data/conf_data_31.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_32.dat -o Data/conf_data_32.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_33.dat -o Data/conf_data_33.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_34.dat -o Data/conf_data_34.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_35.dat -o Data/conf_data_35.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_36.dat -o Data/conf_data_36.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_37.dat -o Data/conf_data_37.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_38.dat -o Data/conf_data_38.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_39.dat -o Data/conf_data_39.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_40.dat -o Data/conf_data_40.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_41.dat -o Data/conf_data_41.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_42.dat -o Data/conf_data_42.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_43.dat -o Data/conf_data_43.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_44.dat -o Data/conf_data_44.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_45.dat -o Data/conf_data_45.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_46.dat -o Data/conf_data_46.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_47.dat -o Data/conf_data_47.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_48.dat -o Data/conf_data_48.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_49.dat -o Data/conf_data_49.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_50.dat -o Data/conf_data_50.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_51.dat -o Data/conf_data_51.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_52.dat -o Data/conf_data_52.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_53.dat -o Data/conf_data_53.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_54.dat -o Data/conf_data_54.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_55.dat -o Data/conf_data_55.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_56.dat -o Data/conf_data_56.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_57.dat -o Data/conf_data_57.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_58.dat -o Data/conf_data_58.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_59.dat -o Data/conf_data_59.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_60.dat -o Data/conf_data_60.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_61.dat -o Data/conf_data_61.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_62.dat -o Data/conf_data_62.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_63.dat -o Data/conf_data_63.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_64.dat -o Data/conf_data_64.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_65.dat -o Data/conf_data_65.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_66.dat -o Data/conf_data_66.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_67.dat -o Data/conf_data_67.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_68.dat -o Data/conf_data_68.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_69.dat -o Data/conf_data_69.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_70.dat -o Data/conf_data_70.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_71.dat -o Data/conf_data_71.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_72.dat -o Data/conf_data_72.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_73.dat -o Data/conf_data_73.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_74.dat -o Data/conf_data_74.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_75.dat -o Data/conf_data_75.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_76.dat -o Data/conf_data_76.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_77.dat -o Data/conf_data_77.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_78.dat -o Data/conf_data_78.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_79.dat -o Data/conf_data_79.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_80.dat -o Data/conf_data_80.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_81.dat -o Data/conf_data_81.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_82.dat -o Data/conf_data_82.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_83.dat -o Data/conf_data_83.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_84.dat -o Data/conf_data_84.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_85.dat -o Data/conf_data_85.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_86.dat -o Data/conf_data_86.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_87.dat -o Data/conf_data_87.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_88.dat -o Data/conf_data_88.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_89.dat -o Data/conf_data_89.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_90.dat -o Data/conf_data_90.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_91.dat -o Data/conf_data_91.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_92.dat -o Data/conf_data_92.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_93.dat -o Data/conf_data_93.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_94.dat -o Data/conf_data_94.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_95.dat -o Data/conf_data_95.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_96.dat -o Data/conf_data_96.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_97.dat -o Data/conf_data_97.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_98.dat -o Data/conf_data_98.npz &
srun -n 1 python cy_ising_for_cluster.py -i struct_conf_99.dat -o Data/conf_data_99.npz &

wait

