import numpy as np
import os
# dafault_params = {
#     "dimensionality": 2,
#     "lattice_size": 250,
#     "L": 50,
#     "U" : 1.0,
#     "J": 1.0,
#     "new_conf": 0, # 0 if new configuration, 1 if old one
#     "new_stats": 0, # 0 if new statistics,    1 if old one
#     "initial U": 0.0,
#     "step_print": 5e5,
#     "step_write": 1e7,
#     "step for measuring": 1000,
#     "move_prob": 0.4,
#     "reconnect_prob": 0.2,
#     "time_limit": 0.01,
#     "seed_rng": '4836 2748'
# }


def write_params_to_dir(dir_name, L, U, N):
    num_file = open(dir_name + 'num_of_files.txt', 'tw')
    num_file.write(str(N) + '\n')
    num_file.close()
    if dir_name[-1] != '/':
        fname = dir_name + '/par_conf_'
    else:
        fname = dir_name + 'par_conf_'

    for i in range(N):
        origin = open('par_1', 'r')
        dest = open(fname + str(i), 'w')
        for j in range(2):
            dest.write(origin.readline())

        for j in range(2):
            origin.readline()

        dest.write(str(L*5)+'\n')
        dest.write(str(L)+'\n')
        dest.write(str(U)+'d0\n')
        origin.readline()
        for j in range(11):
            dest.write(origin.readline())

        dest.write(str(2000 + i) + ' ' + str(2000 + i) + '\n')

        origin.close()
        dest.close()


param_dir = "./parameters/"

L = [250, 500, 1000, 2000]

N_rep = 1000
L_size = [l*5 for l in L]
U = [0.2, 0.4, 0.6, 0.8]

print('U = ', U)
print('L = ', L)


for u in U:
    U_dir = param_dir + f"U={u}/"
    if not os.path.exists(U_dir):
        os.makedirs(U_dir, exist_ok=True)
    for l in L:
        dir_name = U_dir + f"L={l}/"
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        write_params_to_dir(dir_name, l, u, N_rep)



