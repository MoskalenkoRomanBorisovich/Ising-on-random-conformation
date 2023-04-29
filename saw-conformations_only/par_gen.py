import numpy as np

fname = 'parameters/par_conf_'

L = [250]
N = len(L)

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



N_rep = 10000
L_size = [l*5 for l in L]
U = np.linspace(1.0, 1.0, N)

print('U = ', U)
print('L = ', L)

num_file = open('parameters/num_of_files.txt', 'tw')
num_file.write(str(N * N_rep) + '\n')
num_file.write(str(N_rep) + ' # number of replicas per set of parametrs')
num_file.close()

for i in range(N * N_rep):
    origin = open('par_1', 'r')
    dest = open(fname + str(i), 'w')
    for j in range(2):
        dest.write(origin.readline())

    for j in range(2):
        origin.readline()

    t = i // N_rep
    dest.write(str(L_size[t])+'\n')
    dest.write(str(L[t])+'\n')
    dest.write(str(U[t])+'d0\n')
    origin.readline()
    for j in range(11):
        dest.write(origin.readline())

    dest.write(str(2000 + i) + ' ' + str(2000 + i) + '\n')

    origin.close()
    dest.close()
