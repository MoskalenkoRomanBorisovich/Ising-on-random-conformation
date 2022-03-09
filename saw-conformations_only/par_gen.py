import numpy as np

fname = 'parameters/par_conf_'

L = [1000]
N = len(L)

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
    dest.write(str(U[t])+'\n')
    origin.readline()
    for j in range(11):
        dest.write(origin.readline())
    
    dest.write(str(2000 + i) + ' ' + str(2000 + i) + '\n')
    
    origin.close()
    dest.close()