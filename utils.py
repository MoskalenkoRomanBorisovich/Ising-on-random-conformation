import numpy as np

def read_conformation(fname):
    struct_conf = []
    
    f = open(fname, 'r')
    f.readline()
    for line in f:
        line_ar = line.split()
        struct_conf += [[int(line_ar[0]), int(line_ar[1])]]
        
    f.close()
    return struct_conf

def tabulate_neighbors(struct):
    neighb = np.zeros((len(struct), 5), dtype=int)
    for site in range(len(struct)):
        coordinate = struct[site]
        c1 = [coordinate[0] + 1, coordinate[1]]
        c2 = [coordinate[0] - 1, coordinate[1]]
        c3 = [coordinate[0], coordinate[1] + 1]
        c4 = [coordinate[0], coordinate[1] - 1]
        C_arr = [c1, c2, c3, c4]
        for coord in C_arr:
            try:
                site1 = struct.index(coord)
                neighb[site, 0] += 1
                neighb[site, neighb[site, 0]] = site1
            except:
                continue
    return neighb

def radius_of_gyration(struct):
    x = 0
    y = 0
    for coord in struct:
        x += coord[0]
        y += coord[1]
    x = x / len(struct)
    y = y / len(struct)
    # print(x, y)
    r2_sum = 0
    for coord in struct:
        r2_sum += (coord[0]-x)**2 + (coord[1]-y)**2
    return np.sqrt(r2_sum / len(struct))

def mean_conections(neighbors):
    sum = 0
    for line in neighbors:
        sum += line[0]
    return sum / len(neighbors)

def generate_1D(L):
    struct = []
    for i in range(L):
        struct += [[i, 0]]
    
    return struct

def square_1D(side_L):
    struct = []
    for i in range(side_L):
        struct += [[i, 0]]
        
    for i in range(1, side_L):
        struct += [[side_L-1, i]]
    
    for i in range(side_L - 2, -1, -1):
        struct += [[i, side_L-1]]
        
    for i in range(side_L-2, 0, -1):
        struct += [[0, i]]
    
    return struct

def square_2D(side_L):
    struct = [[0, 0]] * side_L**2
    for i in range(side_L):
        for j in range(side_L):
            if(i%2 == 0):
                struct[i*side_L+j] = [i, j]
            
            else:
                struct[i*side_L+j] = [i, side_L - 1 - j]
    
    return struct

def R_to_norm(R, L):
    return R / np.sqrt(L)

def norm_to_R(N, L):
    return N * np.sqrt(L)

def ising_1D_true_value(beta, L):
    th = np.tanh(beta)
    return  -th * (1 + th**(L-2)) / (1 + th**L)