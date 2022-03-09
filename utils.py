import numpy as np
import matplotlib.pyplot as plt
from mc_lib.observable import RealObservable
from numba import njit
import h5py

def draw_conformation(struct, spins=None, bridges=None):
    if spins is not None and bridges is not None:
        raise Exception('Cant draw spins and bridges at the same time')
    plt.plot(struct[:, 0], struct[:, 1], '-g', color='gray')
    if spins is not None:
        plt.scatter(struct[spins==1, 0], struct[spins==1, 1], color='red')
        plt.scatter(struct[spins==-1, 0], struct[spins==-1, 1], color='blue')
    elif bridges is not None:
        plt.scatter(struct[np.logical_not(bridges), 0], struct[np.logical_not(bridges), 1], color='purple')
        plt.scatter(struct[bridges, 0], struct[bridges, 1], color='orange')
    else:
        plt.scatter(struct[:, 0], struct[:, 1], color='green')
        
    plt.axis('off')
    plt.show()

def read_conformation(fname):
    struct_conf = []
    
    f = open(fname, 'r')
    f.readline()
    for line in f:
        line_ar = line.split()
        struct_conf += [[int(line_ar[0]), int(line_ar[1])]]
        
    f.close()
    return np.array(struct_conf, dtype=int)

@njit(cache=True)
def tabulate_neighbors(struct): 
    neighb = np.zeros((struct.shape[0], 5), dtype=np.int32)
    for site in range(struct.shape[0]):
        c = struct[site]
        C_arr = np.empty((4, 2), dtype=np.int32)
        C_arr[0, 0] = c[0]
        C_arr[0, 1] = c[1] + 1
        
        C_arr[1, 0] = c[0]
        C_arr[1, 1] = c[1] - 1
        
        C_arr[2, 0] = c[0] + 1
        C_arr[2, 1] = c[1]
        
        C_arr[3, 0] = c[0] - 1
        C_arr[3, 1] = c[1]
        for j in range(C_arr.shape[0]):
            for i in range(struct.shape[0]):
                if struct[i, 0] == C_arr[j, 0] and struct[i, 1] == C_arr[j, 1]:
                    neighb[site, 0] += 1
                    neighb[site, neighb[site, 0]] = i
                    
    return neighb

@njit(cache=True)
def radius_of_gyration(struct):
    x = 0
    y = 0
    for i in range(struct.shape[0]):
        x += struct[i, 0]
        y += struct[i, 1]
    x = x / struct.shape[0]
    y = y / struct.shape[0]
    
    r2_sum = 0
    for i in range(struct.shape[0]):
        r2_sum += (struct[i, 0]-x)**2 + (struct[i, 1]-y)**2
    return np.sqrt(r2_sum / struct.shape[0])

def mean_conections(neighbors):
    sum = 0
    for line in neighbors:
        sum += line[0]
    return sum / len(neighbors)

def generate_1D(L):
    struct = []
    for i in range(L):
        struct += [[i, 0]]
    
    return np.array(struct, dtype=int)

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
    
    return np.array(struct, dtype=int)

def square_2D(side_L):
    struct = [[0, 0]] * side_L**2
    for i in range(side_L):
        for j in range(side_L):
            if(i%2 == 0):
                struct[i*side_L+j] = [i, j]
            
            else:
                struct[i*side_L+j] = [i, side_L - 1 - j]
    
    return np.array(struct, dtype=int)

def R_to_norm(R, L):
    return R / np.sqrt(L)

def norm_to_R(N, L):
    return N * np.sqrt(L)

def ising_1D_true_value(beta, L):
    th = np.tanh(beta)
    return  -th * (1 + th**(L-2)) / (1 + th**L)

class Conformation():
    def __init__(self):
        self.struct = np.empty((0, 2), dtype=int)
        self.ene = np.empty(0, dtype=RealObservable)
        self.mag2 = np.empty(0, dtype=RealObservable)
        self.mag4 = np.empty(0, dtype=RealObservable)
        self.U = np.zeros(0)
        self.betas = np.zeros(0)
        self.R = 0
        self.nBetas = 0
        
    def set_betas(self, _betas):
        self.betas = _betas.copy()
        self.nBetas = self.betas.shape[0]
        
    def load_data(self, data_f_name):
        data = np.load(data_f_name, allow_pickle=True)
        self.ene = data['ene']
        self.mag2 = data['mag2']
        self.mag4 = data['mag4']
        self.nBetas = self.ene.shape[0]
        self.U = np.zeros(self.nBetas)
        for i in range(self.nBetas):
            self.U[i] = 1 - self.mag4[i].mean / (3 * self.mag2[i].mean ** 2)
            
        try:
            self.betas = data['betas']
        except:
            pass
    
    def load_struct(self, struct_f_name):
        self.struct = read_conformation(struct_f_name)
        self.R = radius_of_gyration(self.struct)
        self.L = self.struct.shape[0]
        self.R_norm = R_to_norm(self.R, self.L)
        
#     TODO:
#     def load_hdf5(self, fname):
#         with h5py.File(fname, "r") as f:
#             self.struct = f['structure'][()]
#             self.ene = f['ene'][()]
#             self.mag2 = f['mag2'][()]
#             self.mag4 = f['mag4'][()]
#             self.betas = f['betas'][()]
#         self.nBetas = self.ene.shape[0]
        
def load_Conformations_from_dir(dir_name, load_struct=True, load_data=True):
    if dir_name[-1] != '/' and dir_name[-1] != '\\':
        dir_name += '/'
        
    num_file = open(dir_name + 'num_of_files.txt', 'r')
    N_conf = int(num_file.readline())
    num_file.close()
    
    ret_arr = []
    data_f = dir_name + 'conf_data_'
    struct_f = dir_name + 'struct_conf_'
    for i in range(N_conf):
        cur_conf = Conformation()
        if load_data:
            cur_conf.load_data(data_f + str(i) + '.npz')
        if load_struct:
            cur_conf.load_struct(struct_f + str(i) + '.dat')
        ret_arr.append(cur_conf)
    
    return ret_arr

def save_conformation(struct, fname):
    f = open(fname, 'w')
    f.write('\n');
    for c in struct:
        f.write(str(c[0]) + ' ' + str(c[1]) + '\n')
        
    f.close()

if __name__ == "__main__":
    dir_name = 'Conformations/L250_beta0.1_1_10/'
    conformations = load_Conformations_from_dir(dir_name)
    assert(len(conformations) == 100)
    print(conformations[0].ene[0].mean)
    struct1 = square_1D(3)
    test_file = 'Conformations/Tests/save_conf_test.txt'
    print('Structure before saveing:\n', struct1)
    save_conformation(struct1, test_file)
    struct2 = read_conformation(test_file)
    print('Structure after loading:\n', struct2)
    np.testing.assert_array_equal(struct1, struct2)
    print('raduis of gyration:', radius_of_gyration(struct1))
    
    struct = square_2D(3)
    print(struct)
    tbn = tabulate_neighbors(struct)
    print(tbn)
    