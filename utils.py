import numpy as np
import matplotlib.pyplot as plt
from mc_lib.observable import RealObservable
from numba import njit
import h5py


def draw_conformation(struct, spins=None, bridges=None, draw_path=True):
    if spins is not None and bridges is not None:
        raise Exception('Cant draw spins and bridges at the same time')
    
    plt.figure()
    if draw_path:
        plt.plot(struct[:, 0], struct[:, 1], '-g', color='gray')
    if spins is not None:
        plt.scatter(struct[spins == 1, 0],
                    struct[spins == 1, 1],
                    color='red')
        plt.scatter(struct[spins == -1, 0],
                    struct[spins == -1, 1],
                    color='blue')
    elif bridges is not None:
        plt.scatter(struct[np.logical_not(bridges), 0],
                    struct[np.logical_not(bridges), 1], color='purple')
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
    neighb = np.zeros((struct.shape[0], 5), dtype=np.int64)
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
            if(i % 2 == 0):
                struct[i*side_L+j] = [i, j]

            else:
                struct[i*side_L+j] = [i, side_L - 1 - j]

    return np.array(struct, dtype=int)


def R_to_norm(R, L):
    return R / np.sqrt(L)


def norm_to_R(N, L):
    return N * np.sqrt(L)


def ising_1D_true_value(beta, L):
    """ exact value of internal energy of 1D Ising model"""
    th = np.tanh(beta)
    return -th * (1 + th**(L-2)) / (1 + th**L)


class Conformation():
    def __init__(self):
        self.struct = np.empty((0, 2), dtype=int)
        self.ene = np.empty(0, dtype=RealObservable)
        self.mag_abs = np.empty(0, dtype=RealObservable)
        self.mag2 = np.empty(0, dtype=RealObservable)
        self.mag4 = np.empty(0, dtype=RealObservable)
        self.U = np.zeros(0, dtype=float)
        self.betas = np.zeros(0, dtype=float)
        self.R = 0
        self.nBetas = 0

    def set_betas(self, _betas):
        self.betas = _betas.copy()
        self.nBetas = self.betas.shape[0]

    def load_data(self, data_f_name):
        data = np.load(data_f_name, allow_pickle=True)
        self.ene = data['ene']
        if 'mag_abs' in list(data.keys()):
            self.mag_abs = data['mag_abs']
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


class Conformation_lite():
    def __init__(self):
        self.struct = np.empty((0, 2), dtype=int)
        self.ene = np.empty(0, dtype=float)
        self.ene_er = np.empty(0, dtype=float)
        self.mag_abs = None
        self.mag_abs_er = None
        self.mag2 = np.empty(0, dtype=float)
        self.mag2_er = np.empty(0, dtype=float)
        self.mag4 = np.empty(0, dtype=float)
        self.mag4_er = np.empty(0, dtype=float)
        self.U = np.zeros(0, dtype=float)
        self.betas = np.zeros(0, dtype=float)
        self.R = 0
        self.nBetas = 0

    def set_betas(self, _betas):
        self.betas = _betas.copy()
        self.nBetas = self.betas.shape[0]

    def load_data(self, data_f_name):
        data = np.load(data_f_name, allow_pickle=True)
        ene = data['ene']
        self.ene = np.array([e.mean for e in ene], dtype=float)
        self.ene_er = np.array([e.errorbar for e in ene], dtype=float)
        
        if 'mag_abs' in list(data.keys()):
            mag_abs = data['mag_abs']
            self.mag_abs = np.array([m.mean for m in mag_abs], dtype=float)
            self.mag_abs_er = np.array([m.errorbar for m in mag_abs], dtype=float)
        
        mag2 = data['mag2']
        self.mag2 = np.array([m.mean for m in mag2], dtype=float)
        self.mag2_er = np.array([m.errorbar for m in mag2], dtype=float)
        mag4 = data['mag4']
        self.mag4 = np.array([m.mean for m in mag4], dtype=float)
        self.mag4_er = np.array([m.errorbar for m in mag4], dtype=float)

        self.nBetas = self.ene.shape[0]
        self.U = np.zeros(self.nBetas)
        for i in range(self.nBetas):
            self.U[i] = 1 - self.mag4[i] / (3 * self.mag2[i] ** 2)

        try:
            self.betas = data['betas']
        except:
            pass

    def load_struct(self, struct_f_name):
        self.struct = read_conformation(struct_f_name)
        self.R = radius_of_gyration(self.struct)
        self.L = self.struct.shape[0]
        self.R_norm = R_to_norm(self.R, self.L)


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


def load_Conformations_lite_from_dir(dir_name: str, load_struct=True, load_data=True, load_count=None):
    if dir_name[-1] != '/' and dir_name[-1] != '\\':
        dir_name += '/'

    num_file = open(dir_name + 'num_of_files.txt', 'r')
    N_conf = int(num_file.readline())
    num_file.close()

    if load_count is not None:
        if N_conf < load_count:
            raise Exception("number of conformations in directory is lower than requested ")
        N_conf = load_count
    ret_arr = []
    data_f = dir_name + 'conf_data_'
    struct_f = dir_name + 'struct_conf_'
    for i in range(N_conf):
        cur_conf = Conformation_lite()
        if load_data:
            cur_conf.load_data(data_f + str(i) + '.npz')
        if load_struct:
            cur_conf.load_struct(struct_f + str(i) + '.dat')
        ret_arr.append(cur_conf)

    return ret_arr


def save_conformation(struct, fname):
    f = open(fname, 'w')
    f.write('\n')
    for c in struct:
        f.write(str(c[0]) + ' ' + str(c[1]) + '\n')

    f.close()


# Bridges and clusters


@njit()
def find_1D_sigments_2(neighbors):
    """
    Marks vertexes witch have less then 3 neighbors

    Parameters
    ----------
    neihbors: np.array[:, :] int
        list of neighbors of vertexes

    Returns
    -------
    bridges: np.array[:] bool
        True - if number of neighbors < 3
    """
    bridges = np.zeros(neighbors.shape[0], np.bool_)
    for i in range(neighbors.shape[0]):
        if neighbors[i, 0] < 3:
            bridges[i] = True

    return bridges


@njit()
def dfs_size_clusters_2(v, neighbors, used, bridges_spins, cb_id, cluster):
    """
    Calculates size of a cluster using dfs algorithm


    Parameters
    ----------
    v: int
        id of curent vertex

    neihbors: np.array[:, :] int
        list of neighbors of vertexes

    used: np.array[:] bool
        True if vertex has already been visited

    bridges_spins: np.array[:] bool
        True - if vertex curently marked as bridge

    cb_id: np.array[:] int
        Ids of clusters to which the vertices belong,
        here it will write new Id for visited vertexes

    cluster: int
        Id of the curent cluster

    Return
    ------
    s: int
        number of points in cluster
    """
    s = 1
    cb_id[v] = cluster
    used[v] = True
    for i in range(1, neighbors[v, 0]+1):
        to = neighbors[v, i]
        if used[to] or bridges_spins[to]:
            continue

        s += dfs_size_clusters_2(to, neighbors, used,
                                 bridges_spins, cb_id, cluster)

    return s


@njit()
def clusters_and_bridges(neighbors):
    """
    Calculates number and size of clusters and bridges (1D sigments)
    Bridge - set of onnected spins with 2 or less neighbors
    Cluster - set of connected spins that are not in bridges

    bridges must connect different clusters

    Parameters
    ----------
    neighbors: np.array(L, :) int
        list of neighbors
        neighbors[i, 0] = number of neghbors

    Returns
    -------
    clusters: np.array(:) int
        sizes of clusters

    bridges: np.array(:) int
        sizes of bridges

    bridges_spins: np.array(L) bool
        True if spin is in a bridge
    """
    bridges_spins = find_1D_sigments_2(neighbors)

    n = neighbors.shape[0]

    used = np.zeros(n, dtype=np.bool_)

    cluster_count = 0
    bridges_count = 0
    clusters = np.zeros(n, dtype=np.int32)
    bridges = np.zeros(n, dtype=np.int32)

    cl_id = np.zeros(n, dtype=np.int32)

    bridge_len = 0
    last_c = -1
    for i in range(n):
        if bridges_spins[i]:
            used[i] = True
            bridge_len += 1
        else:
            if used[i]:
                if bridge_len > 0:
                    if cl_id[i] == last_c:  # bridge connects the same cluster
                        for k in range(i-bridge_len, i):
                            bridges_spins[k] = False
                            cl_id[k] = cl_id[i]
                        clusters[cluster_count-1] += bridge_len
                    else:
                        bridges[bridges_count] = bridge_len
                        bridges_count += 1
            else:
                if bridge_len > 0:
                    bridges[bridges_count] = bridge_len
                    bridges_count += 1

                clusters[cluster_count] = dfs_size_clusters_2(
                    i, neighbors, used, bridges_spins, cl_id, cluster_count)
                cluster_count += 1

            last_c = cl_id[i]
            bridge_len = 0

    if cluster_count > 0:
        # add last bridge

        if bridge_len > 0:
            flg = False
            for i in range(1, neighbors[-1, 0]+1):
                if not bridges_spins[neighbors[-1, i]]:
                    if cl_id[neighbors[-1, i]] != last_c:
                        flg = True

            if flg:
                bridges[bridges_count] = bridge_len
                bridges_count += 1
            else:
                for i in range(1, bridge_len+1):
                    bridges_spins[-i] = False
                    cl_id[-i] = last_c
                clusters[last_c] += bridge_len

        # add first bridge
        flg = False
        for i in range(1, neighbors[0, 0]+1):
            if not bridges_spins[neighbors[0, i]]:
                if cl_id[neighbors[0, i]] != 0:
                    flg = True

        if not flg:
            i = 0
            while bridges_spins[i]:
                bridges_spins[i] = False
                cl_id[i] = 0
                clusters[0] += 1
                i += 1

            if i > 0:
                bridges_count -= 1
                bridges = bridges[1:]
    else:
        if bridge_len > 0:
            bridges[bridges_count] = bridge_len
            bridges_count += 1

    clusters = clusters[:cluster_count]
    bridges = bridges[:bridges_count]

    return clusters, bridges, bridges_spins

def clusters_and_bridges_from_list(conformations):
    """
    Finds clusters and bridges for every conformation in list
    """
    k = 0
    clusters = []
    bridges = []
    for conf in conformations:
        k += 1
        if k % 1000 == 0:
            print(k)
        neighbors = tabulate_neighbors(conf.struct)
        c, b, b_s = clusters_and_bridges(neighbors)
        clusters.append(c)
        bridges.append(b)
        
    return clusters, bridges

def magnetic_susceptibility(conf):
    """ calculates magnetic susceptibility of conformation """
    return conf.betas * (conf.mag2 - conf.mag_abs**2)


def find_ms_peaks(confs):
    p = []
    for c in confs:
        ms = magnetic_susceptibility(c)
        p.append(np.argmax(ms))
        
    return p

def generate_cluster_conformation(W, H, N, L):
    """
    generates clusterizd conformation structure
    example: W=2, H=5, N=2, L=3
        ####   ####
        ####   ####
        ####   ####
        ####   ####
        ###########
    
    Parameters
    ----------
    W: int
        half ol the cluster width
    H int
        hight of cluster
    N: int
        number of clusters
    L: np.array
        length of bridges between clusters
    """
    struct = []
    for x in range(W * 2):
        for y in range(H):
            if x % 2 == 0:
                struct.append([x, y])
            else:
                struct.append([x, H-1-y])
                
    for i in range(N-1):
        offset = i * (W * 2 + L)
        for j in range(L):
            struct.append([offset+W*2+j, 0])
        
        for x in range(offset+W*2+L, offset+W*2+L+W*2):
            for y in range(H):
                if (x - W*2-L - offset) % 2 == 0:
                    struct.append([x, y])
                else:
                    struct.append([x, H-1-y])
        
    return np.array(struct, dtype=int)


def generate_clusters_simple(W, H, N, L, pos, make_bridges=True):
    """
    Generates clusterized graph structure, without requireing it being a conformation

    if make bridges == True, clusters will not be connected 
    Parameters
    ----------
    W : int
        width of cluster
    H : int
        hight of cluster
    N : int
        Number of clusters
    L : int
        lengths of bridges
    P : list[int]
        positions of bridges (hights)
    """
    struct = []
    # add clusters
    offset = 0
    for cl in range(N):
        for i in range(W):
            for j in range(H):
                struct.append([offset+i, j])
        
        offset += W + L
    
    # add bridges
    if make_bridges:
        offset = W
        for br in range(N-1):
            for i in range(L):
                struct.append([offset+i, pos[br]])
            offset += W + L

    return np.array(struct)


def mag_sus1D(b: np.ndarray, N: int, J=1.0):
    """
    Calculates magnetic susceptibility of 1D Ising model with open boundary
    conditions.
    Result is not scaled. Divide by N**2 to scale values.

    Parameters
    ----------
    b : np.array
        beta array
    N : int
        number of spins
    J : float, optional
        magnetic interaction coefficient, by default 1.0

    Returns
    -------
    np.array
        magnetic susceptibility
    """
    e2 = np.exp(2*b*J)
    e4 = e2 ** 2
    X = b/2 * ((2*N*e2 - e4+1) + np.tanh(b*J)**(N-1) *(e4 - 2*e2 + 1))
    return X

def mag_sus1D_2(b: np.ndarray, N: int, J=1.0):
    """_summary_

    Parameters
    ----------
    b : np.array
        beta array
    N : int
        number of spins
    J : float, optional
        magnetic interaction coefficient, by default 1.0

    Returns
    -------
    np.array
        magnetic susceptibility
    """
    tanJb = np.tanh(J*b)
    X = b * (N*(1+2*tanJb/(1-tanJb)) - 2*tanJb*(1-tanJb**N)/(1-tanJb)**2)
    return X

def check_dataset(dataset, betas, err:float, length=None):
    """
    checks if dataset is correct
    checks:
    - number of files
    - beta values
    - real observables convergence

    Parameters
    ----------
    directory : list
        Conformations array
    betas : numpy.array
        expected beta values
    err : float
        observables error threshold value
    length : int(optional)
        length of conformations

    Returns
    -------
    numpy.array
        indexes of conformations with wrong values of beta
    numpy.array
        indexes of conformations with not converged magnetization
    """
    wrong_betas = [] # list of indexes 
    wrong_mag = []
    max_mag_abs = 0
    max_mag_2 = 0
    max_mag_4 = 0
    for i, c in enumerate(dataset):
        if not np.array_equal(c.betas, betas):
            wrong_betas.append(i)
        
        if length is not None:
            l = len(c.struct)
            if length != l:
                raise Exception(f'wong conformations length: expected {length}, got {l} in conformation #{i}')

        for m in c.mag_abs:
            max_mag_abs = max(max_mag_abs, m.errorbar)
            if m.errorbar > err:
                wrong_mag.append(i)
                break
        for m in c.mag2:
            max_mag_2 = max(max_mag_2, m.errorbar)
            if m.errorbar > err:
                wrong_mag.append(i)
                break
        for m in c.mag4:
            max_mag_4 = max(max_mag_4, m.errorbar)
            if m.errorbar > err:
                wrong_mag.append(i)
                break
    
    wrong_betas = np.array(wrong_betas) 
    wrong_mag = np.array(wrong_mag)

    ret = {
        "wrong_betas": wrong_betas,
        "wrong_mag": wrong_mag,
        "max_abs_er": max_mag_abs,
        "max_2_er": max_mag_2,
        "max_4_er": max_mag_4
    }
    return ret
            
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

    conformations = load_Conformations_lite_from_dir(dir_name)
    print(conformations[1].mag2)