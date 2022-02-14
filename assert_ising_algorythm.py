import numpy as np

import cy_ising_cluster
from utils import read_conformation
from utils import square_1D
from utils import square_2D
from utils import tabulate_neighbors
from utils import ising_1D_true_value
import exact_ising

from mc_lib.rndm import RndmWrapper
from mc_lib.observable import RealObservable


class Assert_claster_upd:
    '''
    betas - inverse temperature values at witch tsts will be done
    N - number of values of beta
    '''
    
    def __init__(self,
                 betas_v=np.linspace(0.1, 1, 10),
                 cluste_upd_prob_v = [0.0, 0.5, 1.0]):
        self.betas = betas_v
        self.nBetas = len(self.betas)
        self.cluster_upd_prob = cluste_upd_prob_v.copy()
        self.nProb = len(self.cluster_upd_prob)
        
    def cmp_to_true_1D_values(self, side_L = 10, n_sweeps = 1000000, n_therm=1000000):
        print('Comparison to true 1D values')
        print('betas:', self.betas)
        conf = square_1D(side_L)
        L = len(conf) # length of the conformation = (L - 1) * 4
        neighbors = tabulate_neighbors(conf)
        ene = np.empty(self.nBetas, dtype=RealObservable)
        for i in range(self.nBetas):
            ene[i], _, _, _ = cy_ising_cluster.simulate(neighbors = neighbors,
                                                        beta = self.betas[i],
                                                        num_sweeps = n_sweeps,
                                                        num_therm = n_therm)
        
        true_ene = ising_1D_true_value(self.betas, L)
        ans = True
        for i in range(self.nBetas):
            if abs(ene[i].mean - true_ene[i]) > ene[i].errorbar:
                ans = False
                print('Cluster update results do not match the true values')
                print('beta =', self.betas[i])
                print('algorythm result ene =', ene[i].mean, '+-', ene[i].errorbar)
                print('true value ene =', true_ene[i])
        if ans:
            print('everything is fine')
        return ans
    
    def cmp_to_1_spin_upd(self, conf, n_sweeps = 1000000, n_therm = 1000000):
        print('Compare results with different probabilities of cluster and 1 spin update')
        L = len(conf)
        neighbors = tabulate_neighbors(conf)
        ene = np.empty((self.nProb, self.nBetas), dtype=RealObservable)
        
        for i in range(self.nProb):
            for j in range(self.nBetas):
                ene[i, j], _, _, _ = cy_ising_cluster.simulate(neighbors = neighbors,
                                                                              beta = self.betas[j],
                                                                              num_sweeps = n_sweeps,
                                                                              num_therm = n_therm,
                                                                              cluster_upd_prob = self.cluster_upd_prob[i])
        
        ans = True
        for i1 in range(self.nProb):
            for i2 in range(i1 + 1, self.nProb):
                for j in range(self.nBetas):
                    if abs(ene[i1, j].mean - ene[i2, j].mean) > ene[i1, j].errorbar + ene[i2, j].errorbar:
                        ans = False
                        print('Results does not match')
                        print('beta =', self.betas[i])
                        print('cluster update probability 1 =', self.cluster_upd_prob[i1])
                        print('ene1 =', ene[i1, j].mean, '+-', ene[i1, j].errorbar)
                        print('cluster update probability 2 =', self.cluster_upd_prob[i2])
                        print('ene2 =', ene[i2, j].mean, '+-', ene[i2, j].errorbar)
                        
        if ans:
            print('everything is fine')
            
        return ans
                    
    def cmp_to_exact_solution(self, conf, n_sweeps = 1000000, n_therm = 1000000):
        print('Compare reults algorithm to exact solution of iseng model')
        L = len(conf)
        neighbors = tabulate_neighbors(conf)
        ene = np.empty(self.nBetas, dtype=RealObservable)
        ene_exact = np.empty(self.nBetas, dtype=float)
        
        for i in range(self.nBetas):
            ene[i], _, _, _ = cy_ising_cluster.simulate(neighbors = neighbors,
                                                        beta = self.betas[i],
                                                        num_sweeps = n_sweeps,
                                                        num_therm = n_therm)
        
        for i in range(self.nBetas):
            ene_exact[i] = exact_ising.calculate(L, neighbors, self.betas[i])
            
        ans = True
        for i in range(self.nBetas):
            if abs(ene[i].mean - ene_exact[i]) > ene[i].errorbar:
                ans = False
                print('Cluster update results do not match the exact values')
                print('beta =', self.betas[i])
                print('algorythm result ene =', ene[i].mean, '+-', ene[i].errorbar)
                print('exact value ene =', ene_exact[i])
        if ans:
            print('everything is fine')
        return ans
        

if __name__ == "__main__":
    cluster_assert = Assert_claster_upd()
    cluster_assert.cmp_to_true_1D_values(10, 1000, 1000)
    print('=======================================================================================')
    cluster_assert.cmp_to_1_spin_upd(square_1D(10), 1000, 1000)
    print('=======================================================================================')
    cluster_assert.cmp_to_exact_solution(square_2D(3), 1000, 1000)
    