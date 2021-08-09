# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 20:42:05 2021

@author: Роман
"""

import numpy as np
import cy_ising_cluster

import sys, getopt

from mc_lib.observable import RealObservable


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


def mean_conections(neighbors):
    sum = 0
    for line in neighbors:
        sum += line[0]
        
    return sum / len(neighbors)


def measure_structure(struct_1, b_min=0, b_max=1, N=10, n_sweeps=100000, num_therm=100000):
    #draw_conformation(struct_1)
    print("conformation length: ", len(struct_1))
    
    ene = np.empty(N, dtype=RealObservable)
    mag2 = np.empty(N, dtype=RealObservable)
    mag4 = np.empty(N, dtype=RealObservable)
    ene_arr = [0] * N
    neighbors = tabulate_neighbors(struct_1)
    
    # mean_nghbrs = mean_conections(neighbors)
    L = len(struct_1)
    betas = np.linspace(b_min+(b_max-b_min)/N, b_max, N)
    
    for i in range(N):
        #upd_per_sweep = max(1, int(10 * betas[i]))
        ene[i], mag2[i], mag4[i], ene_arr[i] = cy_ising_cluster.simulate(L=L,
                                                                 neighbors=neighbors,
                                                                 beta=betas[i],
                                                                 num_sweeps=n_sweeps,
                                                                 num_therm=num_therm,
                                                                 sampl_frequency = 100000,
                                                                 do_intermediate_measure = 1)
    
    
    
    return ene, mag2, mag4, np.array(ene_arr)


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"i:o:",["ifile=","ofile="])
    
    except getopt.GetoptError:
        print('cy_ising_for_cluster.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-i':
            inputfile = arg
        elif opt == '-o':
            outputfile = arg
    
    struct = read_conformation(inputfile)
    # L = len(struct)
    
    # therm = 300000 is usually enough for L <= 1000
    N = 10
    b_max = 1
    b_min = 0
    
    ene, mag2, mag4, ene_arr = measure_structure(struct, b_min, b_max, N,
                                                        n_sweeps=2000000,
                                                        num_therm=500000)
    betas = np.linspace(b_min+(b_max-b_min)/N, b_max, N)
    print('number of energy measurments: ', ene_arr.shape)
    np.savez(outputfile, ene=ene, mag2=mag2, mag4=mag4, ene_arr=ene_arr, betas=betas)
    

if __name__ == "__main__":
    main(sys.argv[1:])
