#cython: language_level=3

import numpy as np

cimport cython
from libc.math cimport exp, tanh
from libcpp.vector cimport vector
from libcpp cimport bool

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double energy(long[::1] spins, 
                   long[:, ::1] neighbors,
                   double J=1.0):
    """Ising model energy of a spin state.
    """
    cdef:
        double ene = 0.0
        Py_ssize_t site, site1, num_neighb

    for site in range(spins.shape[0]):
        num_neighb = neighbors[site, 0]
        for j in range(1, num_neighb+1):
            site1 = neighbors[site, j]
            ene += -J * spins[site] * spins[site1]
    
    # each bond is counted twice, hence divide by two
    return ene / 2.0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double magnetization(long[::1] spins):
    cdef:
        double mag = 0.0
        Py_ssize_t site
    
    for site in range(spins.shape[0]):
        mag += spins[site]

    return mag

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int next_spins(long[::1] spins):
    cdef:
        int coef = -1
        Py_ssize_t site
    
    for site in range(spins.shape[0]):
        spins[site] *= coef
        coef *= spins[site]
        if coef == 1:
            break
    
    return coef

@cython.boundscheck(False)
@cython.wraparound(False)
def calculate(Py_ssize_t L,
             long[:, ::1] neighbors,
             double beta,
             int verbose = 0):
    
    cdef:   
        double T = 1./beta
        double ene
        double ene_total = 0.0
        double mag
        double mag2 = 0.0
        double mag4 = 0.0
        double Z = 0.0
        double p
        long cnt = 0

    # print(np.asarray(neighbors))
    if verbose >= 1:
        print("beta = ", beta, "  T = ", 1./beta)

    
    # initialize spins
    cdef long[::1] spins =  np.ones(L, dtype=int)
    
    while True:
        ene = energy(spins, neighbors)
        p = exp(-beta * ene)
        Z += p
        ene_total += ene * p
        if next_spins(spins) == -1:
            break
        
    ene_total = ene_total / Z
    
    
    return ene_total / L
    
    