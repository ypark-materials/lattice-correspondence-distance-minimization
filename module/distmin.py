'''
distmin.py

Distmin module contains functions for finding the correspondance matrix with the transformation stretch tensor with the minimum distance. 

Author: Yunsu Park
Created: April 9 2025
Affiliation: University of California, Santa Barbara
Contact: yunsu@ucsb.edu
'''

import numpy as np

def calcDist(E_ref, E_def, P_ref, P_def):
    '''
    Calculates the distance function for a given lattice vector (E) and lattice correspondance matrix (P)

    Parameters:
        E_ref (ndarray [shape (3, 3)]): 
            lattice vector for the reference configuration
        E_def (ndarray [shape (3, 3)]): 
            lattice vector for the deformed configuration
        P_ref (ndarray [shape (3, 3)]): 
            lattice correspondance matrix for the reference configuration
        P_def (ndarray [shape (3, 3)]): 
            lattice correspondance matrix for the deformed configuration

    Returns:
        distFunc (float): 
            distance function for the given E and P
        U (ndarray [shape (3, 3)]):
            transformation stretch tensor for the given E and P
    '''
    # Calculates the deformation gradient from the lattice vectors and its correspondance matrix
    F = np.dot(E_def @ P_def, np.linalg.inv(E_ref @ P_ref))

    # Polar Decomposes the deformation gradiet to calculate the stretch tensor U
    C = F.T@F
    eig_val, eig_vec = np.linalg.eig(C)
    lam = np.sqrt(eig_val)
    U = lam[0]*np.outer(eig_vec[0],eig_vec[0]) + lam[1]*np.outer(eig_vec[1],eig_vec[1]) + lam[2]*np.outer(eig_vec[2],eig_vec[2])

    # Calculates the distance function defined by Chen et al. 
    U_inv = np.linalg.inv(U) 
    U_inv_sq = np.dot(U_inv, U_inv) 
    A = U_inv_sq - np.eye(3)

    distFunc = np.trace(A.T @ A)

    return distFunc, U