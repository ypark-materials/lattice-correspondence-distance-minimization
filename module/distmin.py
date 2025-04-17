'''
distmin.py

Distmin module contains functions for finding the correspondance matrix with the transformation stretch tensor with the minimum distance. 

Author: Yunsu Park
Created: April 9 2025
Affiliation: University of California, Santa Barbara
Contact: yunsu@ucsb.edu
'''

import numpy as np
import os as os

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
    C = F.T @ F
    eig_val, eig_vec = np.linalg.eig(C)
    lam = np.sqrt(eig_val)
    U = lam[0]*np.outer(eig_vec[0],eig_vec[0]) + lam[1]*np.outer(eig_vec[1],eig_vec[1]) + lam[2]*np.outer(eig_vec[2],eig_vec[2])

    # Calculates the distance function defined by Chen et al. 
    U2 = U.T @ U
    F2 = F.T @ F
    A = np.linalg.inv(F2)- np.eye(3)
    distFunc = np.trace(A.T @ A)

    return distFunc, U

def loopDist(file_path, ref_files, def_files, reflat, deflat):
    '''
    Loops through each possible combination of correspondance matrix for each configuration and finds the three minimum distance

    Parameters:
        file_path (string): 
            directory for the folder with files with d
        ref_files (list [shape (*, 1)]): 
            directory for the files with files with d and p
        def_files (ndarray [shape (*, 1)]): 
            directory for the files with files with d and q
        reflat (ndarray [shape (3, 3)]): 
            lattice vector for the reference configuration
        reflat (ndarray [shape (3, 3)]): 
            lattice vector for the deformed configuration

    Returns:
        distFunc_stored (ndarray [shape (3, 1)]):
            distance function that is top three min
        U_stored (ndarray [shape (3, 3, 3)]): 
            stretch tensor with distance function that is top three min
        P_ref_stored (ndarray [shape (3, 3, 3)]): 
            correspondance matrix with distance function that is top three min
        P_def_stored (ndarray [shape (3, 3, 3)]): 
            correspondance matrix with distance function that is top three min
    '''
    print('Calculating minimum distance function')

    # Initialization
    P_ref_stored = np.zeros((3, 3, 3))
    P_def_stored = np.zeros((3, 3, 3))
    distFunc_stored = np.ones((3,1)) * 1e100
    U_stored = np.zeros((3, 3, 3))

    # Loops through the correspondance files
    for i in range(len(ref_files)):
        P_ref_file = np.load(os.path.join(file_path, ref_files[i]))
        
        for j in range(len(def_files)):
            P_def_file = np.load(os.path.join(file_path, def_files[j]))

            # Loops through the correspondance matrix
            for k in range(P_ref_file.shape[0]):
                P_ref = P_ref_file[k]
                #print(f'Calculating File: {i+1}/{len(ref_files)} --- Completed Matrix: {k+1}/{P_ref_file.shape[0]}')

                for l in range(P_def_file.shape[0]):
                    P_def = P_def_file[l]
                    distFunc, U = calcDist(reflat, deflat, P_ref, P_def)

                    # Find if this distFunc is better than one of the top 3
                    worst_idx = np.argmax(distFunc_stored)
                    if distFunc < distFunc_stored[worst_idx]:
                        distFunc_stored[worst_idx] = distFunc
                        P_ref_stored[worst_idx] = P_ref
                        P_def_stored[worst_idx] = P_def
                        U_stored[worst_idx] = U

    print('   Complete: the three lowest distance function is')
    print(' ' * 12, np.array2string(distFunc_stored, prefix=' ' * 12), '\n')
    return distFunc_stored, U_stored, P_ref_stored, P_def_stored

def saveDist(distFunc, U, P_ref, P_def):
    '''
    Saves the three minimum distance function and its associated stretch tensor and correspondance matrix for each configuration

    Parameters:
        distFunc_stored (ndarray [shape (3, 1)]):
            distance function that is top three min
        U_stored (ndarray [shape (3, 3, 3)]): 
            stretch tensor with distance function that is top three min
        P_ref_stored (ndarray [shape (3, 3, 3)]): 
            correspondance matrix with distance function that is top three min
        P_def_stored (ndarray [shape (3, 3, 3)]): 
            correspondance matrix with distance function that is top three min

    Returns:
        int: Always returns 0 to indicate completion.
    '''
    # Prints out top 3 min distance data
    print('Saving distance minimization data')
    print('   Distance function saved')
    print(' ' * 12, np.array2string(distFunc, prefix=' ' * 12), '\n')
    print('   Stretch Tensor saved')
    print(' ' * 12, np.array2string(U, prefix=' ' * 12), '\n')
    print('   Correspondance matrix of the reference configuration saved')
    print(' ' * 12, np.array2string(P_ref, prefix=' ' * 12), '\n')
    print('   Correspondance matrix of the deformed configuration saved')
    print(' ' * 12, np.array2string(P_def, prefix=' ' * 12), '\n')

    # Saves the data with unique name
    counter = 0
    folder_path = os.path.join(os.getcwd(), 'calc_data')

    while True:
        file_name = f'stored_results{counter}.npz'
        file_path = os.path.join(folder_path, file_name)
        if not os.path.exists(file_path):
            break
        counter += 1

    # Makes save directory if it doesnt exist
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
        
    np.savez(file_path,
            P_ref=P_ref,
            P_def=P_def,
            U=U,
            dist=distFunc)

    print(f'   COMPLETE: Calculated data stored as "stored_results{counter}.npz" \n')

    return 0