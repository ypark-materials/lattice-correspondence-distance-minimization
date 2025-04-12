'''
corrmat.py

Corrmat module contains functions for generating all possible correspondance matrix

Author: Yunsu Park
Created: April 4 2025
Affiliation: University of California, Santa Barbara
Contact: yunsu@ucsb.edu
'''

import numpy as np
import os as os

def genCorMat(dir_path, d):
    '''
    Generates all possible correspondance matrix given d and saves it in a file

    Parameters:
        dir_path (string): 
            path that will save the generated correspondance matrix
        d (integer): 
            maximum integer difference between the length between the unit cell parameter of the reference to transformed configuration

    Returns:
        int: Always returns 0 to indicate completion.
    '''
    
    # Initialize the correspondance matrix file
    elements = list(range(-d, d + 1))                   #all posible range of d
    n = 100000                                          #number of correspondance matrix in 1 file
    Pdata = np.tile(np.eye(3), (n, 8, 1, 1))            #initialized correspondance matrix 

    # Loops through all possible combinations of d in the elements of the correspondance matrix
    count = np.zeros(8, dtype=int)
    countSave = np.zeros(8, dtype=int)
    for P11 in elements:
        for P12 in elements:
            for P13 in elements:
                for P21 in elements:
                    for P22 in elements:
                        for P23 in elements:
                            for P31 in elements:
                                for P32 in elements:
                                    for P33 in elements:

                                        Pmat = np.array([[P11, P12, P13],
                                                         [P21, P22, P23],
                                                         [P31, P32, P33]])

                                        detP = np.linalg.det(Pmat)
                                        detPint = int(round(detP))

                                        if detP < 8.1 and np.isclose(detP, detPint) and detP > 0:
                                            Pdata[count[detPint-1], detPint-1, :, :] = Pmat.astype(np.int8)
                                            count[detPint-1] += 1

                                            if count[detPint-1] > n-1:
                                                filename = f'Pmat_d{d}_det{detPint}_{countSave[detPint-1]}.npy'
                                                file_path = os.path.join(dir_path, filename)
                                                np.save(file_path, Pdata[:, detPint-1, :, :].astype(np.int8))

                                                countSave[detPint-1] += 1
                                                count[detPint-1] = 0

    for i in range(8):
        Pdata_unique = np.unique(Pdata[:, i, :, :], axis=0)

        filename = f'Pmat_d{d}_det{i+1}_{countSave[i]}.npy'
        file_path = os.path.join(dir_path, filename)
        np.save(file_path, Pdata_unique.astype(np.int8))

    print('COMPLETE: Correspondance matricies for d = ', d, 'saved') 
    print('          Saved at file path "', dir_path, '"')
    print('          Total number of files = ', sum(countSave))
    
    return 0

def saveCorMat(d):
    '''
    Determines the directory path the correspondance matrix file will be saved in which will be "./data/data_d*"

    Parameters:
        d (integer): 
            maximum integer difference between the length between the unit cell parameter of the reference to transformed configuration

    Returns:
        int: Always returns 0 to indicate completion.
    '''

    # Fetches the directory path that the file will be saved in
    dir_path = os.path.join(os.getcwd(), 'data', f'data_d{d}')

    # Checks if correspondance matrix file exists
    if not os.path.exists(dir_path):
        # Makes directory and generates correspondance matrix file if path doesnt exist
        os.mkdir(dir_path)
        genCorMat(dir_path, d)
    else:
        # Skips file generation if path exists
        print('COMPELTE: Correspondance matrix files for d = ', d, ' already exits at "', dir_path, '"')

    return 0

def readCorMat(d, det, n):
    '''
    Reads the correspondance matrix file with the given value of d

    Parameters:
        d (integer):
            maximum integer difference between the length between the unit cell parameter of the reference to transformed configuration
        det (integer):
            
    Returns:
        corrmat (narray [shape (*, 3, 3)]):
            correspondance matrix where * is the index of the matrix with elements 3, 3
        check (integer):
            checker for when all files for a given d and det are met
    '''

    dir_path = os.getcwd()
    file_path = os.path.join(dir_path, 'data', f'data_d{d}', f'Pmat_d{d}_det{det}_{n}.npy')

    if os.path.exists(file_path):
        print('READING:')
        print(' ', f'Pmat_d{d}_det{det}_{n}.npy')

        corrmat = np.load(file_path)
        check = 1
    else:
        print('READING:')
        print(' completed')
        corrmat = 0
        check = 0

    return corrmat, check