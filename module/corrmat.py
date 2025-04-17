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
    
    Pdata0 = np.array([[[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]])
    filename = f'Pmat_d{d}_det{0}_{0}.npy'
    file_path = os.path.join(dir_path, filename)
    np.save(file_path, Pdata0.astype(np.int8))

    print(f'   COMPLETE: Correspondance matricies for d = {d} saved') 
    print(f'             Saved at file path "{dir_path}"')
    print(f'             Total number of files = {sum(countSave)+9} \n')
    
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
    print('Checking correspondance matrix file')
    if not os.path.exists(dir_path):
        # Makes directory and generates correspondance matrix file if path doesnt exist
        print('File does not exist')
        print('Generating correspondance matrix file')
        os.mkdir(dir_path)
        genCorMat(dir_path, d)
    else:
        # Skips file generation if path exists
        print('   COMPELTE: Correspondance matrix files for d = ', d, ' already exits at "', dir_path, '" \n')

    return 0

def readCorMat(d, p, q):
    '''
    Reads the correspondance matrix file with the given value of d

    Parameters:
        d (integer):
            maximum integer difference between the length between the unit cell parameter of the reference to transformed configuration
        p (integer):
            determinant of reference configuration
        q (integer):
            determinant of deformed configuration
            
    Returns:
        file_path (string): 
            directory for the folder with files with d
        ref_files (list [shape (*, 1)]): 
            directory for the files with files with d and p
        def_files (ndarray [shape (*, 1)]): 
            directory for the files with files with d and q
    '''

    print('Reading correspondance matrix data files')

    # Finds files for d
    file_path = os.path.join(os.getcwd(), 'data', f'data_d{d}')
    files = os.listdir(file_path)
    # Checks which files to read
    m = p/q
    if int(m) == m:
        p = int(0)
        q = int(m)
        print(f'   READING: m = {m}, only q read (reference has more atoms/molecules)')
    elif int(1/m) == 1/m:
        p = int(1/m)
        q = int(0)
        print(f'   READING: 1/m = {1/m}, only p read (deformed has more atoms/molecules)')
    else:
        print(f'   READING: m = {m}, so both determinant of p and q read')

    # Finds files with the corresponding determinate for the reference and deformed
    ref_files = [file for file in files if f'det{p}' in file]
    def_files = [file for file in files if f'det{q}' in file]

    print(f'   COMPLETE: reference files read {ref_files}')
    print(f'             deformed files read {def_files}\n')
    return file_path, ref_files, def_files