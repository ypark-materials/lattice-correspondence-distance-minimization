'''
crystallo.py

Crystallo module contains functions for doing crystallographic operations

Author: Yunsu Park
Created: Mar 6 2025
Affiliation: University of California, Santa Barbara
Contact: yunsu@ucsb.edu
'''

import numpy as np

def unit2vect(abc, angle):
    """
    Converts unit cell parameters (fractional coordinate) to lattice vectors (Cartesian coordinate)

    Parameters:
        abc (ndarray [shape (1, 3)]): 
            unit cell parameter lengths a, b, c
        angle (ndarray [shape (1, 3)]): 
            unit cell parameter angles alpha, beta, gamma

    Returns:
        latvec (ndarray [shape (3, 3)]):
            lattice matrix with a as the parallel vector
            [[a1, b1, c1], 
             [0,  b2, c2], 
             [0,  0,  c3]]
    """

    # Unit cell parameters length a, b, c
    a, b, c = [l for l in abc]

    # Unit cell parameter angles alpha, beta, gamma
    alp, bet, gam = [deg * np.pi/180 for deg in angle]

    # Lattaice matrix elements
    a1 = a
    b1 = b*np.cos(gam)
    b2 = b*np.sin(gam)
    c1 = c*np.cos(bet)
    c2 = (c/np.sin(gam)) * (np.cos(alp)-np.cos(bet)*np.cos(gam))
    c3 = (c/np.sin(gam)) * np.sqrt(1-(np.cos(alp))**2-(np.cos(bet))**2-(np.cos(gam))**2 + 2*np.cos(alp)*np.cos(bet)*np.cos(gam))

    # Lattaice matrix 
    latvec = np.array([[a1, b1, c1],
                       [0,  b2, c2],
                       [0,  0,  c3]])
    
    if np.isnan(latvec).any():
        print('ERROR: there is a nan value in the lattice vector')
        print(latvec)

    return latvec

def frac2cart(milfrac, latvec):
    """
    Converts fractional coordinates to cartesian coordinates

    Parameters:
        milfrac (ndarray [shape (1, 3)]): 
            miller indice with fractional coordinate (h, k, l)
        latvec (ndarray [shape (3, 3)]):
            lattice matrix with a as the parallel vector
            [[a1, b1, c1], 
             [0,  b2, c2], 
             [0,  0,  c3]]

    Returns:
        milvec (ndarray [shape (3, 1)]):
            miller indice in vector form in cartesian coordinate
    """
    # Millar indice (h k l)
    h, k, l = [ind for ind in milfrac]

    # latticce vectors
    e1 = latvec[:,0]
    e2 = latvec[:,1]
    e3 = latvec[:,2]

    milvec = h*e1 + k*e2 + l*e3

    return milvec

def cart2frac(milvec, latvec):
    """
    Converts Cartesian coordinates to fractional coordinates

    Parameters:
        milvec (ndarray [shape (3, 1)]):
            miller indice in vector form in cartesian coordinate    

        latvec (ndarray [shape (3, 3)]):
            lattice matrix with a as the parallel vector
            [[a1, b1, c1], 
             [0,  b2, c2], 
             [0,  0,  c3]]

    Returns:
        milfrac (ndarray [shape (1, 3)]): 
            miller indice with fractional coordinate (*, *, *)
    """

def deltang(vec1, vec2):
    """
    Calculates the angle between two vectors

    Parameters:
        vec1 (ndarray [shape (3, 1)]):
            vector 1

        vec2 (ndarray [shape (3, 1)]):
            vector 2

    Returns:
        theta (float): 
            angle in degrees between vec1 and vec2
    """
        
    # Normalize the vectors
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)

    # Checks for zero length vectors
    if norm1 == 0 or norm2 == 0:
        raise ValueError("One of the vectors is zero-length.")
    
    cos = np.dot(vec1, vec2) / (norm1*norm2)

    # Clip to avoid numerical issues at theta ~ n*(pi/2)
    cos = np.clip(cos, -1.0, 1.0)

    theta = np.arccos(cos) * 180/np.pi

    return theta

def minang(n, n_n, milvec):
    """
    Finds the minimum angle between normal vectors

    Parameters:
        n (ndarray [shape (3, 1)]):
            analytical vector normal to habit plane from kinematic compatibility with kappa = 1

        n_n (ndarray [shape (3, 1)]):
            analytical vector normal to habit plane from kinematic compatibility with kappa = -1

        milvec (ndarray [shape (3, 1)]):
            experimental vector normal to habit plane    
        
    Returns:
        theta (float): 
            minimum angle between the analytical and experimental normal vector
    """

    # 4 possible theta values
    theta1 = deltang(n, milvec)
    theta2 = deltang(n, -milvec)
    theta3 = deltang(n_n, milvec)
    theta4 = deltang(n_n, -milvec)

    # Finds the minimun of the 4 theta values
    theta = np.min([theta1, theta2, theta3, theta4])

    return theta