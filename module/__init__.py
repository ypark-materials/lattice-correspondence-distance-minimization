'''
fincorrmat

A library for finding the lattice correspondance with the shortest distance during transformation.
Distance minimization lgorithm used in Chen et al. and Zhang et al. utilized.

Author: Yunsu Park
Version: 0.1.0
'''

__author__ = 'Yunsu Park'
__version__ = '0.1.0'

# --- Crystallography tools ---
from .crystallo import unit2vect

# --- Correspondance matrix generation tools ---
from .corrmat import saveCorMat

# --- Distance minimization tools ---
from .distmin import calcDist

__all__ = ['unit2vect',
           'saveCorMat',
           'calcDist']
