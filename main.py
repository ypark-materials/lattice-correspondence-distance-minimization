import numpy as np
import os as os

from module import corrmat as cm
from module import distmin as dm
from module import crystallo as cr

'''
Input
'''
d = 1
p = 5
q = 5

# Unit cell parameters of the reference configuration
abc_ref = np.array([7.6, 8.58, 17.23])
angle_ref = np.array([78.22, 86.71, 72.1])

# Unit cell parameters of the deformed configuration
abc_def = np.array([7.76, 7.74, 16.94])
angle_def = np.array([77.80, 88.5, 82.2])

'''
Calculation
'''
# Generates correspondance matricies if it doesnt exist in file
cm.saveCorMat(d)

# Converts unit cell parameters (fractional coordinate) to lattice vectors (Cartesian coordinate)
reflat = cr.unit2vect(abc_ref, angle_ref)
deflat = cr.unit2vect(abc_def, angle_def)

# Reads the correspondance matrix files
file_path, ref_files, def_files = cm.readCorMat(d, p, q)

# Calculates the three minimum distance functions
distFunc, U, P_ref, P_def = dm.loopDist(file_path, ref_files, def_files, reflat, deflat)

'''
Output
'''
dm.saveDist(distFunc, U, P_ref, P_def)



