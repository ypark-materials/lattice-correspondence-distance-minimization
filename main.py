import numpy as np
import os as os

from module import corrmat as cm
from module import distmin as dm
from module import crystallo as cr

'''
Input
'''
# Number of atoms/molecules in unit cell
p = 4 # reference phase
q = 2 # deformed phase

# Unit cell parameters of the reference configuration
abc_ref = np.array([7.381, 11.755, 15.94])
angle_ref = np.array([102.912, 92.025, 100.595])

# Unit cell parameters of the deformed configuration
abc_def = np.array([6.0552, 7.0297, 15.969])
angle_def = np.array([96.315, 93.979, 90.279])

# Largest edge ratio
d = round(max(np.array([*abc_ref, *abc_def])) / min(np.array([*abc_ref, *abc_def])))

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



