import numpy as np
import os as os

from module import corrmat as cm
from module import distmin as dm
from module import crystallo as cr

'''
Input
'''
# Number of atoms/molecules in unit cell
p = 2 # reference phase
q = 1 # deformed phase

# Unit cell parameters of the reference configuration
abc_ref = np.array([15.7380, 9.2352, 15.7040])
angle_ref = np.array([90, 109.1209, 90])

# Unit cell parameters of the deformed configuration
abc_def = np.array([12.8946, 9.4837, 9.3384])
angle_def = np.array([90, 90, 90])

'''
Calculation
'''
# Largest edge ratio
d = round(max(np.array([*abc_ref, *abc_def])) / min(np.array([*abc_ref, *abc_def])))

# Generates correspondance matricies if it doesnt exist in file
cm.saveCorMat(d)

# Converts unit cell parameters (fractional coordinate) to lattice vectors (Cartesian coordinate)
reflat = cr.unit2vect(abc_ref, angle_ref)
deflat = cr.unit2vect(abc_def, angle_def)

# Reads the correspondance matrix files
file_path, ref_files, def_files = cm.readCorMat(d, p, q)

# Calculates the three minimum distance functions
distFunc, U, P_ref, P_def = dm.loopDist(file_path, ref_files, def_files, reflat, deflat)

# Calculates the new unit cell parameters
abc_ref_new, abc_def_new, angle_ref_new, angle_def_new = dm.newLatt(reflat, deflat, P_ref[1,:,:], P_def[1,:,:])

'''
Output
'''
# Saves the correspondance matrix, stretch tensor, distance function value, and unit cell parameters of the reduced parameters
dm.saveDist(distFunc, U, P_ref, P_def, abc_ref_new, abc_def_new, angle_ref_new, angle_def_new)



