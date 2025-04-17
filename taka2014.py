import numpy as np
import os as os

from module import corrmat as cm
from module import distmin as dm
from module import crystallo as cr

'''
Input
'''
# Number of atoms/molecules in unit cell
p = 1 # reference phase
q = 2 # deformed phase

# Unit cell parameters of the reference configuration
abc_ref = np.array([5.03, 5.395, 7.202])
angle_ref = np.array([103.413, 100.269, 92.382])

# Unit cell parameters of the deformed configuration
abc_def = np.array([5.3663, 7.268, 10.16])
angle_def = np.array([104.149, 97.699, 92.382])

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

# Calculates the new unit cell parameters
abc_ref_n, abc_def_n, angle_ref_n, angle_def_n = dm.newLatt(reflat, deflat, P_ref, P_def)

'''
Output
'''
dm.saveDist(distFunc, U, P_ref, P_def, abc_ref_n, abc_def_n, angle_ref_n, angle_def_n)



