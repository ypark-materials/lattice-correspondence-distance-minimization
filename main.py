import numpy as np
import os as os

from module import corrmat as cm
from module import distmin as dm
from module import crystallo as cr
from module import micromech as mm

'''
Input
'''
d = 1
p = 3
q = 2

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

cor_dir_path = os.getcwd()+'\\data_d' + str(d)
files = os.listdir(cor_dir_path)

ref_det_file = 'det' + str(p)
def_det_file = 'det' + str(q)

ref_files = [file for file in files if ref_det_file in file]
def_files = [file for file in files if def_det_file in file]

P_ref_stored = np.zeros(3,3,3)
P_def_stored = np.zeros(3,3,3)
distFunc_stored  = np.ones(1)*1e100
U_stored = np.zeros(3,3,3)

for i in range(len(ref_files)):
    P_ref = np.load(cor_dir_path +'\\' + ref_files[i])
    
    for i in range(len(def_files)):
        P_def = np.load(cor_dir_path +'\\' + def_files[i])

        distFunc, U = dm.calcDist(reflat, deflat, P_ref, P_def)

        if distFunc < distFunc_stored[0]:
            distFunc_stored[0] = distFunc
            U_stored[0,:,:] = U
            P_ref_stored[0,:,:] = P_ref
            P_def_stored[0,:,:] = P_def
        elif distFunc < distFunc_stored[1]:
            distFunc_stored[1] = distFunc
            U_stored[1,:,:] = U
            P_ref_stored[1,:,:] = P_ref
            P_def_stored[1,:,:] = P_def
        elif distFunc < distFunc_stored[2]:
            distFunc_stored[2] = distFunc
            U_stored[2,:,:] = U
            P_ref_stored[2,:,:] = P_ref
            P_def_stored[2,:,:] = P_def

print(distFunc_stored)
print(U_stored)
print(P_ref_stored)
print(P_def_stored)


#P = np.load(os.getcwd()+'\\data_d2\\Pmat_d2_det4_0.npy')
#F = np.dot(E_d @ P_d, np.linalg.inv(E_i @ P_i))
#F = np.linalg.inv(E_i @ P_i).dot(E_d @ P_d)
#print(F@P_d)




