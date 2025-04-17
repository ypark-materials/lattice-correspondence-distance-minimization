import numpy as np
import os as os

from module import corrmat as cm

'''
Input
'''
d = 3

'''
Calculation
'''
# Generates correspondance matricies if it doesnt exist in file
cm.saveCorMat(d)

# Loads one of the the generated files
file = os.path.join(os.getcwd(),'data', f'data_d{d}', f'Pmat_d{d}_det1_0.npy')
P = np.load(file)