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

data = np.load('your_file.npy')
print(data)