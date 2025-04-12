# %%
from module import corrmat as cm

d = 3
det = 1

check = 1
n = 0
while check == 1:
    data, check = cm.readCorMat(d, det, n)
    n += 1
    print(data)


