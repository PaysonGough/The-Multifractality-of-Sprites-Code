

import matplotlib.pyplot as plt
import numpy as np

mat = plt.imread("sprite2.tif")

x, y = np.shape(mat)

big_pix = 0

for i in range(x):
    for j in range(y):
        if mat[i, j] > big_pix:
            big_pix = mat[i, j]
            
print(big_pix)