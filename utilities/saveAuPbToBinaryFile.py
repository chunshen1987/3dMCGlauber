"""
    This is an example code to save a python data array to a binary file
    and c++ program can read from it. The data is saved with single float
    precision (float32).
"""

import array
import numpy as np

data = np.loadtxt("../tables/pb208-1.dat")
A = 208

data = data[:, :4]
data = data.reshape(-1, 4*A)
print(data.shape)

nev = 0
outputFile = open("Pb208.bin.in", "wb")
for config_i in data:
    float_array = array.array('f', config_i)
    float_array.tofile(outputFile)
    nev += 1
    if nev > 3000: break
outputFile.close()
