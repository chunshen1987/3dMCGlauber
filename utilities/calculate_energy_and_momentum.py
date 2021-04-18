#!/usr/bin/evn python3

from sys import argv
from numpy import *

filename = str(argv[1])
data = loadtxt(filename)

E = data[:, 0]*(data[:, 15]*cosh(data[:, 17]) + data[:, 16]*cosh(data[:, 18]))
Pz = data[:, 0]*(data[:, 15]*sinh(data[:, 17]) + data[:, 16]*sinh(data[:, 18]))
print(sum(E), sum(Pz))

