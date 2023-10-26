#!/usr/bin/evn python3

from sys import argv
from numpy import *

filename = str(argv[1])
data = loadtxt(filename)

# string energy and Pz
StringE = data[:, 0]*(cosh(data[:, 17]) - cosh(data[:, 13])
                      + cosh(data[:, 18]) - cosh(data[:, 14]))
StringPz = data[:, 0]*(sinh(data[:, 17]) - sinh(data[:, 13])
                       + sinh(data[:, 18]) - sinh(data[:, 14]))
# Remnant energy and Pz
RemnantE = data[:, 0]*(data[:, 15]*cosh(data[:, 13])
                       + data[:, 16]*cosh(data[:, 14]))
RemnantPz = data[:, 0]*(data[:, 15]*sinh(data[:, 13])
                        + data[:, 16]*sinh(data[:, 14]))
print(sum(StringE) + sum(RemnantE), sum(StringPz) + sum(RemnantPz))

