#!/usr/bin/env python3

import matplotlib.pyplot as plt
from numpy import *



sampled_quark = loadtxt("check_deformed_Woods_Saxon_sampling.dat")

# d quark
fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.83, 0.83])
plt.plot(sampled_quark[:, 0], sampled_quark[:, 1], "--b",
         label="d samples with 0.95 < sum < 1")
plt.plot(sampled_quark[:, 0], sampled_quark[:, 2], "r+",
         label="d samples with 0.95 < sum < 1")
plt.legend()
plt.savefig("ddd")
