#!/usr/bin/env python3

import matplotlib.pyplot as plt
from numpy import *


d_quark = loadtxt("check_sampled_single_d_quark_distribution.dat")
u_quark = loadtxt("check_sampled_single_u_quark_distribution.dat")

sampled_quark = loadtxt("check_sampled_valence_quarks_distribution.dat")

# d quark
fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.83, 0.83])
plt.plot(d_quark[:, 0], d_quark[:, 1], '-k', label="d pdf")
plt.plot(d_quark[:, 0], d_quark[:, 2], "+r", label="d samples")
plt.plot(sampled_quark[:, 0], sampled_quark[:, 1], "--b",
         label="d samples with 0.95 < sum < 1")
plt.legend()
plt.show()

# u quark
fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.83, 0.83])
plt.plot(u_quark[:, 0], u_quark[:, 1], '-k', label="u pdf")
plt.plot(u_quark[:, 0], u_quark[:, 2], "+r", label="u samples")
plt.plot(sampled_quark[:, 0], sampled_quark[:, 2], "--b",
         label="u samples with 0.95 < sum < 1")
plt.legend()
plt.show()
