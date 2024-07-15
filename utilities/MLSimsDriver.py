#!/usr/bin/env python

"""
    This script drives the multi-stage simulations for
    3D-Glauber + ML emulators
    to predict high statistics final-state observables
"""

import shutil
import subprocess, sys
import numpy as np

paraDict = {
    "Projectile": "Au",
    "Target": "Au",
    "BG": 6,
    "b_min": 7.,
    "b_max": 10.,
}


def printHelp() -> None:
    print("Usage: {} nev".format(sys.argv[0]))


try:
    nev = int(sys.argv[1])
except:
    printHelp()
    exit(0)

command = "./3dMCGlb.e {} input -1 batch_density_output=1".format(nev)
for ikey in paraDict.keys():
    command += " {}={}".format(ikey, paraDict[ikey])

subprocess.run(command, shell = True, executable="/bin/bash")

edArr = np.fromfile("ed_etas_distribution_N_72.dat", dtype="float32")
nBArr = np.fromfile("nB_etas_distribution_N_72.dat", dtype="float32")
nQArr = np.fromfile("nQ_etas_distribution_N_72.dat", dtype="float32")

edArr = edArr.reshape(-1, 72)
nBArr = nBArr.reshape(-1, 72)
nQArr = nQArr.reshape(-1, 72)

