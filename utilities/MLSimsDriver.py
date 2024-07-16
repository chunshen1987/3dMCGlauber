#!/usr/bin/env python3

"""
    This script drives the multi-stage simulations for
    3D-Glauber + ML emulators
    to predict high statistics final-state observables
"""
import shutil
import subprocess, sys
import numpy as np
import pickle
import sys
import tensorflow as tf

def save_data(X, Arr, fname="out") -> None:
    """ Data is saved as a 2D numpy array.
        The first row is the X values
        The Nev remaining rows are the model predictions
    """
    OUT = np.vstack([X, Arr])
    fi = open(fname, 'wb')
    pickle.dump(OUT, fi)
    fi.close()

def Predict(ARR, path_to_emulator_object, out_fname="Predictions", Nmax=-1) -> None:
    """ 1) Load emulator pickle object. 
        The emulator object needs:
        - a predict method taking data as input and output predictions 
        2) Perform the prediction of the desired 2D data ARR
    """
    # Load emulator
    try:
        with open(path, "rb") as pf:
            EM = pickle.load(pf)
    except FileNotFoundError:
        print("Can't find emulator object at path: "+path)
        exit(0)

    N = len(ARR) - 1 if Nmax == -1 else Nmax
    X, INPUT_DATA = ARR[0], ARR[1:N]
    save_data(X, EM.predict(INPUT_DATA), fname=out_fname)

def generate_3DMCGlauber_input(nev, paraDict):
    """ This function generates nev initial stage 
        events using 3DMCGlauber set on paraDict parameters
    """
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
    return edArr, nBArr, nQArr


def printHelp() -> None:
    print("Usage: {} nev path_to_emulator_object output_file_name".format(sys.argv[0]))
    print("nev: the number of initial stage events from 3DMCGlauber")
    print("path to emulator: the path to the emulator object")
    print("output_file_name: output name for the output pickle dictionnary containing the predictions.")


###############################################################################
### analysis starts from here ...
###############################################################################

para_dict = {
    'Projectile':  "Au",         # projectile nucleus name
    'Target'    :  "Au",         # target nucleus name
    'roots'     :  200,         # collision energy (GeV)
    'seed'      :   -1,          # random seed (-1: system)
    'baryon_junctions': 1,       # 0: baryon number assumed to be at string end
                                 # 1: baryon number transported assuming baryon
                                 # junctions (at smaller x)
                                 # see arXiv:nucl-th/9602027
    'lambdaB': 1.0,              # parameter the controls the strength of
    'lambdaQ': 1.0,
                                 # the baryon junction stopping
    'lambdaBs': 1.0,             # fraction of single-to-double string junction stopping
    'lambdaQs': 1.0,
    'baryonInStringProb': 0.1,
    'electric_junctions':  1,
    'integer_electric_charge': 1,
    'electricChargeInStringProb': 0.0,


    'BG': 17.095,
    'shadowing_factor':  0.145,     # a shadowning factor for producing strings from multiple scatterings
    'rapidity_loss_method': 4,
    'remnant_energy_loss_fraction': 0.611,      # nucleon remnants energy loss fraction (fraction of string's y_loss) [0, 1]
    'ylossParam4At2': 1.467,
    'ylossParam4At4': 1.759,
    'ylossParam4At6': 2.260,
    'ylossParam4At10': 3.262,
    'ylossParam4var': 0.356,
    'evolve_QCD_string_mode': 4,        # string evolution mode
                                        # 1: deceleration with fixed rapidity loss (m/sigma = 1 fm, dtau = 0.5 fm)
                                        # 2: deceleration with LEXUS sampled rapidit loss (both dtau and sigma fluctuate)
                                        # 3: deceleration with LEXUS sampled rapidit loss (m/sigma = 1 fm, dtau fluctuates)
                                        # 4: deceleration with LEXUS sampled rapidit loss (dtau = 0.5 fm, m/sigma fluctuates)
}

try:
    nev, path_to_emulator = int(sys.argv[1]), str(sys.argv[2])
except IndexError:
    printHelp()
    exit(0)
try:
    fname_out = str(sys.argv[3])
except IndexError:
    fname_out = "OUT_DICT.dat"

e, B, Q = generate_3DMCGlauber_input(nev, paraDict)
Predict(B, path_to_emulator, out_fname=fname_out, Nmax=-1)
