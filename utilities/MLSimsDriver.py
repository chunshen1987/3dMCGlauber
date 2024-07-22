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
import importlib
import string
import os
import random

def printHelp() -> None:
    print("Usage: {} nev path_to_parameters path_to_emulator_folder [output_file_name] [output_folder_name]".format(sys.argv[0]))
    print("nev: the number of initial stage events from 3DMCGlauber")
    print("path_to_parameters: The parameter dictionnary for 3DMCGlauber")
    print("path_to_emulator_folder: the path to the emulator folder containing the emulator object and the object definition.")
    print("output_file_name: [optional] output name for the output pickle dictionnary containing the predictions.")
    print("Default output_file_name: OUTPUT_{Proj}{Targ}_{roots}_{Nev}.dat")
    print("output_folder: [optional] Folder in which the initial stage data and predictions are sent to.")
    print("Default output_folder: {Proj}{Targ}_{roots}_{Nev}")

def dynamic_import(module_name, file_path):
    try:
        # Load the module from the given file path
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        return module
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def create_folder(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        return folder_path
    else:
        random_char = random.choice(string.ascii_letters)
        folder_path_ = folder_path+"_"+str(random_char)
        return create_folder(folder_path_)

def save_data(X, Arr, fname="out") -> None:
    """ Data is saved as a 2D numpy array.
        The first row is the X values
        The Nev remaining rows are the model predictions
    """
    OUT = np.vstack([X, Arr])
    fi = open(fname, 'wb')
    pickle.dump(OUT, fi)
    fi.close()

def Predict(ARR, path_to_emulator_object, out_fname="Predictions", out_Fname="PRED", Nmax=-1) -> None:
    """ 1) Load emulator pickle object. 
        The emulator object needs:
        - a predict method taking data as input and output predictions 
        2) Perform the prediction of the desired 2D data ARR
    """
    # Load emulator
    h = os.getcwd()
    os.chdir(path_to_emulator_object)
    try:
        with open("EM.dat", "rb") as pf:
            EM = pickle.load(pf)
    except FileNotFoundError:
        print("Can't find emulator object at path: "+path)
        exit(0)

    N = len(ARR) - 1 if Nmax == -1 else Nmax
    X, INPUT_DATA = ARR[0], ARR[1:N]
    pred = EM.predict(INPUT_DATA)
    save_data(np.linspace(np.min(X), np.max(X), pred.shape[1]), pred, fname=out_Fname+"/"+out_fname)
    os.chdir(h)

def Generate_3DMCGlauber_input(nev, paraDict, Folder_path):
    """ This function generates nev initial stage 
        events using 3DMCGlauber set on paraDict parameters
    """
    h = os.getcwd()
    os.chdir(os.path.abspath(Folder_path))
    command = "./3DMCG {} input -1 batch_density_output=1".format(nev)
    for ikey in paraDict.keys():
        command += " {}={}".format(ikey, paraDict[ikey])
    subprocess.run(command, shell = True, executable="/bin/bash")

    edArr = np.fromfile("ed_etas_distribution_N_72.dat", dtype="float32")
    nBArr = np.fromfile("nB_etas_distribution_N_72.dat", dtype="float32")
    nQArr = np.fromfile("nQ_etas_distribution_N_72.dat", dtype="float32")

    edArr = edArr.reshape(-1, 72)
    nBArr = nBArr.reshape(-1, 72)
    nQArr = nQArr.reshape(-1, 72)
    os.chdir(h)
    return edArr, nBArr, nQArr

def PrepareFolder(path_to_emulator, Fname_out):
    Folder_path = create_folder(path_to_emulator+"/"+Fname_out) 
    Fname_out = Folder_path.split("/")[-1]
    h = os.getcwd()
    shutil.copyfile(path_to_parameters, Folder_path+"/parameters.py")
    shutil.copyfile(h+"/input", Folder_path+"/input")
    os.system("ln -s "+h+"/3dMCGlb.e " +os.path.abspath(Folder_path)+"/3DMCG")
    os.system("ln -s "+h+"/Metropolis.e " +os.path.abspath(Folder_path)+"/Metropolis.e")
    os.system("ln -s "+h+"/tables " +os.path.abspath(Folder_path)+"/tables")
    return Folder_path, Fname_out

def CleanTemporary(Folder_path) -> None:
    os.system("rm -rf __pycache__")
    os.system("rm -rf "+os.path.abspath(Folder_path)+"/3DMCG")
    os.system("rm -rf "+os.path.abspath(Folder_path)+"/Metropolis.e")
    os.system("rm -rf "+os.path.abspath(Folder_path)+"/tables")
    os.system("rm -rf "+os.path.abspath(Folder_path)+"/input")
###############################################################################
### analysis starts from here ...
###############################################################################


# Parse user input
try:
    nev, path_to_parameters, path_to_emulator = int(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3])
except IndexError:
    printHelp()
    exit(0)
try:
    fname_out = str(sys.argv[4])
    need_fname = False 
except IndexError:
    need_fname = True 
try:
    Fname_out = str(sys.argv[5])
    need_Fname = False 
except IndexError:
    need_Fname = True 

# Dynamical import of the custom EMULATOR class
Emulator = dynamic_import("EMULATOR", os.path.abspath(path_to_emulator+"/EMULATOR.py"))

# Dynamical import of the 3DMCGlauber model parameters dictionnary 
parameters = dynamic_import("para_dict", path_to_parameters)
para_dict = parameters.para_dict

# Define output folder and output file name
if need_fname:
    fname_out = "OUTPUT_"+para_dict["Projectile"]+ para_dict["Target"]+"_"+str(para_dict["roots"])+"_"+str(nev)+".dat"
if need_Fname:
    Fname_out = para_dict["Projectile"]+ para_dict["Target"]+"_"+str(para_dict["roots"])+"_"+str(nev)

# Generate 3DMCGlauber initial conditions in Folder_path
# and generate predictions using Emulator Folder_path/EM.py 
# Defined by Class in Folder_path/EMULATOR.py <- where the predict function should be defined. 
# predict function should take 2D array (Nev, initial Neta = 72) and output 2D array (Nev, final Neta)
Folder_path, Fname_out = PrepareFolder(path_to_emulator, Fname_out)
e, B, Q = Generate_3DMCGlauber_input(nev, para_dict, Folder_path)
Predict(B, path_to_emulator, out_Fname=Fname_out, out_fname=fname_out, Nmax=-1)
CleanTemporary(Folder_path)
