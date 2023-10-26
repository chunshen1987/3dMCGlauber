#! /usr/bin/env python3
"""
     This script performs event averaging for particle 
     spectra and anisotropic flow coefficients calculated 
     from event-by-event simulations

     v_n is analyzed up to n = 6

     Format for particle_XXX_vndata.dat file:
     n_order  real_part  real_part_err  imag_part  imag_part_err

     Format for particle_XXX_vndata_diff.dat file:
     pT(GeV)  pT_err(GeV)  dN/(2pi dy pT dpT)(GeV^-2)  dN/(2pi dy pT dpT)_err(GeV^-2)
     vn_real  vn_real_err  vn_imag  vn_imag_err

     All the errors are only statistic errors
"""

from sys import argv, exit
from os import path, mkdir
from glob import glob
import numpy as np
import h5py
import shutil
from scipy.optimize import curve_fit


def LinearFunc(x, a, b):
    return a * x + b


def calculate_rn_eta(eta_array, eta_min, eta_max, dN_array, vn_array,
                     outputFileName):
    """
        This function computes the longitudinal factorization breaking ratios
        for all n passed from vn_array
            eta, rn(eta), r_nn(eta)
        the reference flow angles are defined in [eta_min, eta_max]
        and [-eta_max, -eta_min]
    """
    nev, neta = dN_array.shape
    dN_array = dN_array.reshape((nev, 1, neta))
    Qn_array = vn_array
    nQn = Qn_array.shape[1]

    # calculate the reference flow vector for every event
    eta_b_min    = np.abs(eta_min)
    eta_b_max    = np.abs(eta_max)
    eta_ref1_tmp = np.linspace(eta_b_min, eta_b_max, 16)
    eta_ref2_tmp = np.linspace(-eta_b_max, -eta_b_min, 16)
    Qn_ref1      = []
    Qn_ref2      = []
    for iev in range(nev):
        dN1_interp = np.interp(eta_ref1_tmp, eta_array, dN_array[iev, 0, :])
        dN2_interp = np.interp(eta_ref2_tmp, eta_array, dN_array[iev, 0, :])
        Qn_ref1_vec = []
        Qn_ref2_vec = []
        for iorder in range(nQn):
            Qn1_interp = np.interp(eta_ref1_tmp, eta_array,
                                   Qn_array[iev, iorder, :])
            Qn2_interp = np.interp(eta_ref2_tmp, eta_array,
                                   Qn_array[iev, iorder, :])
            Qn_ref1_vec.append(np.sum(dN1_interp*Qn1_interp)
                               /(np.sum(dN1_interp) + 1e-15))
            Qn_ref2_vec.append(np.sum(dN2_interp*Qn2_interp)
                               /(np.sum(dN2_interp) + 1e-15))
        Qn_ref1.append(Qn_ref1_vec)
        Qn_ref2.append(Qn_ref2_vec)
    Qn_ref1 = np.array(Qn_ref1).reshape((nev, nQn, 1))
    Qn_ref2 = np.array(Qn_ref2).reshape((nev, nQn, 1))

    rn_num  = np.real(Qn_array[:, :, ::-1]*np.conj(Qn_ref1))
    rn_den  = np.real(Qn_array*np.conj(Qn_ref1))
    rnn_num = np.real((Qn_ref2*np.conj(Qn_array))
                      *(Qn_array[:, :, ::-1]*np.conj(Qn_ref1)))
    rnn_den = np.real((Qn_ref2*np.conj(Qn_array[:, :, ::-1]))
                      *(Qn_array*np.conj(Qn_ref1)))

    # compute the error using jack-knife
    rn_array  = np.zeros([nev, nQn, neta])
    rnn_array = np.zeros([nev, nQn, neta])
    for iev in range(nev):
        array_idx      = [True]*nev
        array_idx[iev] = False
        array_idx      = np.array(array_idx)
        rn_ev = (  np.mean(rn_num[array_idx], axis=0)
                 /(np.mean(rn_den[array_idx], axis=0) + 1e-15))
        rnn_ev = (  np.mean(rnn_num[array_idx], axis=0)
                  /(np.mean(rnn_den[array_idx], axis=0) + 1e-15))
        rn_array[iev, :, :] = rn_ev
        rnn_array[iev, :, :] = rnn_ev
    rn_mean  = np.mean(rn_array, axis=0)
    rn_err   = np.sqrt((nev - 1.)/nev*np.sum((rn_array - rn_mean)**2., axis=0))
    rnn_mean = np.mean(rnn_array, axis=0)
    rnn_err  = (np.sqrt((nev - 1.)/nev
                        *np.sum((rnn_array - rnn_mean)**2., axis=0)))

    # perform a linear fit for r_n
    rn_slope = []
    rnn_slope = []
    idx = abs(eta_array) < 1.
    for iorder in range(nQn):
        popt, pcov = curve_fit(LinearFunc, eta_array[idx],
                               rn_mean[iorder, idx],
                               sigma = rn_err[iorder, idx],
                               method = "dogbox")
        rn_slope.append([popt[0], np.sqrt(pcov[0, 0])])
        popt, pcov = curve_fit(LinearFunc, eta_array[idx],
                               rnn_mean[iorder, idx],
                               sigma = rnn_err[iorder, idx],
                               method = "dogbox")
        rnn_slope.append([popt[0], np.sqrt(pcov[0, 0])])

    f = open(outputFileName, 'w')
    f.write("#eta  rn(eta)  rn_err(eta)  rnn(eta)  rnn_err(eta) (n=2, 3)\n")
    for ieta in range(len(eta_array)-1):
        f.write("%.10e  " % eta_array[ieta])
        for iorder in range(nQn):
            f.write("%.10e  %.10e  %.10e  %.10e  "
                    % (rn_mean[iorder, ieta], rn_err[iorder, ieta],
                       rnn_mean[iorder, ieta], rnn_err[iorder, ieta]))
        f.write("\n")
    f.close()
    return rn_slope, rnn_slope



###############################################################################
### analysis starts from here ...
###############################################################################

ecc2Data = np.fromfile("ecc_ed_n_2_Neta_72.dat", dtype="float32")
ecc2Data = np.nan_to_num(ecc2Data.reshape(-1, 72))
ecc3Data = np.fromfile("ecc_ed_n_3_Neta_72.dat", dtype="float32")
ecc3Data = np.nan_to_num(ecc3Data.reshape(-1, 72))
edData = np.fromfile("ed_etas_distribution_N_72.dat", dtype="float32")
edData = np.nan_to_num(edData.reshape(-1, 72))

eta_array = ecc2Data[0, :]
ecc2_array = ecc2Data[1::2, :] + 1j*ecc2Data[2::2, :]
ecc3_array = ecc3Data[1::2, :] + 1j*ecc3Data[2::2, :]
nev, neta = ecc2_array.shape
eccn_array = np.zeros([nev, 2, neta]) + 1j*np.zeros([nev, 2, neta])
eccn_array[:, 0, :] = ecc2_array
eccn_array[:, 1, :] = ecc3_array
dEdetas_array = edData[1:, :]

ecc_rnSlopeFile = open("ecc_rnSlope.dat", "w")
ecc_rnSlopeFile.write("#cen  rn_slope  rn_slope_err (n=2, 3)\n")
ecc_rnnSlopeFile = open("ecc_rnnSlope.dat", "w")
ecc_rnnSlopeFile.write("#cen  rnn_slope  rnn_slope_err (n=2, 3)\n")

# calculate the longitudinal flow decorrelation with STAR cut
etaMin = 2.5; etaMax = 4

output_filename = path.join("ecc_rn_eta.dat")
ecc_rnSlope, ecc_rnnSlope = calculate_rn_eta(
        eta_array, etaMin, etaMax, dEdetas_array, eccn_array,
        output_filename)
centralityCenter = 0.
ecc_rnSlopeFile.write("%.2f  " % centralityCenter)
for ir, ir_err in ecc_rnSlope:
    ecc_rnSlopeFile.write("%.6e  %.6e  " % (ir, ir_err))
ecc_rnSlopeFile.write("\n")
ecc_rnnSlopeFile.write("%.2f  " % centralityCenter)
for ir, ir_err in ecc_rnnSlope:
    ecc_rnnSlopeFile.write("%.6e  %.6e  " % (ir, ir_err))
ecc_rnnSlopeFile.write("\n")

ecc_rnSlopeFile.close()
ecc_rnnSlopeFile.close()
print("Analysis is done.")

