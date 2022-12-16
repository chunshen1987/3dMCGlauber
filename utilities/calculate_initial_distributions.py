#!/usr/bin/env python3

import sys
from os import path
from glob import glob

import numpy as np

HBARC = 0.19733

neta = 500
eta = np.linspace(-8, 8, neta)
deta = eta[1] - eta[0]

sigma_eta = 0.5
norm_eta = 1./(np.sqrt(2.*np.pi)*sigma_eta)
etaDis = 10.*sigma_eta

dEdetas = np.zeros([neta])
dNbdetas = np.zeros([neta])
dNbdy = np.zeros([neta])
filename_pattern_list = "strings_event_%d.dat"

work_folder = path.abspath(sys.argv[1])
fileList = glob(path.join(work_folder, "strings_event_*.dat"))
iev = 0; nev = len(fileList)
for file_i in fileList:
    stringList = np.loadtxt(file_i)
    if stringList.ndim == 1:
        stringList = string_list.reshape(1, len(stringList))
    print("processing event %d/%d ..." % (iev+1, nev))
    for istring in range(len(stringList)):
        # calcualte energy distribution
        mass = stringList[istring, 0]
        frac_l = stringList[istring, 15]
        frac_r = stringList[istring, 16]
        y_l = stringList[istring, 13]
        y_r = stringList[istring, 14]
        eta_l = stringList[istring, 11]
        eta_r = stringList[istring, 12]
        y_l_i = stringList[istring, 17]
        y_r_i = stringList[istring, 18]

        Estring = mass*(  np.cosh(y_l_i) + np.cosh(y_r_i)
                        - np.cosh(y_l) - np.cosh(y_r))
        EremL = frac_l*mass*np.cosh(y_l)
        EremR = frac_r*mass*np.cosh(y_r)

        # energy density
        f1_ed = np.zeros([neta])
        idx = ((eta < eta_r) & (eta > eta_l))       # inside the string
        y_eta = y_l + (y_r - y_l)/(eta_r - eta_l)*(eta - eta_l)
        f1_ed[idx] = np.cosh(y_eta[idx])
        idx = ((eta > eta_r) & (eta - eta_r < etaDis))
        f1_ed[idx] += (np.exp(-(eta[idx] - eta_r)**2./(2.*sigma_eta**2.))
                       *np.cosh(y_r))
        idx = ((eta < eta_l) & (eta_l - eta < etaDis))
        f1_ed[idx] += (np.exp(-(eta[idx] - eta_l)**2./(2.*sigma_eta**2.))
                       *np.cosh(y_l))
        dEdetas += Estring*f1_ed/(sum(f1_ed)*deta)

        f1_ed = np.zeros([neta])
        idx = (abs(eta - eta_r) < etaDis)
        f1_ed[idx] = (np.exp(-(eta[idx] - eta_r)**2./(2.*sigma_eta**2.))
                      *np.cosh(y_r))
        dEdetas += EremR*f1_ed/(sum(f1_ed)*deta)

        f1_ed = np.zeros([neta])
        idx = (abs(eta - eta_l) < etaDis)
        f1_ed[idx] = (np.exp(-(eta[idx] - eta_l)**2./(2.*sigma_eta**2.))
                      *np.cosh(y_l))
        dEdetas += EremL*f1_ed/(sum(f1_ed)*deta)

        # calcualte baryon distribution
        nB_eta_l = stringList[istring, 19]
        nB_eta_r = stringList[istring, 20]
        nB_y_l = stringList[istring, 21]
        nB_y_r = stringList[istring, 22]
        nB_frac_l = stringList[istring, 23]
        nB_frac_r = stringList[istring, 24]

        idx = (abs(eta - nB_eta_l) < etaDis)
        dNbdetas[idx] += (
            nB_frac_l*np.exp(-(eta[idx] - nB_eta_l)**2./(2.*sigma_eta**2.)))
        idx = (abs(eta - nB_eta_r) < etaDis)
        dNbdetas[idx] += (
            nB_frac_r*np.exp(-(eta[idx] - nB_eta_r)**2./(2.*sigma_eta**2.)))

        idx = (abs(eta - nB_y_l) < etaDis)
        dNbdy[idx] += (
            nB_frac_l*np.exp(-(eta[idx] - nB_y_l)**2./(2.*sigma_eta**2.)))
        idx = (abs(eta - nB_y_r) < etaDis)
        dNbdy[idx] += (
            nB_frac_r*np.exp(-(eta[idx] - nB_y_r)**2./(2.*sigma_eta**2.)))

dEdetas /= nev
dNbdetas *= norm_eta
dNbdetas /= nev
dNbdy *= norm_eta
dNbdy /= nev

output_ed = np.array([eta, dEdetas]).transpose()
output_nB = np.array([eta, dNbdetas, dNbdy]).transpose()

np.savetxt(path.join(work_folder, "dEdetas.dat"), output_ed,
           fmt='%.10e', delimiter='  ')
np.savetxt(path.join(work_folder, "dNBdis.dat"), output_nB,
           fmt='%.10e', delimiter='  ')
