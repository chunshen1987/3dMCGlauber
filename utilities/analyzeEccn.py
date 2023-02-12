#!/usr/bin/env python3

import sys
import numpy as np
from os import path
import matplotlib.pyplot as plt

NCONFIG = 6

#%% load data
folder = str(sys.argv[1])
dEdeta = np.fromfile(path.join(folder, "ed_etas_distribution_N_72.dat"),
                     dtype="float32")
dEdeta = np.nan_to_num(dEdeta.reshape(-1, 72))
ecc23D = np.fromfile(path.join(folder, "ecc_ed_n_2_Neta_72.dat"),
                    dtype="float32")
ecc23D = np.nan_to_num(ecc23D.reshape(-1, 72))
ecc23D = ecc23D[1::2, :] + 1j*ecc23D[2::2, :]
ecc33D = np.fromfile(path.join(folder, "ecc_ed_n_3_Neta_72.dat"),
                    dtype="float32")
ecc33D = np.nan_to_num(ecc33D.reshape(-1, 72))
ecc33D = ecc33D[1::2, :] + 1j*ecc33D[2::2, :]
ecc22D = np.fromfile(
        path.join(folder, "ecc_ed_n_2_TATB_Nconfig_{}.dat".format(NCONFIG)),
                     dtype="float32")
ecc22D = np.nan_to_num(ecc22D.reshape(-1, NCONFIG))
ecc22D = ecc22D[::2, :] + 1j*ecc22D[1::2, :]
ecc32D = np.fromfile(
        path.join(folder, "ecc_ed_n_3_TATB_Nconfig_{}.dat".format(NCONFIG)),
                     dtype="float32")
ecc32D = np.nan_to_num(ecc32D.reshape(-1, NCONFIG))
ecc32D = ecc32D[::2, :] + 1j*ecc32D[1::2, :]

#%% perform centrality cut

etaArr = dEdeta[0, :]
etaIdx = (etaArr > -0.5) & (etaArr < 0.5)
Emid = np.sum(dEdeta[1:, etaIdx], axis=1)
EmidSorted = -np.sort(-Emid)
Totev = len(EmidSorted)

etaIdx = 36; print("eta_s = {:.1f}".format(etaArr[etaIdx]))
cenList = np.linspace(0., 100, 11)
cenMid = []
ecc2Res = []
ecc2Err = []
ecc3Res = []
ecc3Err = []
PearsonEcc2 = []
PearsonEcc3 = []
for icen in range(len(cenList) - 1):
    if cenList[icen] > cenList[icen + 1]: continue
    cenMid.append((cenList[icen] + cenList[icen+1])/2.)
    EmidHighCut = EmidSorted[int(cenList[icen]/100.*Totev)]
    EmidLowCut = EmidSorted[int(cenList[icen+1]/100.*Totev) - 1]
    cenIdx = (Emid > EmidLowCut) & (Emid <= EmidHighCut)
    nev = len(ecc23D[cenIdx, etaIdx])
    ecc23D_rms = np.sqrt(np.mean(np.abs(ecc23D[cenIdx, etaIdx])**2.))
    ecc23D_err = (np.std(np.abs(ecc23D[cenIdx, etaIdx])**2.)/np.sqrt(nev)
                  /(2.*ecc23D_rms))
    ecc33D_rms = np.sqrt(np.mean(np.abs(ecc33D[cenIdx, etaIdx])**2.))
    ecc33D_err = (np.std(np.abs(ecc33D[cenIdx, etaIdx])**2.)/np.sqrt(nev)
                  /(2.*ecc33D_rms))
    ecc22D_rms = np.sqrt(np.mean(np.abs(ecc22D[cenIdx, :])**2., axis=0))
    ecc22D_err = (np.std(np.abs(ecc22D[cenIdx, :])**2., axis=0)/np.sqrt(nev)
                  /(2.*ecc22D_rms))
    ecc32D_rms = np.sqrt(np.mean(np.abs(ecc32D[cenIdx, :])**2., axis=0))
    ecc32D_err = (np.std(np.abs(ecc32D[cenIdx, :])**2., axis=0)/np.sqrt(nev)
                  /(2.*ecc32D_rms))
    ecc2Res.append([ecc23D_rms] + list(ecc22D_rms))
    ecc2Err.append([ecc23D_err] + list(ecc22D_err))
    ecc3Res.append([ecc33D_rms] + list(ecc32D_rms))
    ecc3Err.append([ecc33D_err] + list(ecc32D_err))

    eccn3D = ecc23D[cenIdx, etaIdx].reshape(nev, 1)
    Pecc2 = (np.mean(np.real(eccn3D*np.conj(ecc22D[cenIdx, :])), axis=0)
             /ecc23D_rms/ecc22D_rms)
    eccn3D = ecc33D[cenIdx, etaIdx].reshape(nev, 1)
    Pecc3 = (np.mean(np.real(eccn3D*np.conj(ecc32D[cenIdx, :])), axis=0)
             /ecc33D_rms/ecc32D_rms)
    PearsonEcc2.append(Pecc2)
    PearsonEcc3.append(Pecc3)

ecc2Res = np.array(ecc2Res)
ecc2Err = np.array(ecc2Err)
ecc3Res = np.array(ecc3Res)
ecc3Err = np.array(ecc3Err)
PearsonEcc2 = np.array(PearsonEcc2)
PearsonEcc3 = np.array(PearsonEcc3)

#%% Make a plot
labelList = ["3D", "$T_A + T_B$", "$\sqrt{T_A T_B}$", "$(T_A T_B)^{2/3}$",
             "$T_A T_B$", "$T_A$", "$T_B$"]

fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.85, 0.85])
ax.errorbar(cenMid, ecc2Res[:, 0], ecc2Err[:, 0],
            marker='s', color='k', linestyle='', label=labelList[0])
for i in range(1, NCONFIG+1):
    ax.errorbar(cenMid, ecc2Res[:, i], ecc2Err[:, i], label=labelList[i])
plt.legend()
plt.xlabel("Centrality (%)")
plt.ylabel(r"$\epsilon_2\{2\}$")
plt.xlim([0, 100])
plt.ylim([0., 0.8])
plt.text(80, 0.05, r"$\eta_s = {:.1f}$".format(etaArr[etaIdx]), fontsize=15)
plt.savefig("ecc2_etas_{:.1f}.pdf".format(etaArr[etaIdx]))

fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.85, 0.85])
ax.errorbar(cenMid, ecc3Res[:, 0], ecc3Err[:, 0],
            marker='o', color='k', linestyle='', label=labelList[0])
for i in range(1, NCONFIG+1):
    ax.errorbar(cenMid, ecc3Res[:, i], ecc3Err[:, i], label=labelList[i])
plt.legend()
plt.xlabel("Centrality (%)")
plt.ylabel(r"$\epsilon_3\{2\}$")
plt.xlim([0, 100])
plt.ylim([0., 0.6])
plt.text(80, 0.05, r"$\eta_s = {:.1f}$".format(etaArr[etaIdx]), fontsize=15)
plt.savefig("ecc3_etas_{:.1f}.pdf".format(etaArr[etaIdx]))

fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.85, 0.85])
ax.plot(cenMid, ecc2Res[:, 0]/ecc3Res[:, 0],
        marker='o', color='k', linestyle='', label=labelList[0])
for i in range(1, NCONFIG+1):
    ax.plot(cenMid, ecc2Res[:, i]/ecc3Res[:, i], label=labelList[i])
plt.legend()
plt.xlabel("Centrality (%)")
plt.ylabel(r"$\epsilon_2\{2\}/\epsilon_3\{2\}$")
plt.xlim([0, 100])
plt.ylim([1., 3])
plt.text(80, 1.05, r"$\eta_s = {:.1f}$".format(etaArr[etaIdx]), fontsize=15)
plt.savefig("ecc2overecc3_etas_{:.1f}.pdf".format(etaArr[etaIdx]))


fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.85, 0.85])
for i in range(0, NCONFIG):
    ax.plot(cenMid, PearsonEcc2[:, i], label=labelList[i+1])
plt.legend()
plt.xlabel("Centrality (%)")
plt.ylabel(r"$P_2$")
plt.xlim([0, 100])
plt.ylim([0., 1.0])
plt.text(80, 0.05, r"$\eta_s = {:.1f}$".format(etaArr[etaIdx]), fontsize=15)
plt.savefig("Pearson_ecc2_etas_{:.1f}.pdf".format(etaArr[etaIdx]))

fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.85, 0.85])
for i in range(0, NCONFIG):
    ax.plot(cenMid, PearsonEcc3[:, i], label=labelList[i+1])
plt.legend()
plt.xlabel("Centrality (%)")
plt.ylabel(r"$P_3$")
plt.xlim([0, 100])
plt.ylim([0., 1.0])
plt.text(80, 0.05, r"$\eta_s = {:.1f}$".format(etaArr[etaIdx]), fontsize=15)
plt.savefig("Pearson_ecc3_etas_{:.1f}.pdf".format(etaArr[etaIdx]))
