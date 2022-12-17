// Copyright @ Chun Shen 2022

#include <fstream>
#include <iostream>
#include <sstream>
#include "MakeDensity.h"

namespace MCGlb {


void MakeDensity::output_energyDensity_xeta_distribution(
        std::string filename, const int eventId) const {
    // compute the 2D local energy density profile in (x, eta_s)
    std::vector<float> eta_arr(gridNx_, 0.);
    std::vector<float> x_arr(gridNy_, 0.);
    std::vector<float> ed_arr(gridNx_*gridNy_, 0.);

    std::vector<float> feC_arr(gridNx_, 0.);
    std::vector<float> feL_arr(gridNx_, 0.);
    std::vector<float> feR_arr(gridNx_, 0.);
    for (int i = 0; i < gridNx_; i++) {
        eta_arr[i] = - gridXSize_/2. + i*gridDx_;
    }
    for (int i = 0; i < gridNy_; i++) {
        x_arr[i] = - gridYSize_/2. + i*gridDy_;
    }

    double stringTransverseShiftFrac = 1.0;
    double sigma_eta = 0.5;
    double two_sigma_eta_sq = 2.*sigma_eta*sigma_eta;
    double sigmaDisEta = 5.*sigma_eta;
    double sigma_x = 0.2;
    double two_sigma_x_sq = 2.*sigma_x*sigma_x;
    double sigmaDisX = 5.*sigma_x;
    double normX = 1./(sqrt(2.*M_PI)*sigma_x);
    for (auto &string_i : QCD_string_output_arr_) {
        double mass = string_i[0];
        double xPerpC = string_i[5];
        double xPerpL = string_i[7];
        double xPerpR = string_i[9];
        double eta_l = string_i[11];
        double eta_r = string_i[12];
        double y_l = string_i[13];
        double y_r = string_i[14];
        double remFrac_l = string_i[15];
        double remFrac_r = string_i[16];
        double y_i_l = string_i[17];
        double y_i_r = string_i[18];

        double Estring = mass*(  cosh(y_i_l) + cosh(y_i_r)
                               - cosh(y_l) - cosh(y_r));
        double EremL = remFrac_l*mass*cosh(y_l);
        double EremR = remFrac_r*mass*cosh(y_r);

        double feCNorm = 0.;
        double feLNorm = 0.;
        double feRNorm = 0.;
        for (int i = 0; i < gridNx_; i++) {
            if (eta_arr[i] > eta_l && eta_arr[i] < eta_r) {
                double y_eta = (y_l + (y_r - y_l)/(eta_r - eta_l)
                                      *(eta_arr[i] - eta_l));
                feC_arr[i] = cosh(y_eta);
            } else if (eta_arr[i] > eta_r
                       && (eta_arr[i] - eta_r) < sigmaDisEta) {
                double dis = eta_arr[i] - eta_r;
                feC_arr[i] = exp(-dis*dis/two_sigma_eta_sq)*cosh(y_r);
            } else if (eta_arr[i] < eta_l
                       && (eta_l - eta_arr[i]) < sigmaDisEta) {
                double dis = std::abs(eta_arr[i] - eta_l);
                feC_arr[i] = exp(-dis*dis/two_sigma_eta_sq)*cosh(y_l);
            }
            double etaDis = std::abs(eta_arr[i] - eta_r);
            if (etaDis < sigmaDisEta) {
                feR_arr[i] = exp(-etaDis*etaDis/two_sigma_eta_sq)*cosh(y_r);
            }
            etaDis = std::abs(eta_arr[i] - eta_l);
            if (etaDis < sigmaDisEta) {
                feL_arr[i] = exp(-etaDis*etaDis/two_sigma_eta_sq)*cosh(y_l);
            }
            feCNorm += feC_arr[i];
            feLNorm += feL_arr[i];
            feRNorm += feR_arr[i];
        }

        for (int i = 0; i < gridNx_; i++) {
            double etaFrac = (eta_arr[i] - eta_l)/(eta_r - eta_l + 1e-16);
            etaFrac = std::max(0., std::min(1., etaFrac));
            for (int j = 0; j < gridNy_; j++) {
                int idx = i*gridNy_ + j;
                double xT = (xPerpC + stringTransverseShiftFrac
                                      *(0.5 - etaFrac)*(xPerpL - xPerpR)/2.);
                double xDis = std::abs(x_arr[j] - xT);
                double fx = 0.;
                if (xDis < sigmaDisX) {
                    fx = normX*exp(-xDis*xDis/two_sigma_x_sq);
                }
                ed_arr[idx] += (  Estring/(feCNorm*gridDx_)*feC_arr[i]
                                + EremL/(feLNorm*gridDx_)*feL_arr[i]
                                + EremR/(feRNorm*gridDx_)*feR_arr[i])*fx;
            }
        }
    }

    // output results
    std::ios_base::openmode modes;
    if (eventId == 0) {
        modes = std::ios::out | std::ios::binary;
    } else {
        modes = std::ios::app | std::ios::binary;
    }

    std::ofstream outFile;
    std::stringstream fileNameDressed;
    fileNameDressed << filename << "_Neta_" << gridNx_
                    << "_Nx_" << gridNy_ << ".dat";
    outFile.open(fileNameDressed.str().c_str(), modes);
    if (eventId == 0) {
        for (int i = 0; i < gridNx_; i++) {
            for (int j = 0; j < gridNy_; j++) {
                outFile.write((char*) &(eta_arr[i]), sizeof(float));
            }
        }
        for (int i = 0; i < gridNx_; i++) {
            for (int j = 0; j < gridNy_; j++) {
                outFile.write((char*) &(x_arr[j]), sizeof(float));
            }
        }
    }
    for (int i = 0; i < gridNx_*gridNy_; i++) {
        outFile.write((char*) &(ed_arr[i]), sizeof(float));
    }
    outFile.close();
}


void MakeDensity::output_energyDensity_eta_distribution(
        std::string filename, const int eventId) const {
    // compute the local energy density profile
    std::vector<float> eta_arr(gridNx_, 0.);
    std::vector<float> ed_arr(gridNx_, 0.);
    std::vector<float> feC_arr(gridNx_, 0.);
    std::vector<float> feL_arr(gridNx_, 0.);
    std::vector<float> feR_arr(gridNx_, 0.);
    for (int i = 0; i < gridNx_; i++) {
        eta_arr[i] = - gridXSize_/2. + i*gridDx_;
    }

    double sigma_eta = 0.2;
    double two_sigma_eta_sq = 2.*sigma_eta*sigma_eta;
    double sigmaDis = 5.*sigma_eta;
    for (auto &string_i : QCD_string_output_arr_) {
        double mass = string_i[0];
        double eta_l = string_i[11];
        double eta_r = string_i[12];
        double y_l = string_i[13];
        double y_r = string_i[14];
        double remFrac_l = string_i[15];
        double remFrac_r = string_i[16];
        double y_i_l = string_i[17];
        double y_i_r = string_i[18];
        double Estring = mass*(  cosh(y_i_l) + cosh(y_i_r)
                               - cosh(y_l) - cosh(y_r));
        double EremL = remFrac_l*mass*cosh(y_l);
        double EremR = remFrac_r*mass*cosh(y_r);

        double feCNorm = 0.;
        double feLNorm = 0.;
        double feRNorm = 0.;
        for (int i = 0; i < gridNx_; i++) {
            if (eta_arr[i] > eta_l && eta_arr[i] < eta_r) {
                double y_eta = (y_l + (y_r - y_l)/(eta_r - eta_l)
                                      *(eta_arr[i] - eta_l));
                feC_arr[i] = cosh(y_eta);
            } else if (eta_arr[i] > eta_r && (eta_arr[i] - eta_r) < sigmaDis) {
                double dis = eta_arr[i] - eta_r;
                feC_arr[i] = exp(-dis*dis/two_sigma_eta_sq)*cosh(y_r);
            } else if (eta_arr[i] < eta_l && (eta_l - eta_arr[i]) < sigmaDis) {
                double dis = std::abs(eta_arr[i] - eta_l);
                feC_arr[i] = exp(-dis*dis/two_sigma_eta_sq)*cosh(y_l);
            }
            double etaDis = std::abs(eta_arr[i] - eta_r);
            if (etaDis < sigmaDis) {
                feR_arr[i] = exp(-etaDis*etaDis/two_sigma_eta_sq)*cosh(y_r);
            }
            etaDis = std::abs(eta_arr[i] - eta_l);
            if (etaDis < sigmaDis) {
                feL_arr[i] = exp(-etaDis*etaDis/two_sigma_eta_sq)*cosh(y_l);
            }
            feCNorm += feC_arr[i];
            feLNorm += feL_arr[i];
            feRNorm += feR_arr[i];
        }

        for (int i = 0; i < gridNx_; i++) {
            ed_arr[i] += (  Estring/(feCNorm*gridDx_)*feC_arr[i]
                          + EremL/(feLNorm*gridDx_)*feL_arr[i]
                          + EremR/(feRNorm*gridDx_)*feR_arr[i]);
        }
    }

    // output results
    std::ios_base::openmode modes;
    if (eventId == 0) {
        modes = std::ios::out | std::ios::binary;
    } else {
        modes = std::ios::app | std::ios::binary;
    }

    std::ofstream outFile;
    std::stringstream fileNameDressed;
    fileNameDressed << filename << "_N_" << gridNx_ << ".dat";
    outFile.open(fileNameDressed.str().c_str(), modes);
    if (eventId == 0) {
        for (int i = 0; i < gridNx_; i++) {
            outFile.write((char*) &(eta_arr[i]), sizeof(float));
        }
    }
    for (int i = 0; i < gridNx_; i++) {
        outFile.write((char*) &(ed_arr[i]), sizeof(float));
    }
    outFile.close();
}


void MakeDensity::output_netBaryon_eta_distribution(std::string filename,
                                                    const int eventId) const {
    // compute the net baryon density profile
    std::vector<float> eta_arr(gridNx_, 0.);
    std::vector<float> nB_arr(gridNx_, 0.);
    for (int i = 0; i < gridNx_; i++) {
        eta_arr[i] = - gridXSize_/2. + i*gridDx_;
    }

    double sigma_eta = 0.2;
    double two_sigma_eta_sq = 2.*sigma_eta*sigma_eta;
    double sigmaDis = 5.*sigma_eta;
    double norm_eta = 1./(sqrt(2.*M_PI)*sigma_eta);
    for (auto &string_i : QCD_string_output_arr_) {
        double nB_eta_l = string_i[19];
        double nB_eta_r = string_i[20];
        double nB_frac_l = string_i[23];
        double nB_frac_r = string_i[24];
        for (int i = 0; i < gridNx_; i++) {
            double dis = std::abs(eta_arr[i] - nB_eta_l);
            if (dis < sigmaDis) {
                nB_arr[i] += nB_frac_l*exp(-dis*dis/two_sigma_eta_sq);
            }
            dis = std::abs(eta_arr[i] - nB_eta_r);
            if (dis < sigmaDis) {
                nB_arr[i] += nB_frac_r*exp(-dis*dis/two_sigma_eta_sq);
            }
        }
    }
    for (int i = 0; i < gridNx_; i++) {
        nB_arr[i] *= norm_eta;
    }

    // output results
    std::ios_base::openmode modes;
    if (eventId == 0) {
        modes = std::ios::out | std::ios::binary;
    } else {
        modes = std::ios::app | std::ios::binary;
    }

    std::ofstream outFile;
    std::stringstream fileNameDressed;
    fileNameDressed << filename << "_N_" << gridNx_ << ".dat";
    outFile.open(fileNameDressed.str().c_str(), modes);
    if (eventId == 0) {
        for (int i = 0; i < gridNx_; i++) {
            outFile.write((char*) &(eta_arr[i]), sizeof(float));
        }
    }
    for (int i = 0; i < gridNx_; i++) {
        outFile.write((char*) &(nB_arr[i]), sizeof(float));
    }
    outFile.close();
}

};
