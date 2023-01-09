// Copyright @ Chun Shen 2022

#include <fstream>
#include <iostream>
#include <sstream>
#include "MakeDensity.h"

namespace MCGlb {


void MakeDensity::compute_energyDensity_3D_distribution(
        std::vector<float> &x_arr, std::vector<float> &y_arr,
        std::vector<float> &eta_arr, std::vector<float> &ed_arr) const {
    // compute the 3D local energy density profile in (x, y, eta_s)
    x_arr.clear();
    y_arr.clear();
    eta_arr.clear();
    ed_arr.clear();
    x_arr.resize(gridNx_, 0.);
    y_arr.resize(gridNy_, 0.);
    eta_arr.resize(gridNeta_, 0.);
    ed_arr.resize(gridNx_*gridNy_*gridNeta_, 0.);

    std::vector<float> feC_arr(gridNeta_, 0.);
    std::vector<float> feL_arr(gridNeta_, 0.);
    std::vector<float> feR_arr(gridNeta_, 0.);
    for (int i = 0; i < gridNx_; i++) {
        x_arr[i] = - gridXSize_/2. + i*gridDx_;
    }
    for (int i = 0; i < gridNy_; i++) {
        y_arr[i] = - gridYSize_/2. + i*gridDy_;
    }
    for (int i = 0; i < gridNeta_; i++) {
        eta_arr[i] = - gridEtaSize_/2. + i*gridDeta_;
    }

    double stringTransverseShiftFrac = 1.0;
    double two_sigma_eta_sq = 2.*sigma_eta_*sigma_eta_;
    double sigmaDisEta = 5.*sigma_eta_;
    double two_sigma_x_sq = 2.*sigma_x_*sigma_x_;
    double sigmaDisXsq = 25.*sigma_x_*sigma_x_;
    double normX = 1./(2.*M_PI*sigma_x_*sigma_x_);
    for (auto &string_i : QCD_string_output_arr_) {
        double mass = string_i[0];
        double xPerpC = string_i[5];
        double yPerpC = string_i[6];
        double xPerpL = string_i[7];
        double yPerpL = string_i[8];
        double xPerpR = string_i[9];
        double yPerpR = string_i[10];
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
        for (int i = 0; i < gridNeta_; i++) {
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
        feCNorm *= gridDeta_;
        feRNorm *= gridDeta_;
        feLNorm *= gridDeta_;

        for (int k = 0; k < gridNeta_; k++) {
            double etaFrac = (eta_arr[k] - eta_l)/(eta_r - eta_l + 1e-16);
            etaFrac = std::max(0., std::min(1., etaFrac));
            double xT = (xPerpC + stringTransverseShiftFrac
                                  *(0.5 - etaFrac)*(xPerpL - xPerpR)/2.);
            double yT = (yPerpC + stringTransverseShiftFrac
                                  *(0.5 - etaFrac)*(yPerpL - yPerpR)/2.);
            for (int i = 0; i < gridNx_; i++) {
                for (int j = 0; j < gridNy_; j++) {
                    int idx = (k*gridNx_ + i)*gridNy_ + j;
                    double rDisSq = (  (x_arr[i] - xT)*(x_arr[i] - xT)
                                     + (y_arr[j] - yT)*(y_arr[j] - yT));
                    double fx = 0.;
                    if (rDisSq < sigmaDisXsq) {
                        fx = normX*exp(-rDisSq/two_sigma_x_sq);
                    }
                    ed_arr[idx] += (  Estring/feCNorm*feC_arr[k]
                                    + EremL/feLNorm*feL_arr[k]
                                    + EremR/feRNorm*feR_arr[k])*fx;
                }
            }
        }
    }
}


void MakeDensity::output_eccentricity(std::string filenameHeader,
                                      const int eventId) const {
    std::vector<float> x_arr, y_arr, eta_arr, ed_arr;
    compute_energyDensity_3D_distribution(x_arr, y_arr, eta_arr, ed_arr);
    std::vector<float> eccnReal(gridNeta_*orderMax_, 0.);
    std::vector<float> eccnImag(gridNeta_*orderMax_, 0.);
    std::vector<float> eccnNorm(gridNeta_*orderMax_, 0.);

    std::vector<int> rPow(orderMax_, 3);
    for (int i = 1; i < orderMax_; i++) {
        rPow[i] = i + 1;
    }

    for (int k = 0; k < gridNeta_; k++) {
        double x_o = 0.;
        double y_o = 0.;
        double norm = 0.;
        for (int i = 0; i < gridNx_; i++) {
            for (int j = 0; j < gridNy_; j++) {
                int idxEd = (k*gridNx_ + i)*gridNy_ + j;
                x_o += x_arr[i]*ed_arr[idxEd];
                y_o += y_arr[j]*ed_arr[idxEd];
                norm += ed_arr[idxEd];
            }
        }
        x_o /= norm;
        y_o /= norm;
        for (int i = 0; i < gridNx_; i++) {
            for (int j = 0; j < gridNy_; j++) {
                int idxEd = (k*gridNx_ + i)*gridNy_ + j;
                double rperp = sqrt(  (x_arr[i] - x_o)*(x_arr[i] - x_o)
                                    + (y_arr[j] - y_o)*(y_arr[j] - y_o));
                double phi = atan2(y_arr[j] - y_o, x_arr[i] - x_o);
                for (int ii = 0; ii < orderMax_; ii++) {
                    int idxEcc = ii*gridNeta_ + k;
                    int iorder = ii + 1;
                    double weight = pow(rperp, rPow[ii])*ed_arr[idxEd];
                    eccnNorm[idxEcc] += weight;
                    eccnReal[idxEcc] += weight*cos(iorder*phi);
                    eccnImag[idxEcc] += weight*sin(iorder*phi);
                }
            }
        }
    }

    for (unsigned idx = 0; idx < eccnNorm.size(); idx++) {
        eccnReal[idx] /= eccnNorm[idx];
        eccnImag[idx] /= eccnNorm[idx];
    }

    // output results
    std::ios_base::openmode modes;
    if (eventId == 0) {
        modes = std::ios::out | std::ios::binary;
    } else {
        modes = std::ios::app | std::ios::binary;
    }

    for (int iorder = 1; iorder <= orderMax_; iorder++) {
        std::ofstream outFile;
        std::stringstream fileNameDressed;
        fileNameDressed << filenameHeader << "_" << iorder
                        << "_Neta_" << gridNeta_ << ".dat";
        outFile.open(fileNameDressed.str().c_str(), modes);
        if (eventId == 0) {
            for (int i = 0; i < gridNeta_; i++) {
                outFile.write((char*) &(eta_arr[i]), sizeof(float));
            }
        }
        for (int i = 0; i < gridNeta_; i++) {
            int idx = (iorder - 1)*gridNeta_ + i;
            outFile.write((char*) &(eccnReal[idx]), sizeof(float));
        }
        for (int i = 0; i < gridNeta_; i++) {
            int idx = (iorder - 1)*gridNeta_ + i;
            outFile.write((char*) &(eccnImag[idx]), sizeof(float));
        }
        outFile.close();
    }
}


void MakeDensity::output_energyDensity_xeta_distribution(
        std::string filename, const int eventId) const {
    // compute the 2D local energy density profile in (x, eta_s)
    std::vector<float> eta_arr(gridNeta_, 0.);
    std::vector<float> x_arr(gridNx_, 0.);
    std::vector<float> ed_arr(gridNx_*gridNeta_, 0.);

    std::vector<float> feC_arr(gridNeta_, 0.);
    std::vector<float> feL_arr(gridNeta_, 0.);
    std::vector<float> feR_arr(gridNeta_, 0.);
    for (int i = 0; i < gridNeta_; i++) {
        eta_arr[i] = - gridEtaSize_/2. + i*gridDeta_;
    }
    for (int i = 0; i < gridNx_; i++) {
        x_arr[i] = - gridXSize_/2. + i*gridDx_;
    }

    double stringTransverseShiftFrac = 1.0;
    double two_sigma_eta_sq = 2.*sigma_eta_*sigma_eta_;
    double sigmaDisEta = 5.*sigma_eta_;
    double two_sigma_x_sq = 2.*sigma_x_*sigma_x_;
    double sigmaDisX = 5.*sigma_x_;
    double normX = 1./(sqrt(2.*M_PI)*sigma_x_);
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
        for (int i = 0; i < gridNeta_; i++) {
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
        feCNorm *= gridDeta_;
        feRNorm *= gridDeta_;
        feLNorm *= gridDeta_;

        for (int i = 0; i < gridNeta_; i++) {
            double etaFrac = (eta_arr[i] - eta_l)/(eta_r - eta_l + 1e-16);
            etaFrac = std::max(0., std::min(1., etaFrac));
            for (int j = 0; j < gridNx_; j++) {
                int idx = i*gridNx_ + j;
                double xT = (xPerpC + stringTransverseShiftFrac
                                      *(0.5 - etaFrac)*(xPerpL - xPerpR)/2.);
                double xDis = std::abs(x_arr[j] - xT);
                double fx = 0.;
                if (xDis < sigmaDisX) {
                    fx = normX*exp(-xDis*xDis/two_sigma_x_sq);
                }
                ed_arr[idx] += (  Estring/feCNorm*feC_arr[i]
                                + EremL/feLNorm*feL_arr[i]
                                + EremR/feRNorm*feR_arr[i])*fx;
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
    fileNameDressed << filename << "_Neta_" << gridNeta_
                    << "_Nx_" << gridNx_ << ".dat";
    outFile.open(fileNameDressed.str().c_str(), modes);
    if (eventId == 0) {
        for (int i = 0; i < gridNeta_; i++) {
            for (int j = 0; j < gridNx_; j++) {
                outFile.write((char*) &(eta_arr[i]), sizeof(float));
            }
        }
        for (int i = 0; i < gridNeta_; i++) {
            for (int j = 0; j < gridNx_; j++) {
                outFile.write((char*) &(x_arr[j]), sizeof(float));
            }
        }
    }
    for (int i = 0; i < gridNeta_*gridNx_; i++) {
        outFile.write((char*) &(ed_arr[i]), sizeof(float));
    }
    outFile.close();
}


void MakeDensity::output_energyDensity_eta_distribution(
        std::string filename, const int eventId) const {
    // compute the local energy density profile
    std::vector<float> eta_arr(gridNeta_, 0.);
    std::vector<float> ed_arr(gridNeta_, 0.);
    std::vector<float> feC_arr(gridNeta_, 0.);
    std::vector<float> feL_arr(gridNeta_, 0.);
    std::vector<float> feR_arr(gridNeta_, 0.);
    for (int i = 0; i < gridNeta_; i++) {
        eta_arr[i] = - gridEtaSize_/2. + i*gridDeta_;
    }

    double two_sigma_eta_sq = 2.*sigma_eta_*sigma_eta_;
    double sigmaDis = 5.*sigma_eta_;
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
        for (int i = 0; i < gridNeta_; i++) {
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

        for (int i = 0; i < gridNeta_; i++) {
            ed_arr[i] += (  Estring/(feCNorm*gridDeta_)*feC_arr[i]
                          + EremL/(feLNorm*gridDeta_)*feL_arr[i]
                          + EremR/(feRNorm*gridDeta_)*feR_arr[i]);
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
    fileNameDressed << filename << "_N_" << gridNeta_ << ".dat";
    outFile.open(fileNameDressed.str().c_str(), modes);
    if (eventId == 0) {
        for (int i = 0; i < gridNeta_; i++) {
            outFile.write((char*) &(eta_arr[i]), sizeof(float));
        }
    }
    for (int i = 0; i < gridNeta_; i++) {
        outFile.write((char*) &(ed_arr[i]), sizeof(float));
    }
    outFile.close();
}


void MakeDensity::output_netBaryon_eta_distribution(std::string filename,
                                                    const int eventId) const {
    // compute the net baryon density profile
    std::vector<float> eta_arr(gridNeta_, 0.);
    std::vector<float> nB_arr(gridNeta_, 0.);
    for (int i = 0; i < gridNeta_; i++) {
        eta_arr[i] = - gridEtaSize_/2. + i*gridDeta_;
    }

    double two_sigma_eta_sq = 2.*sigma_eta_*sigma_eta_;
    double sigmaDis = 5.*sigma_eta_;
    double norm_eta = 1./(sqrt(2.*M_PI)*sigma_eta_);
    for (auto &string_i : QCD_string_output_arr_) {
        double nB_eta_l = string_i[19];
        double nB_eta_r = string_i[20];
        double nB_frac_l = string_i[23];
        double nB_frac_r = string_i[24];
        for (int i = 0; i < gridNeta_; i++) {
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
    for (int i = 0; i < gridNeta_; i++) {
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
    fileNameDressed << filename << "_N_" << gridNeta_ << ".dat";
    outFile.open(fileNameDressed.str().c_str(), modes);
    if (eventId == 0) {
        for (int i = 0; i < gridNeta_; i++) {
            outFile.write((char*) &(eta_arr[i]), sizeof(float));
        }
    }
    for (int i = 0; i < gridNeta_; i++) {
        outFile.write((char*) &(nB_arr[i]), sizeof(float));
    }
    outFile.close();
}

};
