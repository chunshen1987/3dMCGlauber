// Copyright @ Chun Shen 2022

#include <fstream>
#include <iostream>
#include "MakeDensity.h"

namespace MCGlb {

void MakeDensity::output_netBaryon_eta_distribution(std::string filename,
                                                    int eventId) {
    // compute the net baryon density profile
    std::vector<real> eta_arr(gridNx_, 0.);
    std::vector<float> nB_arr(gridNx_, 0.);
    for (int i = 0; i < gridNx_; i++) {
        eta_arr[i] = - gridSize_/2. + i*gridDx_;
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
    outFile.open(filename.c_str(), modes);
    for (int i = 0; i < gridNx_; i++) {
        outFile.write((char*) &(nB_arr[i]), sizeof(float));
    }
    outFile.close();
}

};
