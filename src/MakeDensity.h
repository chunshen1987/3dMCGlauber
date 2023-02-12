// Copyright @ Chun Shen 2022

#ifndef SRC_MAKEDENSITY_H_
#define SRC_MAKEDENSITY_H_

#include "Glauber.h"
#include "pretty_ostream.h"

#include<vector>
#include<string>

namespace MCGlb {

class MakeDensity {

 private:
    std::vector<std::vector<real>> QCD_string_output_arr_;
    std::vector<std::vector<real>> participantList_;
    pretty_ostream messager_;
    const int orderMax_ = 3;    // the maximum order of eccentricity
    int gridNx_, gridNy_, gridNeta_;
    double gridDx_, gridDy_, gridDeta_;
    double gridXSize_, gridYSize_, gridEtaSize_;
    double sigma_eta_, sigma_x_;
    double stringTransverseShiftFrac_;

 public:
    MakeDensity() = default;
    ~MakeDensity() {};

    void set_QCD_string_output_arr(
            std::vector<std::vector<real>> QCD_string_output_arr) {
        QCD_string_output_arr_ = QCD_string_output_arr;
    }


    void setParticipantList(std::vector<std::vector<real>> participantList) {
        participantList_ = participantList;
    }

    void setGaussianWidths(double sigma_x, double sigma_eta) {
        sigma_x_ = sigma_x;
        sigma_eta_ = sigma_eta;
    }

    void setStringTransShiftFrac(double stringTransverseShiftFrac) {
        stringTransverseShiftFrac_ = stringTransverseShiftFrac;
    }

    void set_1D_grid_info(int gridNx, double gridDx) {
        gridNx_ = gridNx;
        gridDx_ = gridDx;
        gridXSize_ = gridNx_*gridDx_;
    }

    void set_1D_grid_info_eta(int gridNeta, double gridDeta) {
        gridNeta_ = gridNeta;
        gridDeta_ = gridDeta;
        gridEtaSize_ = gridNeta_*gridDeta_;
    }

    void set_2D_grid_info(int gridNx, double gridDx,
                          int gridNy, double gridDy) {
        gridNx_ = gridNx;
        gridDx_ = gridDx;
        gridXSize_ = gridNx_*gridDx_;
        gridNy_ = gridNy;
        gridDy_ = gridDy;
        gridYSize_ = gridNy_*gridDy_;
    }

    void set_2D_grid_info_etax(int gridNeta, double gridDeta,
                               int gridNx, double gridDx) {
        gridNeta_ = gridNeta;
        gridDeta_ = gridDeta;
        gridEtaSize_ = gridNeta_*gridDeta_;
        gridNx_ = gridNx;
        gridDx_ = gridDx;
        gridXSize_ = gridNx_*gridDx_;
    }

    void set_3D_grid_info(int gridNx, double gridDx,
                          int gridNy, double gridDy,
                          int gridNeta, double gridDeta) {
        gridNx_ = gridNx;
        gridDx_ = gridDx;
        gridXSize_ = gridNx_*gridDx_;
        gridNy_ = gridNy;
        gridDy_ = gridDy;
        gridYSize_ = gridNy_*gridDy_;
        gridNeta_ = gridNeta;
        gridDeta_ = gridDeta;
        gridEtaSize_ = gridNeta_*gridDeta_;
    }

    int getIdx3D(int i, int j, int k) const {
        return((k*gridNx_ + i)*gridNy_ + j);
    }

    void output_netBaryon_eta_distribution(std::string filename,
                                           const int eventId) const;
    void output_energyDensity_eta_distribution(std::string filename,
                                               const int eventId) const;
    void output_energyDensity_xeta_distribution(std::string filename,
                                                const int eventId) const;
    void compute_energyDensity_3D_distribution(
        std::vector<float> &x_arr, std::vector<float> &y_arr,
        std::vector<float> &eta_arr, std::vector<float> &ed_arr) const;
    void computeTATB(
        std::vector<float> &x_arr, std::vector<float> &y_arr,
        std::vector<float> &TA_arr, std::vector<float> &TB_arr) const;

    void output_eccentricity(std::string filenameHeader,
                             const int eventId) const;
    void outputTATBEccentricity(std::string filenameHeader,
                                const int eventId) const;
};

};

#endif  // SRC_MAKEDENSITY_H_
