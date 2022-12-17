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
    pretty_ostream messager_;
    int gridNx_, gridNy_;
    double gridDx_, gridDy_;
    double gridXSize_, gridYSize_;

 public:
    MakeDensity() = default;
    ~MakeDensity() {};

    void set_QCD_string_output_arr(
            std::vector<std::vector<real>> QCD_string_output_arr) {
        QCD_string_output_arr_ = QCD_string_output_arr;
    }

    void set_1D_grid_info(int gridNx, double gridDx) {
        gridNx_ = gridNx;
        gridDx_ = gridDx;
        gridXSize_ = gridNx_*gridDx_;
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

    void output_netBaryon_eta_distribution(std::string filename,
                                           const int eventId) const;
    void output_energyDensity_eta_distribution(std::string filename,
                                               const int eventId) const;
    void output_energyDensity_xeta_distribution(std::string filename,
                                                const int eventId) const;

};

};

#endif  // SRC_MAKEDENSITY_H_
