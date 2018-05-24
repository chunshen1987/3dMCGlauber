// Copyright (C) 2018 Chun Shen
#include "doctest.h"
#include "Nucleus.h"
#include <fstream>
#include <memory>
#include <vector>
#include <iostream>

using MCGlb::Nucleus;
using MCGlb::real;
using MCGlb::SpatialVec;
using MCGlb::MomentumVec;
using MCGlb::WoodsSaxonParam;

TEST_CASE("Test random seed") {
    int seed = 23;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
                                            new RandomUtil::Random(seed));
    Nucleus test_nucleus("Au", ran_gen_ptr);
    CHECK(test_nucleus.get_random_seed() == seed);
}

TEST_CASE("Test set nucleus parameters") {
    Nucleus test_nucleus;
    test_nucleus.set_nucleus_parameters("p");
    CHECK(test_nucleus.get_nucleus_A() == 1);
    CHECK(test_nucleus.get_nucleus_Z() == 1);
    auto WS_params = test_nucleus.get_woods_saxon_parameters();
    WoodsSaxonParam WS_params_p = {0.17, 0.0, 1.0, 1.0, 0.0, 0.0};
    CHECK(WS_params == WS_params_p);
    
    test_nucleus.set_nucleus_parameters("Au");
    CHECK(test_nucleus.get_nucleus_A() == 197);
    CHECK(test_nucleus.get_nucleus_Z() == 79);
    WS_params = test_nucleus.get_woods_saxon_parameters();
    WoodsSaxonParam WS_params_Au = {
                                0.17, 0.0, 6.38, 0.505, -0.13, -0.03};
    CHECK(WS_params == WS_params_Au);
}

TEST_CASE("Test generate nucleus configuratin") {
    Nucleus test_nucleus("p");
    test_nucleus.generate_nucleus_3d_configuration();
    CHECK(test_nucleus.get_number_of_nucleons() == 1);
    auto test_nucleon = (*test_nucleus.get_nucleon(0));
    auto nucleon_x = test_nucleon.get_x();
    SpatialVec x_check = {0.0, 0.0, 0.0, 0.0};
    CHECK(nucleon_x == x_check);
    
    test_nucleus.set_nucleus_parameters("d");
    test_nucleus.generate_nucleus_3d_configuration();
    CHECK(test_nucleus.get_number_of_nucleons() == 2);
    auto nucleon1 = *test_nucleus.get_nucleon(0);
    auto nucleon2 = *test_nucleus.get_nucleon(1);
    auto nucleon1_x = nucleon1.get_x();
    auto nucleon2_x = nucleon2.get_x();
    CHECK(nucleon1_x[0] ==  nucleon2_x[0]);
    CHECK(nucleon1_x[1] == -nucleon2_x[1]);
    CHECK(nucleon1_x[2] == -nucleon2_x[2]);
    CHECK(nucleon1_x[3] == -nucleon2_x[3]);
    
    test_nucleus.set_nucleus_parameters("Au");
    test_nucleus.generate_nucleus_3d_configuration();
    CHECK(test_nucleus.get_number_of_nucleons()
            == test_nucleus.get_nucleus_A());
    CHECK(test_nucleus.get_nucleon_minimum_distance() == 0.9);
    test_nucleus.set_dmin(1.5);
    CHECK(test_nucleus.get_nucleon_minimum_distance() == 1.5);
}

TEST_CASE("Test shift the nucleus") {
    Nucleus test_nucleus("p");
    test_nucleus.generate_nucleus_3d_configuration();
    SpatialVec x_shift = {0.0, 1.0, 0.0, -1.0};
    test_nucleus.shift_nucleus(x_shift);
    CHECK(test_nucleus.get_nucleon(0)->get_x() == x_shift);
}

TEST_CASE("Test recenter the nucleus") {
    Nucleus test_nucleus("Au");
    test_nucleus.recenter_nucleus();
    auto nucleon_list = test_nucleus.get_nucleon_list();
    real meanx = 0., meany = 0., meanz = 0.;
    for (auto const &nucleon_i : (*nucleon_list)) {
        auto x_vec = nucleon_i->get_x();
        meanx += x_vec[1];
        meany += x_vec[2];
        meanz += x_vec[3];
    }
    CHECK(meanx == 0.0);
    CHECK(meany == 0.0);
    CHECK(meanz == 0.0);
}

TEST_CASE("Test Woods-Saxon sampling") {
    std::cout << "Testing the Woods-Saxon sampling routine..." << std::endl;
    Nucleus test_nucleus("Au");
    auto WS_params = test_nucleus.get_woods_saxon_parameters();
    auto a_WS = WS_params[3];
    auto R_WS = WS_params[2];

    const real r_min = 0.0, r_max = 20.0, dr = 0.1;
    const int n_r = static_cast<int>((r_max - r_min)/dr) + 1;
    std::vector<real> r(n_r, 0.);
    std::vector<real> rho_r(n_r, 0.);
    std::vector<real> WS(n_r, 0.);
    real norm_WS = 0.;
    for (int i = 0; i < n_r; i++) {
        real r_local  = r_min + i*dr;
        norm_WS += r_local*r_local/(exp((r_local - R_WS)/a_WS) + 1.)*dr;
    }

    int n_samples = 1000000;
    auto weight   = 1./(n_samples*dr);
    for (int i = 0; i < n_samples; i++) {
        auto r_sample = test_nucleus.sample_r_from_woods_saxon();
        int idx       = static_cast<int>((r_sample - r_min)/dr);
        if (idx >= 0 && idx < n_r) {
            r[idx]     += weight*r_sample;
            rho_r[idx] += weight;
        }
    }
    std::ofstream of("check_Woods_Saxon_sampling.dat");
    of << "# r  WS  Sampled" << std::endl;
    for (int i = 0; i < n_r; i++) {
        real r_mean = r_min + i*dr;
        if (rho_r[i] > 0) {
            r_mean = r[i]/rho_r[i];
        }
        real WS_local = r_mean*r_mean/(exp((r_mean - R_WS)/a_WS) + 1.)/norm_WS;
        of << r_mean << "   " << WS_local << "  " << rho_r[i] << std::endl;
    }
    std::cout << "please check the output file check_Woods_Saxon_sampling.dat"
              << std::endl;
}


TEST_CASE("Test deformed nucleus") {
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
                                            new RandomUtil::Random(-1));
    Nucleus test_nucleus("U", ran_gen_ptr, false, 0.9, true);
    test_nucleus.generate_nucleus_3d_configuration();
    CHECK(test_nucleus.get_number_of_nucleons() == 238);
    CHECK(test_nucleus.is_deformed() == true);
    
    Nucleus test_nucleus1("Au");
    CHECK(test_nucleus1.is_deformed() == true);
    Nucleus test_nucleus2("Au", nullptr, false, 0.9, false);
    CHECK(test_nucleus2.is_deformed() == false);
}


TEST_CASE("Test sampled nuclear density distribution") {
    std::cout << "Testing the sampling routine..." << std::endl;
    //Nucleus test_nucleus("Pb");
    Nucleus test_nucleus("Au", nullptr, false, 0.9, false);
    auto WS_params = test_nucleus.get_woods_saxon_parameters();
    auto a_WS = WS_params[3];
    auto R_WS = WS_params[2];

    const real r_min = 0.0, r_max = 20.0, dr = 0.1;
    const int n_r = static_cast<int>((r_max - r_min)/dr) + 1;
    std::vector<real> r(n_r, 0.);
    std::vector<real> rho_r(n_r, 0.);
    std::vector<real> WS(n_r, 0.);
    real norm_WS = 0.;
    for (int i = 0; i < n_r; i++) {
        real r_local  = r_min + i*dr;
        norm_WS += r_local*r_local/(exp((r_local - R_WS)/a_WS) + 1.)*dr;
    }

    int n_samples = 10000;
    auto weight   = 1./(n_samples*dr*test_nucleus.get_nucleus_A());
    for (int i = 0; i < n_samples; i++) {
        test_nucleus.generate_nucleus_3d_configuration();
        auto nucleon_list = test_nucleus.get_nucleon_list();
        for (auto const& it: (*nucleon_list)) {
            auto xvec = it->get_x();
            auto r_sample = sqrt(xvec[1]*xvec[1] + xvec[2]*xvec[2]
                                 + xvec[3]*xvec[3]);
            int idx       = static_cast<int>((r_sample - r_min)/dr);
            if (idx >= 0 && idx < n_r) {
                r[idx]     += weight*r_sample;
                rho_r[idx] += weight;
            }
        }
    }
    std::ofstream of("check_sampled_nucleon_distribution.dat");
    of << "# Nucleus: " << test_nucleus.get_name() << std::endl;
    of << "# r  WS  Sampled" << std::endl;

    for (int i = 0; i < n_r; i++) {
        real r_mean = r_min + i*dr;
        if (rho_r[i] > 0) {
            r_mean = r[i]/rho_r[i];
        }
        real WS_local = r_mean*r_mean/(exp((r_mean - R_WS)/a_WS) + 1.)/norm_WS;
        of << r_mean << "   " << WS_local << "  " << rho_r[i] << std::endl;
    }

    std::cout << "please check the output file "
              << "check_sampled_nucleon_distribution.dat" << std::endl;
}


TEST_CASE("Test get_z_max and get_z_min") {
    Nucleus test_nucleus("Au");
    test_nucleus.generate_nucleus_3d_configuration();
    test_nucleus.accelerate_nucleus(20., 1);
    auto z_max = test_nucleus.get_z_max();
    SpatialVec x_shift = {0.0, 0.0, 0.0, -z_max};
    test_nucleus.shift_nucleus(x_shift);
    for (auto const&it: (*test_nucleus.get_nucleon_list())) {
        auto xvec = it->get_x();
        CHECK(xvec[3] <= 0);
    }
    auto z_min = test_nucleus.get_z_min();
    x_shift = {0.0, 0.0, 0.0, -z_min};
    test_nucleus.shift_nucleus(x_shift);
    for (auto const&it: (*test_nucleus.get_nucleon_list())) {
        auto xvec = it->get_x();
        CHECK(xvec[3] >= 0);
    }
}

TEST_CASE("Test accelerate_nucleus()") {
    Nucleus test_nucleus1("Au");
    Nucleus test_nucleus2("Pb");
    test_nucleus1.generate_nucleus_3d_configuration();
    test_nucleus2.generate_nucleus_3d_configuration();
    test_nucleus1.accelerate_nucleus(20.,  1);
    test_nucleus2.accelerate_nucleus(20., -1);
    for (auto const&it: (*test_nucleus1.get_nucleon_list())) {
        auto pvec = it->get_p();
        CHECK(pvec[3]/pvec[0] > 0.);
    }
    for (auto const&it: (*test_nucleus2.get_nucleon_list())) {
        auto pvec = it->get_p();
        CHECK(pvec[3]/pvec[0] < 0.);
    }
}

TEST_CASE("Test get_number_of_wounded_nucleons()") {
    Nucleus test_nucleus1("Au");
    test_nucleus1.generate_nucleus_3d_configuration();
    auto list = test_nucleus1.get_nucleon_list();
    for (auto &it: (*list)) {
        it->set_wounded(true);
    }
    CHECK(test_nucleus1.get_number_of_wounded_nucleons()
            == test_nucleus1.get_nucleus_A());
    auto nucleon_i = test_nucleus1.get_nucleon(0);
    nucleon_i->set_wounded(false);
    CHECK(test_nucleus1.get_number_of_wounded_nucleons()
            == (test_nucleus1.get_nucleus_A() - 1));

}

TEST_CASE("Test rotate_nucleus") {
}


TEST_CASE("Test sampled valence quark spatial distribution") {
    std::cout << "Testing the valence quark sampling routine..." << std::endl;
    Nucleus test_nucleus("Au");
    const real r_min = 0.0, r_max = 4.0, dr = 0.05;
    const int n_r = static_cast<int>((r_max - r_min)/dr) + 1;
    const real a = 3.87;
    std::vector<real> r(n_r, 0.);
    std::vector<real> rho_r(n_r, 0.);
    std::vector<real> WS(n_r, 0.);
    real norm_WS = 0.;
    for (int i = 0; i < n_r; i++) {
        real r_local = r_min + (i + 0.5)*dr;
        norm_WS += test_nucleus.ExponentialDistribution(a, r_local)*dr;
    }

    int n_samples = 1000000;
    auto weight   = 1./(n_samples*dr);
    for (int i = 0; i < n_samples; i++) {
        auto xvec     = test_nucleus.sample_valence_quark_position();
        auto r_sample = sqrt(xvec[1]*xvec[1] + xvec[2]*xvec[2]
                             + xvec[3]*xvec[3]);
        int idx       = static_cast<int>((r_sample - r_min)/dr);
        if (idx >= 0 && idx < n_r) {
            r[idx]     += weight*r_sample;
            rho_r[idx] += weight;
        }
    }
    std::ofstream of("check_sampled_valence_quark_spatial_distribution.dat");
    of << "# r  ExponentialDistribution  Sampled" << std::endl;

    for (int i = 0; i < n_r; i++) {
        real r_mean = r[i]/rho_r[i];
        of << r_mean << "   "
           << test_nucleus.ExponentialDistribution(a, r_mean)/norm_WS << "  "
           << rho_r[i] << std::endl;
    }

    std::cout << "please check the output file "
              << "check_sampled_valence_quark_spatial_distribution.dat"
              << std::endl;
}


TEST_CASE("Test sample quark momentum fraction") {
    std::cout << "Testing the sampling routine for valence quarks..."
              << std::endl;
    Nucleus test_nucleus1("Au", nullptr, true);

    const real x_min = 0.0, x_max = 1.0, dx = 0.02;
    const int n_x = static_cast<int>((x_max - x_min)/dx);
    std::vector<real> x(n_x, 0.);
    std::vector<real> Px_0(n_x, 0.);
    std::vector<real> Px_1(n_x, 0.);
    std::vector<real> Px_2(n_x, 0.);
    for (int i = 0; i < n_x; i++) {
        x[i]  = x_min + i*dx;
    }

    int n_samples = 10000;
    for (int i = 0; i < n_samples; i++) {
        std::vector<real> xQuark;
        test_nucleus1.sample_quark_momentum_fraction(xQuark, 3);
        int x_idx = static_cast<int>((xQuark[0] - x_min)/dx);
        if (x_idx >= 0 && x_idx < n_x)
            Px_0[x_idx]++;
        x_idx = static_cast<int>((xQuark[1] - x_min)/dx);
        if (x_idx >= 0 && x_idx < n_x)
            Px_1[x_idx]++;
        x_idx = static_cast<int>((xQuark[2] - x_min)/dx);
        if (x_idx >= 0 && x_idx < n_x)
            Px_2[x_idx]++;
    }
    std::ofstream of("check_sampled_valence_quarks_distribution.dat");
    of << "# x  P_0(x)  P_1(x)  P_2(x)" << std::endl;

    for (int i = 0; i < n_x; i++) {
        of << x[i] << "   "
           << Px_0[i]/n_samples/dx << "  "
           << Px_1[i]/n_samples/dx << "  "
           << Px_2[i]/n_samples/dx << "  "
           << std::endl;
    }
    of.close();
    
    std::cout << "please check the output file "
              << "check_sampled_valence_quarks_distribution.dat" << std::endl;
}
