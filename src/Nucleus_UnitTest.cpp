// Copyright (C) 2018 Chun Shen
#include "doctest.h"
#include "Nucleus.h"
#include <fstream>
#include <iostream>


TEST_CASE("Test random seed") {
    MCGlb::Nucleus test_nucleus;
    test_nucleus.set_random_seed(1);
    CHECK(test_nucleus.get_random_seed() == 1);
}

TEST_CASE("Test set nucleus parameters") {
    MCGlb::Nucleus test_nucleus;
    test_nucleus.set_nucleus_parameters("p");
    CHECK(test_nucleus.get_nucleus_A() == 1);
    CHECK(test_nucleus.get_nucleus_Z() == 1);
    auto WS_params = test_nucleus.get_woods_saxon_parameters();
    MCGlb::WoodsSaxonParam WS_params_p = {0.17, 0.0, 1.0, 1.0};
    CHECK(WS_params == WS_params_p);
    
    test_nucleus.set_nucleus_parameters("Au");
    CHECK(test_nucleus.get_nucleus_A() == 197);
    CHECK(test_nucleus.get_nucleus_Z() == 79);
    WS_params = test_nucleus.get_woods_saxon_parameters();
    MCGlb::WoodsSaxonParam WS_params_Au = {
                                0.17, 0.0, 6.38, 0.505, -0.13, -0.03};
    CHECK(WS_params == WS_params_Au);
}

TEST_CASE("Test generate nucleus configuratin") {
    MCGlb::Nucleus test_nucleus("p", 1);
    test_nucleus.generate_nucleus_3d_configuration();
    CHECK(test_nucleus.get_number_of_nucleons() == 1);
    auto test_nucleon = test_nucleus.get_nucleon(0);
    auto nucleon_x = test_nucleon.get_x();
    MCGlb::SpatialVec x_check = {0.0, 0.0, 0.0, 0.0};
    CHECK(nucleon_x == x_check);
    
    test_nucleus.set_nucleus_parameters("d");
    test_nucleus.generate_nucleus_3d_configuration();
    CHECK(test_nucleus.get_number_of_nucleons() == 2);
    auto nucleon1 = test_nucleus.get_nucleon(0);
    auto nucleon2 = test_nucleus.get_nucleon(1);
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
}

TEST_CASE("Test shift the nucleus") {
    MCGlb::Nucleus test_nucleus("p");
    test_nucleus.generate_nucleus_3d_configuration();
    MCGlb::SpatialVec x_shift = {0.0, 1.0, 0.0, -1.0};
    test_nucleus.shift_nucleus(x_shift);
    CHECK(test_nucleus.get_nucleon(0).get_x() == x_shift);
}

TEST_CASE("Test recenter the nucleus") {
    MCGlb::Nucleus test_nucleus("Au");
    test_nucleus.recenter_nucleus();
    auto nucleon_list = test_nucleus.get_nucleon_list();
    MCGlb::real meanx = 0., meany = 0., meanz = 0.;
    for (auto const &nucleon_i : nucleon_list) {
        auto x_vec = nucleon_i.get_x();
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
    MCGlb::Nucleus test_nucleus("Au");
    auto WS_params = test_nucleus.get_woods_saxon_parameters();
    auto a_WS = WS_params[3];
    auto R_WS = WS_params[2];

    const MCGlb::real r_min = 0.0, r_max = 20.0, dr = 0.1;
    const int n_r = static_cast<int>((r_max - r_min)/dr) + 1;
    std::vector<MCGlb::real> r(n_r, 0.);
    std::vector<MCGlb::real> rho_r(n_r, 0.);
    for (int i = 0; i < n_r; i++)
        r[i] = r_min + i*dr;

    int n_samples = 100000;
    auto weight   = 1./(n_samples*dr);
    for (int i = 0; i < n_samples; i++) {
        auto r_sample = test_nucleus.sample_r_from_woods_saxon();
        int idx       = static_cast<int>((r_sample - r_min)/dr);
        if (idx >= 0 && idx < n_r) {
            rho_r[idx] += weight;
        }
    }
    std::ofstream of("check_Woods_Saxon_sampling.dat");
    of << "# r  WS  Sampled" << std::endl;
    for (int i = 0; i < n_r; i++) {
        auto WS = r[i]*r[i]/(exp((r[i] - R_WS)/a_WS) + 1.)/92.;
        of << r[i] << "   " << WS << "  " << rho_r[i] << std::endl;
    }
    std::cout << "please check the output file check_Woods_Saxon_sampling.dat"
              << std::endl;
}

TEST_CASE("Test deformed nucleus") {
    MCGlb::Nucleus test_nucleus("U", -1, 0.9, true);
    test_nucleus.generate_nucleus_3d_configuration();
}

