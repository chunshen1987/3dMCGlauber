// Copyright @ Chun Shen 2018

#include "Nucleus.h"
#include "PhysConsts.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <array>
#include <tuple>
#include <cmath>
#include <cstdlib>
#include <utility>

#include "eps09.h"

using std::cout;
using std::endl;

namespace MCGlb {

Nucleus::Nucleus(std::string nucleus_name,
                 std::shared_ptr<RandomUtil::Random> ran_gen,
                 bool sample_valence_quarks_in, real BG,
                 real d_min, bool deformed, bool confFromFile) {
    d_min_    = d_min;
    deformed_ = deformed;
    confFromFile_ = confFromFile;
    BG_ = BG;
    ran_gen_ptr = ran_gen;
    set_nucleus_parameters(nucleus_name);

    sample_valence_quarks = sample_valence_quarks_in;
    if (sample_valence_quarks) {
        number_of_valence_quark_samples_ = readin_valence_quark_samples();
    }
    nucleon_configuration_loaded_ = false;
}


Nucleus::~Nucleus() {
    participant_list_.clear();
    nucleon_list_.clear();
    if (sample_valence_quarks) {
        proton_valence_quark_x_.clear();
        neutron_valence_quark_x_.clear();
    }
}


void Nucleus::set_woods_saxon_parameters(int A_in, int Z_in,
                                         real rho, real w, real R, real a,
                                         real beta2, real beta3, real beta4,
                                         real gamma, real da, real dR,
                                         int density_function_type_in) {
    A_ = A_in;
    Z_ = Z_in;
    density_function_type = density_function_type_in;
    setWoodsSaxonParameters(rho, w, R, a, beta2, beta3, beta4, gamma, da, dR);
}


void Nucleus::setWoodsSaxonParameters(real rho, real w, real R, real a,
                                      real beta2, real beta3, real beta4,
                                      real gamma, real da, real dR) {
    WS_param_vec[0] = rho;
    WS_param_vec[1] = w;
    WS_param_vec[2] = R;
    WS_param_vec[3] = a;
    WS_param_vec[4] = beta2;
    WS_param_vec[5] = beta3;
    WS_param_vec[6] = beta4;
    WS_param_vec[7] = gamma;
    WS_param_vec[8] = da;
    WS_param_vec[9] = dR;

    if (gamma > 0. && d_min_ > 0.) {
        cout << "[Nucleus]: With a non-zero WS_gamma, we do not support "
             << "a non-zero d_min between nucleon pairs! "
             << "resetting d_min to 0." << endl;
        d_min_ = 0.;
    }
}


void Nucleus::set_nucleus_parameters(std::string nucleus_name) {
    name = nucleus_name;
    if (nucleus_name.compare("p") == 0) {
        set_woods_saxon_parameters(
                1, 1, 0.17, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 3);
    }else if (nucleus_name.compare("dipole") == 0) {
        set_woods_saxon_parameters(
                0, 0, 0.17, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 8);
    }else if (nucleus_name.compare("d") == 0) {
        set_woods_saxon_parameters(
                2, 1, 0.17, 1.18, 1.0, 0.228, 0.0, 0.0, 0.0, 0.0, 0, 0, 8);
    } else if (nucleus_name.compare("He3") == 0) {
        set_woods_saxon_parameters(
                3, 2, 0.17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 1);
    } else if (nucleus_name.compare("He4") == 0) {
        set_woods_saxon_parameters(
                4, 2, 0.17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1);
    } else if (nucleus_name.compare("C") == 0) {
        set_woods_saxon_parameters(
                12, 6, 0.17, 1.403, 2.44, 1.635, 0, 0, 0, 0, 0, 0, 1);
    } else if (nucleus_name.compare("O") == 0) {
        set_woods_saxon_parameters(
                16, 8, 0.17, -0.051, 2.608, 0.513, 0, 0, 0, 0, 0, 0, 3);
    } else if (nucleus_name.compare("Al") == 0) {
        set_woods_saxon_parameters(
                27, 13, 0.17, 0.0, 3.07, 0.519, 0.0, 0.0, 0.0, 0.0, 0, 0, 3);
    } else if (nucleus_name.compare("Cu") == 0) {
        set_woods_saxon_parameters(
                63, 29, 0.17, 0.0, 4.163, 0.606, 0.162, 0, 0.006, 0, 0, 0, 3);
    } else if (nucleus_name.compare("Zr") == 0) {
        set_woods_saxon_parameters(
                96, 40, 0.17, 0.0, 5.020, 0.520, 0.06, 0.16, 0, 0, 0, 0, 3);
    } else if (nucleus_name.compare("Ru") == 0) {
        set_woods_saxon_parameters(
                96, 44, 0.17, 0.0, 5.090, 0.460, 0.16, 0, 0, 0, 0, 0, 3);
    } else if (nucleus_name.compare("In") == 0) {
        set_woods_saxon_parameters(
                115, 49, 0.17, 0.0, 5.35, 0.55, 0, 0, 0, 0, 0, 0, 3);
    } else if (nucleus_name.compare("Xe") == 0) {
        set_woods_saxon_parameters(
                129, 54, 0.17, 0.0, 5.36, 0.590, 0.162, 0, -0.003, 0, 0, 0, 3);
    } else if (nucleus_name.compare("Au") == 0) {
        set_woods_saxon_parameters(
                197, 79, 0.17, 0.0, 6.38, 0.505, -0.13, 0, -0.03, 0, 0, 0, 3);
    } else if (nucleus_name.compare("Pb") == 0) {
        set_woods_saxon_parameters(
                208, 82, 0.17, 0.0, 6.62, 0.546, 0, 0, 0, 0, 0, 0, 3);
    } else if (nucleus_name.compare("U") == 0) {
        set_woods_saxon_parameters(
                238, 92, 0.17, 0.0, 6.874, 0.556, 0.28, 0, 0.093, 0, 0, 0, 3);
    } else {
        cout << "[Error] Unknown_nucleus: " << nucleus_name << endl;
        cout << "Exiting... " << endl;
        exit(1);
    }
}


void Nucleus::generate_nucleus_3d_configuration() {
    if (nucleon_list_.size() > 0)
        nucleon_list_.clear();

    if (participant_list_.size() > 0)
        participant_list_.clear();

    int status = 2;
    // sample the nucleons' positions
    if (A_ == 1) {  // p
        SpatialVec  x = {0.0};
        MomentumVec p = {0.0};
        std::shared_ptr<Nucleon> nucleon_ptr(new Nucleon(x, p));
        nucleon_ptr->set_electric_charge(1.);
        nucleon_list_.push_back(std::move(nucleon_ptr));
        status = 0;
    } else if (A_ == 2) {   // deuteron
        generate_deuteron_configuration();
        status = 0;
    } else if (confFromFile_) {
        status = sample_nucleon_configuration();
    }

    if (status != 0) {
        if (!deformed_) {
            generate_nucleus_configuration_with_woods_saxon();
        } else {
            generate_nucleus_configuration_with_deformed_woods_saxon();
        }
    }

    recenter_nucleus();

    sample_fermi_momentum();

    real phi   = 2.*M_PI*ran_gen_ptr->rand_uniform();
    real theta = acos(1. - 2.*ran_gen_ptr->rand_uniform());
    rotate_nucleus(phi, theta);
}


void Nucleus::recenter_nucleus() {
    // compute the center of mass position and shift it to (0, 0, 0)
    real meanx = 0., meany = 0., meanz = 0.;
    for (auto const &nucleon_i : nucleon_list_) {
        auto x_vec = nucleon_i->get_x();
        meanx += x_vec[1];
        meany += x_vec[2];
        meanz += x_vec[3];
    }

    meanx /= static_cast<real>(A_);
    meany /= static_cast<real>(A_);
    meanz /= static_cast<real>(A_);

    SpatialVec x_shift = {0, -meanx, -meany, -meanz};
    shift_nucleus(x_shift);
}


void Nucleus::shift_nucleus(SpatialVec x_shift) {
    for (auto &nucleon_i : nucleon_list_) {
        auto x_vec = nucleon_i->get_x();
        for (int i = 0; i < 4; i++)
            x_vec[i] += x_shift[i];
        nucleon_i->set_x(x_vec);
    }
}


void Nucleus::rotate_nucleus(real phi, real theta) {
    auto cth  = cos(theta);
    auto sth  = sin(theta);
    auto cphi = cos(phi);
    auto sphi = sin(phi);
    for (auto &nucleon_i : nucleon_list_) {
        auto x_vec = nucleon_i->get_x();
        auto x_new = cth*cphi*x_vec[1] - sphi*x_vec[2] + sth*cphi*x_vec[3];
        auto y_new = cth*sphi*x_vec[1] + cphi*x_vec[2] + sth*sphi*x_vec[3];
        auto z_new = -sth    *x_vec[1] + 0.  *x_vec[2] + cth     *x_vec[3];
        x_vec[1] = x_new; x_vec[2] = y_new; x_vec[3] = z_new;
        nucleon_i->set_x(x_vec);
    }
}


void Nucleus::sample_fermi_momentum() {
    const real pi2_3 = 3.*M_PI*M_PI;
    for (auto &nucleon_i : nucleon_list_) {
        auto x_vec = nucleon_i->get_x();
        real r = sqrt(x_vec[1]*x_vec[1] + x_vec[2]*x_vec[2]
                      + x_vec[3]*x_vec[3]);
        real avgDensity = getAvgWoodsSaxonDensity(r);
        real denFrac = static_cast<real>(Z_)/static_cast<real>(A_);
        if (std::abs(nucleon_i->get_electric_charge()) < 1e-8)
            denFrac = 1. - denFrac;
        real p_F = PhysConsts::HBARC*pow(pi2_3*avgDensity*denFrac
                                         *ran_gen_ptr->rand_uniform(), 1./3.);
        real phi = 2.*M_PI*ran_gen_ptr->rand_uniform();
        real theta = acos(1. - 2.*ran_gen_ptr->rand_uniform());
        nucleon_i->set_fermi_momentum(p_F*sin(theta)*cos(phi),
                                      p_F*sin(theta)*sin(phi),
                                      p_F*cos(theta));
    }
}


void Nucleus::sample_valence_quarks_inside_nucleons(real ecm, int direction) {
    const int number_of_quarks = PhysConsts::NumValenceQuark;
    for (auto &nucleon_i: nucleon_list_) {
        if (nucleon_i->is_wounded()
            && nucleon_i->get_number_of_quarks() == 0) {
            std::vector<real> xQuark;
            std::vector<real> eQuark;
            sample_quark_momentum_fraction(xQuark, eQuark, number_of_quarks,
                                           nucleon_i->get_electric_charge(),
                                           ecm);
            for (int i = 0; i < number_of_quarks; i++) {
                auto xvec = sample_valence_quark_position();
                std::shared_ptr<Quark> quark_ptr(new Quark(xvec, xQuark[i],
                                                           eQuark[i]));
                nucleon_i->push_back_quark(quark_ptr);
            }
            nucleon_i->accelerate_quarks(ecm, direction);
        }
    }
}


void Nucleus::add_soft_parton_ball(real ecm, int direction) {
    real E_rem_min = 0.1;  // GeV
    real beam_rapidity = direction*acosh(ecm/(2.*PhysConsts::MProton));
    real Pz_rem_min = E_rem_min*tanh(beam_rapidity);
    for (auto &nucleon_i: nucleon_list_) {
        if (nucleon_i->is_wounded()
            && nucleon_i->get_number_of_quarks() != 0) {
            auto soft_pvec = nucleon_i->get_p();
            auto valence_quark_list = nucleon_i->get_quark_list();
            for (const auto & q_i: valence_quark_list) {
                auto quark_pvec = q_i->get_p();
                for (int i = 0; i < 4; i++) {
                    soft_pvec[i] -= quark_pvec[i];
                }
            }
            soft_pvec[0] -= E_rem_min;
            soft_pvec[3] -= Pz_rem_min;
            real mass = PhysConsts::MQuarkValence;
            if (soft_pvec[0] > mass) {
                // assuming the soft parton ball has valence quark mass
                // only add a soft parton when the leftover energy is
                // larger than mq
                real rapidity = direction*acosh(soft_pvec[0]/mass);
                soft_pvec[3] = mass*sinh(rapidity);
                auto xvec = sample_valence_quark_position();
                std::shared_ptr<Quark> quark_ptr(new Quark(xvec, soft_pvec));
                quark_ptr->set_rapidity(rapidity);
                nucleon_i->push_back_quark(quark_ptr);
            }
        }
    }
}

void Nucleus::generate_deuteron_configuration() {
    // sample the distance between the two nucleons
    real r_dis = get_inverse_CDF_hulthen_function(ran_gen_ptr->rand_uniform());
    // sample the solid angles
    real phi = 2.*M_PI*ran_gen_ptr->rand_uniform();
    real theta = acos(1. - 2.*ran_gen_ptr->rand_uniform());

    // calculate the spatial poisition in the center of mass frame
    real x = 0.5*r_dis*sin(theta)*cos(phi);
    real y = 0.5*r_dis*sin(theta)*sin(phi);
    real z = 0.5*r_dis*cos(theta);
    SpatialVec  x_1 = {0.0,  x,  y,  z};
    SpatialVec  x_2 = {0.0, -x, -y, -z};
    MomentumVec p_1 = {0.0};
    MomentumVec p_2 = p_1;

    std::shared_ptr<Nucleon> nucleon1_ptr(new Nucleon(x_1, p_1));
    nucleon1_ptr->set_electric_charge(1.);
    nucleon_list_.push_back(std::move(nucleon1_ptr));
    std::shared_ptr<Nucleon> nucleon2_ptr(new Nucleon(x_2, p_2));
    nucleon2_ptr->set_electric_charge(0.);
    nucleon_list_.push_back(std::move(nucleon2_ptr));
}


real Nucleus::get_inverse_CDF_hulthen_function(real y) const {
    if (y < 0. || y > 1.) {
        cout << "[Error]: Glauber::get_inverse_CDF_hulthen_function: "
             << "input y < 0 or y > 1, y = " << y << endl;
        cout << "can not find the value from inverse CDF!" << endl;
        exit(1);
    }
    real x_min = 0.0;
    real x_max = 100.0;
    real y_max = hulthen_function_CDF(x_max);
    real abs_tol = 1e-8;
    while (y > y_max) {
        x_max += 100.0;
        y_max = hulthen_function_CDF(x_max);
    }
    real x_mid = (x_min + x_max)/2.;
    real y_mid = hulthen_function_CDF(x_mid);
    while (fabs(y - y_mid) > abs_tol) {
        if (y > y_mid) {
            x_min = x_mid;
        } else {
            x_max = x_mid;
        }
        x_mid = (x_min + x_max)/2.;
        y_mid = hulthen_function_CDF(x_mid);
    }
    return(x_mid);
}


real Nucleus::hulthen_function_CDF(real r) const {
    real alpha = 0.228;
    real beta  = 1.18;
    real c     = alpha*beta*(alpha + beta)/((alpha - beta)*(alpha - beta));
    real res = (2.*c*(2.*exp(-r*(alpha + beta))/(alpha + beta)
                 - 1./2.*exp(-2.*alpha*r)/alpha
                 - 1./2.*exp(-2.*beta*r)/beta
                 + 1./(2.*alpha) + 1./(2.*beta) - 2./(alpha + beta)));
    return(res);
}


int Nucleus::readin_valence_quark_samples() {
    std::stringstream of_p_name, of_n_name;
    of_p_name << "tables/proton_valence_quark_samples";
    of_n_name << "tables/neutron_valence_quark_samples";
    if (A_ == 197) {
        of_p_name << "_NPDFAu.dat";
        of_n_name << "_NPDFAu.dat";
    } else if (A_ == 208) {
        of_p_name << "_NPDFPb.dat";
        of_n_name << "_NPDFPb.dat";
    } else {
        of_p_name << ".dat";
        of_n_name << ".dat";
    }
    std::ifstream of_p_test(of_p_name.str().c_str(), std::ios::binary);
    std::ifstream of_n_test(of_n_name.str().c_str(), std::ios::binary);
    if (!of_p_test.good() || !of_n_test.good()) {
        std::cout << "Generating files " << of_p_name.str()
                  << " and " << of_n_name.str() << std::endl;
        std::stringstream command;
        command << "./Metropolis.e " << A_;
        system_status_ = std::system(command.str().c_str());
    } else {
        of_p_test.close();
        of_n_test.close();
    }

    std::ifstream of_p(of_p_name.str().c_str(), std::ios::binary);
    std::ifstream of_n(of_n_name.str().c_str(), std::ios::binary);
    if (!of_p.good() || !of_n.good()) {
        std::cout << "Can not generate " << of_p_name.str() << " or/and "
                  << of_n_name.str() << std::endl;
        std::cout << "Please check, exiting ... " << std::endl;
        exit(1);
    }
    while (!of_p.eof()) {
        std::array<float, 3> x_array;
        for (int ii = 0; ii < 3; ii++) {
            float temp = 0.;
            of_p.read(reinterpret_cast<char*>(&temp), sizeof(float));
            x_array[ii] = temp;
        }
        proton_valence_quark_x_.push_back(x_array);
    }
    of_p.close();
    while (!of_n.eof()) {
        std::array<float, 3> x_array;
        for (int ii = 0; ii < 3; ii++) {
            float temp = 0.;
            of_n.read(reinterpret_cast<char*>(&temp), sizeof(float));
            x_array[ii] = temp;
        }
        neutron_valence_quark_x_.push_back(x_array);
    }
    of_n.close();
    int size = std::min(proton_valence_quark_x_.size(),
                        neutron_valence_quark_x_.size());
    return(size);
}


void Nucleus::readin_triton_position() {
    // This function reads in spatial configuration for triton
    std::ifstream triton_position("tables/triton_positions.dat");
    if (!triton_position.good()) {
        std::cout << "Triton configurations are not found!" << std::endl;
        std::cout << "Please check file tables/triton_positions.dat."
                  << std::endl;
        exit(1);
    }
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
    triton_position >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
    while (!triton_position.eof()) {
        std::array<double, 9> temp;
        temp = {x1, y1, z1, x2, y2, z2, x3, y3, z3};
        triton_pos_.push_back(temp);
        triton_position >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
    }
    triton_position.close();
}


void Nucleus::generate_triton_configuration() {
    // This function samples the spatial configuration for triton
    if (triton_pos_.size() == 0)
        readin_triton_position();
    const int num_configuration = triton_pos_.size();
    const int rand_num = static_cast<int>(
                            ran_gen_ptr->rand_uniform()*num_configuration);
    auto pos_i = triton_pos_[rand_num];
    SpatialVec  x_1 = {0.0, pos_i[0], pos_i[1], pos_i[2]};
    SpatialVec  x_2 = {0.0, pos_i[3], pos_i[4], pos_i[5]};
    SpatialVec  x_3 = {0.0, pos_i[6], pos_i[7], pos_i[8]};
    MomentumVec p_1 = {0.0};
    MomentumVec p_2 = p_1;
    MomentumVec p_3 = p_1;

    std::shared_ptr<Nucleon> nucleon1_ptr(new Nucleon(x_1, p_1));
    nucleon_list_.push_back(std::move(nucleon1_ptr));
    std::shared_ptr<Nucleon> nucleon2_ptr(new Nucleon(x_2, p_2));
    nucleon_list_.push_back(std::move(nucleon2_ptr));
    std::shared_ptr<Nucleon> nucleon3_ptr(new Nucleon(x_3, p_3));
    nucleon_list_.push_back(std::move(nucleon3_ptr));
}

void Nucleus::readin_nucleon_positions() {
    std::cout << "read in nucleon positions for Nucleus: " << name << "  "
              << std::flush;
    std::ostringstream filename;
    int n_configuration = 0;
    if (A_ == 3) {  // he3
        filename << "tables/he3_plaintext.dat";
        n_configuration = 13699;
    } else if (A_ == 4) {  // he4
        filename << "tables/he4_plaintext.dat";
        n_configuration = 6000;
    } else if (A_ == 12) {  // carbon
        filename << "tables/carbon_plaintext.dat";
        n_configuration = 6000;
    } else if (A_ == 16) {  // oxygen
        filename << "tables/oxygen_plaintext.dat";
        n_configuration = 6000;
    } else if (A_ == 197) {  // Au
        filename << "tables/au197-sw-full_3Bchains-conf1820.dat";
        n_configuration = 1820;
    } else if (A_ == 208) {  // Pb
        int temp = 1;
        filename << "tables/pb208-" << temp << ".dat";
        n_configuration = 10000;
    } else {
        std::cout << "[Warning]: No configuration file for Nucleus: "
                  << name << std::endl;
        std::cout << "Generate configuration with Wood-Saxon distribution"
                  << std::endl;
        return;
    }

    std::ifstream input(filename.str().c_str());
    if (!input.good()) {
        std::cout << "Configuration file not found!" << std::endl;
        std::cout << "Please check file: " << filename.str()
                  << std::endl;
        exit(1);
    }

    double dummy;
    for (int iconf = 0; iconf < n_configuration; iconf++) {
        std::vector< std::array<double, 3> > conf_i;

        if (A_ == 12) {
            input >> dummy >> dummy;
        }
        for (int ia = 0; ia < A_; ia++) {
            double x_local, y_local, z_local;
            int isospin;
            if (A_ == 208) {
                input >> x_local >> y_local >> z_local >> isospin;
            } else if (A_ == 197) {
                input >> x_local >> y_local >> z_local >> isospin >> dummy;
            } else {
                input >> x_local >> y_local >> z_local;
            }
            std::array<double, 3> nucleon_pos = {x_local, y_local, z_local};
            conf_i.push_back(nucleon_pos);
        }
        if (A_ == 3) {
            input >> dummy >> dummy >> dummy >> dummy;
        }
        heavyIon_pos_.push_back(conf_i);
    }
    input.close();
    cout << heavyIon_pos_.size() << " configrations." << endl;
}


int Nucleus::sample_nucleon_configuration() {
    // This function samples the spatial configuration for triton
    if (!nucleon_configuration_loaded_) {
        readin_nucleon_positions();
        nucleon_configuration_loaded_ = true;
    }
    if (heavyIon_pos_.size() == 0)
        return(2);
    const int num_configuration = heavyIon_pos_.size();
    const int rand_num = static_cast<int>(
                            ran_gen_ptr->rand_uniform()*num_configuration);
    auto conf_i = heavyIon_pos_[rand_num];
    for (int iA = 0; iA < A_; iA++) {
        SpatialVec x_i = {0.0, conf_i[iA][0], conf_i[iA][1], conf_i[iA][2]};
        MomentumVec p_i = {0.0};

        std::shared_ptr<Nucleon> nucleon_i_ptr(new Nucleon(x_i, p_i));
        nucleon_list_.push_back(std::move(nucleon_i_ptr));
    }
    return(0);
}


void Nucleus::sample_r_from_woods_saxon(
            std::vector<std::pair<real, real>> &r_array) const {
    real a_WS_p = WS_param_vec[3];
    real a_WS_n = WS_param_vec[3] + WS_param_vec[8];
    real R_WS_p = WS_param_vec[2];
    real R_WS_n = WS_param_vec[2] + WS_param_vec[9];
    for (int i = 0; i < Z_; i++) {
        real rmaxCut = R_WS_p + 10.*a_WS_p;
        real r = 0.;
        do {
            r = rmaxCut*pow(ran_gen_ptr->rand_uniform(), 1.0/3.0);
        } while (ran_gen_ptr->rand_uniform()
                        > fermi_distribution(r, R_WS_p, a_WS_p));
        r_array.push_back(std::make_pair(r, 1.));
    }

    for (int i = Z_; i < A_; i++) {
        real rmaxCut = R_WS_n + 10.*a_WS_n;
        real r = 0.;
        do {
            r = rmaxCut*pow(ran_gen_ptr->rand_uniform(), 1.0/3.0);
        } while (ran_gen_ptr->rand_uniform()
                        > fermi_distribution(r, R_WS_n, a_WS_n));
        r_array.push_back(std::make_pair(r, 0));
    }
}

real Nucleus::sample_r_from_deformed_woods_saxon() const {
    real a_WS = WS_param_vec[3];
    real R_WS = WS_param_vec[2];
    real beta2 = WS_param_vec[4];
    real beta3 = WS_param_vec[5];
    real beta4 = WS_param_vec[6];
    real gamma = WS_param_vec[7];
    real rmaxCut = R_WS + 10.*a_WS;
    real R_WS_theta = R_WS;
    real r = 0.;
    do {
        r = rmaxCut*pow(ran_gen_ptr->rand_uniform(), 1.0/3.0);
        real costheta = 1.0 - 2.0*ran_gen_ptr->rand_uniform();
        real phi  = 2.*M_PI*ran_gen_ptr->rand_uniform();
        real y20  = spherical_harmonics(2, costheta);
        real y30  = spherical_harmonics(3, costheta);
        real y40  = spherical_harmonics(4, costheta);
        real y2_2 = spherical_harmonics_Y22(22, costheta, phi);
        R_WS_theta = R_WS*(1.0
                           + beta2*(cos(gamma)*y20+sin(gamma)*y2_2)
                           + beta3*y30 + beta4*y40);
    } while (ran_gen_ptr->rand_uniform()
             > fermi_distribution(r, R_WS_theta, a_WS));
    return(r);
}

void Nucleus::sample_r_and_costheta_from_deformed_woods_saxon(
        std::vector<std::array<real, 4>> &nucleonPos_array) const {
    real a_WS_p = WS_param_vec[3];
    real a_WS_n = WS_param_vec[3] + WS_param_vec[8];
    real R_WS_p = WS_param_vec[2];
    real R_WS_n = WS_param_vec[2] + WS_param_vec[9];
    real beta2 = WS_param_vec[4];
    real beta3 = WS_param_vec[5];
    real beta4 = WS_param_vec[6];
    real gamma = WS_param_vec[7];
    for (int i = 0; i < Z_; i++) {
        real rmaxCut = R_WS_p + 10.*a_WS_p;
        real R_WS_theta = R_WS_p;
        real r, costheta, phi;
        do {
            r = rmaxCut*pow(ran_gen_ptr->rand_uniform(), 1.0/3.0);
            costheta  = 1.0 - 2.0*ran_gen_ptr->rand_uniform();
            phi       = 2.*M_PI*ran_gen_ptr->rand_uniform();
            real y20  = spherical_harmonics(2, costheta);
            real y30  = spherical_harmonics(3, costheta);
            real y40  = spherical_harmonics(4, costheta);
            real y2_2 = spherical_harmonics_Y22(22, costheta, phi);
            R_WS_theta = R_WS_p*(1.0
                                 + beta2*(cos(gamma)*y20+sin(gamma)*y2_2)
                                 + beta3*y30 + beta4*y40);
        } while (ran_gen_ptr->rand_uniform()
                 > fermi_distribution(r, R_WS_theta, a_WS_p));
        std::array<real, 4> sample_i = {r, costheta, phi, 1.};
        nucleonPos_array.push_back(sample_i);
    }

    for (int i = Z_; i < A_; i++) {
        real rmaxCut = R_WS_n + 10.*a_WS_n;
        real R_WS_theta = R_WS_n;
        real r, costheta, phi;
        do {
            r = rmaxCut*pow(ran_gen_ptr->rand_uniform(), 1.0/3.0);
            costheta  = 1.0 - 2.0*ran_gen_ptr->rand_uniform();
            phi       = 2.*M_PI*ran_gen_ptr->rand_uniform();
            real y20  = spherical_harmonics(2, costheta);
            real y30  = spherical_harmonics(3, costheta);
            real y40  = spherical_harmonics(4, costheta);
            real y2_2 = spherical_harmonics_Y22(22, costheta, phi);
            R_WS_theta = R_WS_n*(1.0
                                 + beta2*(cos(gamma)*y20+sin(gamma)*y2_2)
                                 + beta3*y30 + beta4*y40);
        } while (ran_gen_ptr->rand_uniform()
                 > fermi_distribution(r, R_WS_theta, a_WS_n));
        std::array<real, 4> sample_i = {r, costheta, phi, 0.};
        nucleonPos_array.push_back(sample_i);
    }
}


void Nucleus::generate_nucleus_configuration_with_woods_saxon() {
    std::vector<std::pair<real, real>> r_array;
    sample_r_from_woods_saxon(r_array);
    std::sort(r_array.begin(), r_array.end());

    std::vector<real> x_array(A_, 0.), y_array(A_, 0.), z_array(A_, 0.);
    const real d_min_sq = d_min_*d_min_;
    for (unsigned int i = 0; i < r_array.size(); i++) {
        real r_i = std::get<0>(r_array[i]);
        int reject_flag = 0;
        int iter = 0;
        real x_i, y_i, z_i;
        do {
            iter++;
            reject_flag = 0;
            real phi    = 2.*M_PI*ran_gen_ptr->rand_uniform();
            real theta  = acos(1. - 2.*ran_gen_ptr->rand_uniform());
            x_i = r_i*sin(theta)*cos(phi);
            y_i = r_i*sin(theta)*sin(phi);
            z_i = r_i*cos(theta);
            for (int j = i - 1; j >= 0; j--) {
                real r_j = std::get<0>(r_array[j]);
                if ((r_i - r_j)*(r_i - r_j) > d_min_sq) break;
                real dsq = (  (x_i - x_array[j])*(x_i - x_array[j])
                            + (y_i - y_array[j])*(y_i - y_array[j])
                            + (z_i - z_array[j])*(z_i - z_array[j]));
                if (dsq < d_min_sq) {
                    reject_flag = 1;
                    break;
                }
            }
        } while (reject_flag == 1 && iter < 100);
        //if (iter == 100) {
        //    cout << "[Warning] can not find configuration : "
        //         << "r[i] = " << r_i << ", r[i-1] = " << r_array[i-1] << endl;
        //}
        x_array[i] = x_i;
        y_array[i] = y_i;
        z_array[i] = z_i;
    }
    for (unsigned int i = 0; i < r_array.size(); i++) {
        SpatialVec  x_in = {0.0, x_array[i], y_array[i], z_array[i]};
        MomentumVec p_in = {0.0};
        std::shared_ptr<Nucleon> nucleon_ptr(new Nucleon(x_in, p_in));
        nucleon_ptr->set_electric_charge(std::get<1>(r_array[i]));
        nucleon_list_.push_back(std::move(nucleon_ptr));
    }
    set_nucleons_momentum_with_collision_energy(0.0);
}


void Nucleus::generate_nucleus_configuration_with_deformed_woods_saxon() {
    std::vector<std::array<real, 4>> nucleonPos_array;
    sample_r_and_costheta_from_deformed_woods_saxon(nucleonPos_array);
    std::sort(nucleonPos_array.begin(), nucleonPos_array.end());

    std::vector<real> x_array(A_, 0.), y_array(A_, 0.), z_array(A_, 0.);
    const real d_min_sq = d_min_*d_min_;
    for (int i = 0; i < A_; i++) {
        const real r_i     = nucleonPos_array[i][0];
        const real theta_i = acos(nucleonPos_array[i][1]);
        real phi_i = nucleonPos_array[i][2];
        int reject_flag = 0;
        int iter = 0;
        real x_i, y_i, z_i;
        do {
            iter++;
            reject_flag = 0;
            x_i = r_i*sin(theta_i)*cos(phi_i);
            y_i = r_i*sin(theta_i)*sin(phi_i);
            z_i = r_i*cos(theta_i);
            for (int j = i - 1; j >= 0; j--) {
                const real r_j = nucleonPos_array[j][0];
                if ((r_i - r_j)*(r_i - r_j) > d_min_sq) break;
                real dsq = (  (x_i - x_array[j])*(x_i - x_array[j])
                            + (y_i - y_array[j])*(y_i - y_array[j])
                            + (z_i - z_array[j])*(z_i - z_array[j]));
                if (dsq < d_min_sq) {
                    reject_flag = 1;
                    break;
                }
            }
            phi_i = 2.*M_PI*ran_gen_ptr->rand_uniform();
        } while (reject_flag == 1 && iter < 100);
        //if (iter == 100) {
        //    cout << "[Warning] can not find configuration : "
        //         << "r[i] = " << r_i << ", r[i-1] = " << r_array[i-1] << endl;
        //}
        x_array[i] = x_i;
        y_array[i] = y_i;
        z_array[i] = z_i;
    }
    for (int i = 0; i < A_; i++) {
        SpatialVec  x_in = {0.0, x_array[i], y_array[i], z_array[i]};
        MomentumVec p_in = {0.0};
        std::shared_ptr<Nucleon> nucleon_ptr(new Nucleon(x_in, p_in));
        nucleon_ptr->set_electric_charge(nucleonPos_array[i][3]);
        nucleon_list_.push_back(std::move(nucleon_ptr));
    }
    set_nucleons_momentum_with_collision_energy(0.0);
}


real Nucleus::fermi_distribution(real r, real R_WS, real a_WS) const {
    real f = 1./(1. + exp((r - R_WS)/(a_WS + 1e-16)));
    return (f);
}


real Nucleus::getAvgWoodsSaxonDensity(real r) const {
    const real rho_0 = WS_param_vec[0];  // 1/fm^3
    real density = rho_0*fermi_distribution(r, WS_param_vec[2],
                                            WS_param_vec[3]);
    return(density);  // 1/fm^3
}


real Nucleus::spherical_harmonics(int l, real ct) const {
    // Currently assuming m=0 and available for Y_{20} and Y_{40}
    // "ct" is cos(theta)
    assert(l == 2 || l == 3 || l == 4);

    real ylm = 0.0;
    if (l == 2) {
        ylm = 3.0*ct*ct-1.0;
        ylm *= 0.31539156525252005;  // pow(5.0/16.0/M_PI,0.5);
    } else if (l == 3) {
        ylm  = 5.0*ct*ct*ct;
        ylm -= 3.0*ct;
        ylm *= 0.3731763325901154;  // pow(7.0/16.0/M_PI,0.5);
    } else if (l == 4) {
        ylm  = 35.0*ct*ct*ct*ct;
        ylm -= 30.0*ct*ct;
        ylm += 3.0;
        ylm *= 0.10578554691520431;  // 3.0/16.0/pow(M_PI,0.5);
    }
    return(ylm);
}

real Nucleus::spherical_harmonics_Y22(int l, real ct, real phi) const {
    // Y2,2
    // "ct" is cos(theta)
    assert(l == 22);

    real ylm = 0.0;
    ylm  = 1.0-ct*ct;
    ylm *= cos(2.*phi);
    ylm *= 0.5462742152960397;  // sqrt(2*15)/4/pow(2*M_PI,0.5);
    return(ylm);
}


void Nucleus::accelerate_nucleus(real ecm, int direction) {
    assert(ecm > 2.*PhysConsts::MProton);
    real beam_rapidity = direction*acosh(ecm/(2.*PhysConsts::MProton));
    set_nucleons_momentum_with_collision_energy(beam_rapidity);
    lorentz_contraction(cosh(beam_rapidity));
}


void Nucleus::lorentz_contraction(real gamma) {
    for (auto &it: nucleon_list_) {
        auto xvec = it->get_x();
        xvec[3] /= gamma;
        it->set_x(xvec);
    }

    if (sample_valence_quarks) {
        for (auto &it: nucleon_list_)
            it->lorentz_contraction(gamma);
    }
}


void Nucleus::set_nucleons_momentum_with_collision_energy(real beam_rapidity) {
    MomentumVec p = {PhysConsts::MProton*cosh(beam_rapidity),
                     0.0,
                     0.0,
                     PhysConsts::MProton*sinh(beam_rapidity)};
    for (auto &it: nucleon_list_) {
        it->set_p(p);
        it->set_remnant_p(p);
    }
}


real Nucleus::get_z_min() const {
    real z_min = 100;
    for (auto const &it: nucleon_list_) {
        auto xvec = it->get_x();
        if (xvec[3] < z_min) z_min = xvec[3];
    }
    return(z_min);
}


real Nucleus::get_z_max() const {
    real z_max = -100;
    for (auto const &it: nucleon_list_) {
        auto xvec = it->get_x();
        if (xvec[3] > z_max) z_max = xvec[3];
    }
    return(z_max);
}


void Nucleus::output_nucleon_positions(std::string filename) const {
    std::ofstream of(filename, std::ofstream::out);
    of << "# Nucleus name: " << name << std::endl;
    of << "# x (fm)  y (fm)  z (fm)  rapidity  electric_charge" << std::endl;
    for (auto const& it: nucleon_list_) {
        auto x_vec = it->get_x();
        auto p_vec = it->get_p();
        real rapidity = 0.5*log((p_vec[0] + p_vec[3])/(p_vec[0] - p_vec[3]));
        of << std::scientific << std::setw(10) << std::setprecision(6)
           << x_vec[1] << "  " << x_vec[2] << "  " << x_vec[3] << "  "
           << rapidity << "  " << it->get_electric_charge()
           << std::endl;
    }
    of.close();
}


void Nucleus::sample_quark_momentum_fraction(std::vector<real> &xQuark,
                                             std::vector<real> &eQuark,
                                             const int number_of_quarks,
                                             const real electric_charge,
                                             const real ecm) const {
    if (!sample_valence_quarks) {
        for (int i = 0; i < number_of_quarks; i++) {
            xQuark.push_back(1./number_of_quarks);
            eQuark.push_back(electric_charge/number_of_quarks);
        }
        return;
    }

    const real mq = PhysConsts::MQuarkValence;
    const real mp = PhysConsts::MProton;
    const real ybeam = acosh(ecm/(2.*mp));
    real total_energy = 0.;
    real E_proton = mp*cosh(ybeam);
    do {
        auto sample_idx = static_cast<int>(
            ran_gen_ptr->rand_uniform()*number_of_valence_quark_samples_);
        xQuark.clear();
        if (std::abs(electric_charge - 1) < 1e-8) {
            for (int i = 0; i < number_of_quarks; i++) {
                xQuark.push_back(proton_valence_quark_x_[sample_idx][i]);
                if (i == 0) {
                    eQuark.push_back(-1./3.);
                } else if (i < 3) {
                    eQuark.push_back(2./3.);
                } else {
                    eQuark.push_back(0.);
                }
            }
        } else  {
            for (int i = 0; i < number_of_quarks; i++) {
                xQuark.push_back(neutron_valence_quark_x_[sample_idx][i]);
                if (i < 2) {
                    eQuark.push_back(-1./3.);
                } else if (i == 2) {
                    eQuark.push_back(2./3.);
                } else {
                    eQuark.push_back(0.);
                }
            }
        }
        total_energy = 0.;
        for (int i = 0; i < number_of_quarks; i++) {
            real rap_local = asinh(xQuark[i]*mp/mq*sinh(ybeam));
            real E_local = mq*cosh(rap_local);
            total_energy += E_local;
        }
    } while (total_energy > E_proton);
}


SpatialVec Nucleus::sample_valence_quark_position() const {
    // sample Gaussian distribution for the valence quark position
    // determine x,y,z coordinates of the quark (relative to the nucleon)

    real BG = sqrt(BG_)*PhysConsts::HBARC;
    real x = ran_gen_ptr->rand_normal(0., BG);
    real y = ran_gen_ptr->rand_normal(0., BG);
    real z = ran_gen_ptr->rand_normal(0., BG);

    SpatialVec xq = {0.0, x, y, z};
    return(xq);
}


real Nucleus::ExponentialDistribution(const real a, const real r) const {
    // a = \sqrt{12}/R_p = 3.87, with R_p = 0.895
    // from http://journals.aps.org/rmp/pdf/10.1103/RevModPhys.77.1
    return(r*r*exp(- a*r));
}

}

