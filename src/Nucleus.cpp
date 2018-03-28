// Copyright @ Chun Shen 2018

#include "Nucleus.h"
#include "PhysConsts.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>

using std::cout;
using std::endl;

namespace MCGlb {

Nucleus::Nucleus(std::string nucleus_name, int seed_in, real d_min_in,
                 bool deformed_in) {
    d_min = d_min_in;
    deformed = deformed_in;
    set_random_seed(seed_in);
    set_nucleus_parameters(nucleus_name);
}

Nucleus::~Nucleus() {
    nucleon_list.clear();
}

void Nucleus::set_random_seed(int seed) {
    ran_gen_ptr = (
        std::unique_ptr<RandomUtil::Random>(new RandomUtil::Random(seed)));
}


void Nucleus::set_woods_saxon_parameters(int A_in, int Z_in,
                                         real rho, real w, real R, real a,
                                         real beta2, real beta4,
                                         int density_function_type_in) {
    A                     = A_in;
    Z                     = Z_in;
    WS_param_vec[0]       = rho;
    WS_param_vec[1]       = w;
    WS_param_vec[2]       = R;
    WS_param_vec[3]       = a;
    WS_param_vec[4]       = beta2;
    WS_param_vec[5]       = beta4;
    density_function_type = density_function_type_in;
}


void Nucleus::set_nucleus_parameters(std::string nucleus_name) {
    name = nucleus_name;
    if (nucleus_name.compare("p") == 0) {
        set_woods_saxon_parameters(1, 1, 0.17, 0.0, 1.0, 1.0, 0.0, 0.0, 3);
    } else if (nucleus_name.compare("d") == 0) {
        set_woods_saxon_parameters(2, 1, 0.17, 1.18, 1.0, 0.228, 0.0, 0.0, 8);
     } else if (nucleus_name.compare("He3") == 0) {
        set_woods_saxon_parameters(3, 2, 0.17, 0.0, 0.0, 0.0, 0.0, 0.0, 1);
    } else if (nucleus_name.compare("C") == 0) {
        set_woods_saxon_parameters(
                            12, 6, 0.17, 1.403, 2.44, 1.635, 0.0, 0.0, 1);
    } else if (nucleus_name.compare("O") == 0) {
        set_woods_saxon_parameters(
                            16, 8, 0.17, 2.608, -0.051, 0.513, 0.0, 0.0, 3);
    } else if (nucleus_name.compare("Al") == 0) {
        set_woods_saxon_parameters(
                            27, 13, 0.17, 3.07, 0.0, 0.519, 0.0, 0.0, 3);
    } else if (nucleus_name.compare("Cu") == 0) {
        set_woods_saxon_parameters(
                            63, 29, 0.17, 4.163, 0.0, 0.606, 0.162, 0.006, 3);
    } else if (nucleus_name.compare("Au") == 0) {
        set_woods_saxon_parameters(
                            197, 79, 0.17, 0.0, 6.38, 0.505, -0.13, -0.03, 3);
    } else if (nucleus_name.compare("Pb") == 0) {
        set_woods_saxon_parameters(
                            208, 82, 0.17, 0.0, 6.62, 0.546, 0.0, 0.0, 3);
    } else if (nucleus_name.compare("U") == 0) {
        set_woods_saxon_parameters(
                            238, 92, 0.17, 0.0, 6.874, 0.556, 0.28, 0.093, 3);
    } else {
        cout << "[Error] Unknown_nucleus: " << nucleus_name << endl;
        cout << "Exiting... " << endl;
        exit(1);
    }
}


void Nucleus::generate_nucleus_3d_configuration() {
    if (nucleon_list.size() > 0) {
        nucleon_list.clear();
    }
    if (A == 1) {  // p
        SpatialVec  x = {0.0};
        MomentumVec p = {0.0};
        Nucleon nucleon(x, p);
        nucleon_list.push_back(nucleon);
    } else if (A == 2) {  // deuteron
        generate_deuteron_configuration();
    } else {  // other nucleus
        if (!deformed) {
            generate_nucleus_configuration_with_woods_saxon();
        } else {
            generate_nucleus_configuration_with_deformed_woods_saxon();
        }
    }
    recenter_nucleus();
}


void Nucleus::recenter_nucleus() {
    // compute the center of mass position and shift it to (0, 0, 0)
    real meanx = 0., meany = 0., meanz = 0.;
    for (auto const &nucleon_i : nucleon_list) {
        auto x_vec = nucleon_i.get_x();
        meanx += x_vec[1];
        meany += x_vec[2];
        meanz += x_vec[3];
    }
      
    meanx /= static_cast<real>(A);
    meany /= static_cast<real>(A);
    meanz /= static_cast<real>(A);

    SpatialVec x_shift = {0, -meanx, -meany, -meanz};
    shift_nucleus(x_shift);
}


void Nucleus::shift_nucleus(SpatialVec x_shift) {
    for (auto &nucleon_i : nucleon_list) {
        auto x_vec = nucleon_i.get_x();
        for (int i = 0; i < 4; i++)
            x_vec[i] += x_shift[i];
        nucleon_i.set_x(x_vec);
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

    Nucleon nucleon_1(x_1, p_1);
    Nucleon nucleon_2(x_2, p_2);
    nucleon_list.push_back(nucleon_1);
    nucleon_list.push_back(nucleon_2);
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


real Nucleus::sample_r_from_woods_saxon() const {
    real a_WS = WS_param_vec[3];
    real R_WS = WS_param_vec[2];
    real rmaxCut = R_WS + 10.*a_WS;
    real r = 0.;
    do {
        r = rmaxCut*pow(ran_gen_ptr->rand_uniform(), 1.0/3.0);
    } while (ran_gen_ptr->rand_uniform() > fermi_distribution(r, R_WS, a_WS));
    return(r);
}

void Nucleus::sample_r_and_costheta_from_deformed_woods_saxon(
                                        real &r, real &costheta) const {
    real a_WS = WS_param_vec[3];
    real R_WS = WS_param_vec[2];
    real beta2 = WS_param_vec[4];
    real beta4 = WS_param_vec[5];
    real rmaxCut = R_WS + 10.*a_WS;
    real R_WS_theta = R_WS;
    do {
        r = rmaxCut*pow(ran_gen_ptr->rand_uniform(), 1.0/3.0);
        costheta = 1.0 - 2.0*ran_gen_ptr->rand_uniform();
        real y20 = spherical_harmonics(2, costheta);
        real y40 = spherical_harmonics(4, costheta);
        R_WS_theta = R_WS*(1.0 + beta2*y20 + beta4*y40);
    } while (ran_gen_ptr->rand_uniform()
             > fermi_distribution(r, R_WS_theta, a_WS));
}


void Nucleus::generate_nucleus_configuration_with_woods_saxon() {
    std::vector<real> r_array(A, 0.);
    for (int i = 0; i < A; i++)
        r_array[i] = sample_r_from_woods_saxon();
    std::sort(r_array.begin(), r_array.end());

    std::vector<real> x_array(A, 0.), y_array(A, 0.), z_array(A, 0.);
    const real d_min_sq = d_min*d_min;
    for (unsigned int i = 0; i < r_array.size(); i++) {
        real r_i = r_array[i];
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
                if ((r_i - r_array[j])*(r_i - r_array[j]) > d_min_sq) break;
                real dsq = (  (x_i - x_array[j])*(x_i - x_array[j])
                            + (y_i - y_array[j])*(y_i - y_array[j])
                            + (z_i - z_array[j])*(z_i - z_array[j]));
                if (dsq < d_min_sq) {
                    reject_flag = 1;
                    break;
                }
            }
        } while (reject_flag == 1 && iter < 100);
        if (iter == 100) {
            cout << "[Warning] can not find configuration : "
                 << "r[i] = " << r_i << ", r[i-1] = " << r_array[i-1] << endl;
        }
        x_array[i] = x_i;
        y_array[i] = y_i;
        z_array[i] = z_i;
    }
    for (unsigned int i = 0; i < r_array.size(); i++) {
        SpatialVec  x_in = {0.0, x_array[i], y_array[i], z_array[i]};
        MomentumVec p_in = {0.0};
        Nucleon nucleon(x_in, p_in);
        nucleon_list.push_back(nucleon);
    }
}


void Nucleus::generate_nucleus_configuration_with_deformed_woods_saxon() {
    std::vector<real> r_array(A, 0.);
    std::vector<real> costheta_array(A, 0.);
    for (int i = 0; i < A; i++) {
        sample_r_and_costheta_from_deformed_woods_saxon(r_array[i],
                                                        costheta_array[i]);
    }
    std::sort(r_array.begin(), r_array.end());

    std::vector<real> x_array(A, 0.), y_array(A, 0.), z_array(A, 0.);
    const real d_min_sq = d_min*d_min;
    for (unsigned int i = 0; i < r_array.size(); i++) {
        const real r_i     = r_array[i];
        const real theta_i = acos(costheta_array[i]);
        int reject_flag = 0;
        int iter = 0;
        real x_i, y_i, z_i;
        do {
            iter++;
            reject_flag = 0;
            real phi    = 2.*M_PI*ran_gen_ptr->rand_uniform();
            x_i = r_i*sin(theta_i)*cos(phi);
            y_i = r_i*sin(theta_i)*sin(phi);
            z_i = r_i*cos(theta_i);
            for (int j = i - 1; j >= 0; j--) {
                if ((r_i - r_array[j])*(r_i - r_array[j]) > d_min_sq) break;
                real dsq = (  (x_i - x_array[j])*(x_i - x_array[j])
                            + (y_i - y_array[j])*(y_i - y_array[j])
                            + (z_i - z_array[j])*(z_i - z_array[j]));
                if (dsq < d_min_sq) {
                    reject_flag = 1;
                    break;
                }
            }
        } while (reject_flag == 1 && iter < 100);
        if (iter == 100) {
            cout << "[Warning] can not find configuration : "
                 << "r[i] = " << r_i << ", r[i-1] = " << r_array[i-1] << endl;
        }
        x_array[i] = x_i;
        y_array[i] = y_i;
        z_array[i] = z_i;
    }
    for (unsigned int i = 0; i < r_array.size(); i++) {
        SpatialVec  x_in = {0.0, x_array[i], y_array[i], z_array[i]};
        MomentumVec p_in = {0.0};
        Nucleon nucleon(x_in, p_in);
        nucleon_list.push_back(nucleon);
    }
}


real Nucleus::fermi_distribution(real r, real R_WS, real a_WS) const {
    real f = 1./(1. + exp((r - R_WS)/a_WS));
    return (f);
}

real Nucleus::spherical_harmonics(int l, real ct) const {
    // Currently assuming m=0 and available for Y_{20} and Y_{40}
    // "ct" is cos(theta)
    assert(l == 2 || l == 4);

    real ylm = 0.0;
    if (l == 2) {
        ylm = 3.0*ct*ct-1.0;
        ylm *= 0.31539156525252005;  // pow(5.0/16.0/M_PI,0.5);
    } else if (l == 4) {
        ylm  = 35.0*ct*ct*ct*ct;
        ylm -= 30.0*ct*ct;
        ylm += 3.0;
        ylm *= 0.10578554691520431;  // 3.0/16.0/pow(M_PI,0.5);
    }
    return(ylm);
}


void Nucleus::accelerate_nucleus(real ecm, int direction) {
    assert(ecm > 2.*PhysConsts::MProton);
    real beam_rapidity = direction*acosh(ecm/(2.*PhysConsts::MProton));
    set_nucleons_momentum_with_collision_energy(beam_rapidity);
    lorentz_contraction(cosh(beam_rapidity));
}

void Nucleus::lorentz_contraction(real gamma) {
    for (auto &it: nucleon_list) {
        auto xvec = it.get_x();
        xvec[3] /= gamma;
        it.set_x(xvec);
    }
}

void Nucleus::set_nucleons_momentum_with_collision_energy(real beam_rapidity) {
    MomentumVec p = {PhysConsts::MProton*cosh(beam_rapidity),
                     0.0,
                     0.0,
                     PhysConsts::MProton*sinh(beam_rapidity)};
    for (auto &it: nucleon_list)
        it.set_p(p);
}

real Nucleus::get_z_min() const {
    real z_min = 100;
    for (auto const &it: nucleon_list) {
        auto xvec = it.get_x();
        if (xvec[3] < z_min) z_min = xvec[3];
    }
    return(z_min);
}

real Nucleus::get_z_max() const {
    real z_max = -100;
    for (auto const &it: nucleon_list) {
        auto xvec = it.get_x();
        if (xvec[3] > z_max) z_max = xvec[3];
    }
    return(z_max);
}

void Nucleus::output_nucleon_positions(std::string filename) const {
    std::ofstream of(filename, std::ofstream::out);
    of << "# Nucleus name: " << name << std::endl;
    of << "# x (fm)  y (fm)  z (fm)  rapidity" << std::endl;
    for (auto const& it: nucleon_list) {
        auto x_vec = it.get_x();
        auto p_vec = it.get_p();
        real rapidity = 0.5*log((p_vec[0] + p_vec[3])/(p_vec[0] - p_vec[3]));
        of << std::scientific << std::setw(10) << std::setprecision(6)
           << x_vec[1] << "  " << x_vec[2] << "  " << x_vec[3] << "  "
           << rapidity << std::endl;
    }
}
    
int Nucleus::get_number_of_wounded_nucleons() const {
    int Npart = 0;
    for (auto const& it: nucleon_list)
        if (it.is_wounded()) Npart++;
    return(Npart);
}

}

