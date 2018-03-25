// Copyright @ Chun Shen 2018

#include "Nucleus.h"
#include "PhysConsts.h"
#include <iostream>
#include <algorithm>
#include <cmath>

using std::cout;
using std::endl;

namespace MCGlb {

Nucleus::Nucleus(std::string nucleus_name, int seed_in, real d_min_in) {
    d_min = d_min_in;
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
                                         int density_function_type_in) {
    A                     = A_in;
    Z                     = Z_in;
    WS_param_vec[0]       = rho;
    WS_param_vec[1]       = w;
    WS_param_vec[2]       = R;
    WS_param_vec[3]       = a;
    density_function_type = density_function_type_in;
}


void Nucleus::set_nucleus_parameters(std::string nucleus_name) {
    name = nucleus_name;
    if (nucleus_name.compare("p") == 0) {
        set_woods_saxon_parameters(1, 1, 0.17, 0.0, 1.0, 1.0, 3);
    } else if (nucleus_name.compare("d") == 0) {
        set_woods_saxon_parameters(2, 1, 0.17, 1.18, 1.0, 0.228, 8);
     } else if (nucleus_name.compare("He3") == 0) {
        set_woods_saxon_parameters(3, 2, 0.17, 0.0, 0.0, 0.0, 1);
    } else if (nucleus_name.compare("C") == 0) {
        set_woods_saxon_parameters(12, 6, 0.17, 1.403, 2.44, 1.635, 1);
    } else if (nucleus_name.compare("O") == 0) {
        set_woods_saxon_parameters(16, 8, 0.17, 2.608, -0.051, 0.513, 3);
    } else if (nucleus_name.compare("Al") == 0) {
        set_woods_saxon_parameters(27, 13, 0.17, 3.07, 0.0, 0.519, 3);
    } else if (nucleus_name.compare("Cu") == 0) {
        set_woods_saxon_parameters(63, 29, 0.17, 4.163, 0.0, 0.606, 3);
    } else if (nucleus_name.compare("Au") == 0) {
        set_woods_saxon_parameters(197, 79, 0.17, 0.0, 6.38, 0.505, 3);
    } else if (nucleus_name.compare("Pb") == 0) {
        set_woods_saxon_parameters(208, 82, 0.17, 0.0, 6.62, 0.546, 3);
    } else if (nucleus_name.compare("U") == 0) {
        set_woods_saxon_parameters(238, 92, 0.17, 0.0, 6.874, 0.556, 3);
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
        generate_nucleus_configuration_with_woods_saxon();
    }
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

    shift_nucleus(0, -meanx, -meany, -meanz);
}


void Nucleus::shift_nucleus(real t_shift, real x_shift, real y_shift,
                            real z_shift) {
    for (auto &nucleon_i : nucleon_list) {
        auto x_vec = nucleon_i.get_x();
        x_vec[0] += t_shift;
        x_vec[1] += x_shift;
        x_vec[2] += y_shift;
        x_vec[3] += z_shift;
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


real Nucleus::get_inverse_CDF_hulthen_function(real y) {
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


real Nucleus::hulthen_function_CDF(real r) {
    real alpha = 0.228;
    real beta  = 1.18;
    real res   = 0.0;
    if (r < 0.) {
        res = 0.0;
    } else {
        real c = alpha*beta*(alpha + beta)/((alpha - beta)*(alpha - beta));
        res = (2.*c*(2.*exp(-r*(alpha + beta))/(alpha + beta)
                     - 1./2.*exp(-2.*alpha*r)/alpha
                     - 1./2.*exp(-2.*beta*r)/beta
                     + 1./(2.*alpha) + 1./(2.*beta) - 2./(alpha + beta)));
    }
    return(res);
}


real Nucleus::sample_r_from_woods_saxon() {
    real a_WS = WS_param_vec[3];
    real R_WS = WS_param_vec[2];
    real rmaxCut = R_WS + 10.*a_WS;
    real r = 0.;
    do {
        r = rmaxCut*pow(ran_gen_ptr->rand_uniform(), 1.0/3.0);
    } while (ran_gen_ptr->rand_uniform() > fermi_distribution(r, R_WS, a_WS));
    return(r);
}


void Nucleus::generate_nucleus_configuration_with_woods_saxon() {
    std::vector<real> r_array;
    for (int i = 0; i < A; i++) {
        r_array.push_back(sample_r_from_woods_saxon());
    }
    std::sort(r_array.begin(), r_array.end());

    std::vector<real> x_array, y_array, z_array;
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
        x_array.push_back(x_i);
        y_array.push_back(y_i);
        z_array.push_back(z_i);
    }
    for (unsigned int i = 0; i < r_array.size(); i++) {
        SpatialVec  x_in = {0.0, x_array[i], y_array[i], z_array[i]};
        MomentumVec p_in = {0.0};
        Nucleon nucleon(x_in, p_in);
        nucleon_list.push_back(nucleon);
    }
}


real Nucleus::fermi_distribution(real r, real R_WS, real a_WS) {
    real f = 1./(1. + exp((r - R_WS)/a_WS));
    return (f);
}

}

