// Copyright @ Chun Shen 2018

#include "Nucleon.h"
#include "PhysConsts.h"
#include<cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <utility>

#include "eps09.h"
#include "LHAPDF/LHAPDF.h"
int MCGlb::Nucleon::random_value_ = 0; 

namespace MCGlb {

Nucleon::Nucleon(SpatialVec x_in, MomentumVec p_in,
                 std::shared_ptr<RandomUtil::Random> ran_gen_ptr) {
    ran_gen_ptr_ = ran_gen_ptr;
    set_particle_variables(x_in, p_in);
    number_of_valence_quark_resamples_ = re_readin_valence_quark_samples();
}


Nucleon::Nucleon(SpatialVec x_in, MomentumVec p_in, real mass_in,
                 std::shared_ptr<RandomUtil::Random> ran_gen_ptr) {
    ran_gen_ptr_ = ran_gen_ptr;
    set_particle_variables(x_in, p_in, mass_in);
}

Nucleon::~Nucleon() {
    quark_list.clear();
    proton_resample_quark_x_.clear();
    neutron_resample_quark_x_.clear();
}

bool Nucleon::is_connected_with(std::shared_ptr<Nucleon> targ) {
    bool connected = false;
    for (auto &it: connected_with) {
        if (*(it.lock()) == *targ) {
            connected = true;
            break;
        }
    }
    return(connected);
}


int Nucleon::get_number_of_connections(std::shared_ptr<Nucleon> targ) const {
    int n_connections = 0;
    for (unsigned int idx = 0; idx < connected_with.size(); idx++) {
        if (*(connected_with[idx].lock()) == *targ) {
            n_connections = connected_times_[idx];
            break;
        }
    }
    return(n_connections);
}

void Nucleon::accelerate_quarks(real ecm, int direction) {
    const real mq = PhysConsts::MQuarkValence;
    const real mp = PhysConsts::MProton;
    const real ybeam = acosh(ecm/(2.*mp));
    for (auto &it: quark_list) {
        //real rap_local = direction*asinh(it->get_pdf_x()
        //                                 *sqrt(ecm*ecm/(4.*mq*mq) - 1.));
        real rap_local = direction*asinh(it->get_pdf_x()*mp/mq*sinh(ybeam));
        it->set_rapidity(rap_local);
        MomentumVec p_in = {mq*cosh(rap_local), 0.0, 0.0, mq*sinh(rap_local)};
        it->set_p(p_in);
    }
}

int Nucleon::re_readin_valence_quark_samples() {
    std::stringstream of_p_name_resample, of_n_name_resample;
    of_p_name_resample << "tables/proton_valence_quark_samples.dat";
    of_n_name_resample << "tables/neutron_valence_quark_samples.dat";
    std::ifstream of_p_test(of_p_name_resample.str().c_str(), std::ios::binary);
    std::ifstream of_n_test(of_n_name_resample.str().c_str(), std::ios::binary);
    if (!of_p_test.good() || !of_n_test.good()) {
        std::cout << "Generating files " << of_p_name_resample.str()
                  << " and " << of_n_name_resample.str() << std::endl;
        std::stringstream command;
        command << "./Metropolis.e " << 1;
        nucleon_system_status_ = std::system(command.str().c_str());
    } else {
        of_p_test.close();
        of_n_test.close();
    }

    std::ifstream of_p(of_p_name_resample.str().c_str(), std::ios::binary);
    std::ifstream of_n(of_n_name_resample.str().c_str(), std::ios::binary);
    if (!of_p.good() || !of_n.good()) {
        std::cout << "Can not generate " << of_p_name_resample.str() << " or/and "
                  << of_n_name_resample.str() << std::endl;
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
        proton_resample_quark_x_.push_back(x_array);
    }
    of_p.close();
    while (!of_n.eof()) {
        std::array<float, 3> x_array;
        for (int ii = 0; ii < 3; ii++) {
            float temp = 0.;
            of_n.read(reinterpret_cast<char*>(&temp), sizeof(float));
            x_array[ii] = temp;
        }
        neutron_resample_quark_x_.push_back(x_array);
    }
    of_n.close();
    int size = std::min(proton_resample_quark_x_.size(),
                        neutron_resample_quark_x_.size());
    return(size);
}

void Nucleon::resample_quark_momentum_fraction(std::vector<real> &xQuark,
                                               const int electric_charge,
                                               const real ecm) const {
    const real mq = PhysConsts::MQuarkValence;
    const real mp = PhysConsts::MProton;
    const int number_of_quarks = PhysConsts::NumValenceQuark;
    const real ybeam = acosh(ecm/(2.*mp));
    real total_energy = 0.;
    real E_proton = mp*cosh(ybeam);
    do {
        auto sample_idx = static_cast<int>(
            ran_gen_ptr_->rand_uniform()*number_of_valence_quark_resamples_);
        xQuark.clear();
        if (electric_charge == 1) {
            for (int i = 0; i < number_of_quarks; i++)
                xQuark.push_back(proton_resample_quark_x_[sample_idx][i]);
        } else  {
            for (int i = 0; i < number_of_quarks; i++)
                xQuark.push_back(proton_resample_quark_x_[sample_idx][i]);
        }
        total_energy = 0.;
        for (int i = 0; i < number_of_quarks; i++) {
            real rap_local = asinh(xQuark[i]*mp/mq*sinh(ybeam));
            real E_local = mq*cosh(rap_local);
            total_energy += E_local;
        }
    } while (total_energy > E_proton);
}


SpatialVec Nucleon::resample_valence_quark_position(real BG) const {
    // sample Gaussian distribution for the valence quark position
    // determine x,y,z coordinates of the quark (relative to the nucleon)

    BG = sqrt(BG)*PhysConsts::HBARC;
    real x = ran_gen_ptr_->rand_normal(0., BG);
    real y = ran_gen_ptr_->rand_normal(0., BG);
    real z = ran_gen_ptr_->rand_normal(0., BG);

    SpatialVec xq = {0.0, x, y, z};
    return(xq);
}

void Nucleon::resample_valence_quarks(real ecm, int direction, real charge, 
                   std::vector<double> xvec_q) {
    const int number_of_quarks = PhysConsts::NumValenceQuark;
    erase_quarks();
    std::vector<real> xQuark;
    resample_quark_momentum_fraction(xQuark, charge, ecm);
    for (int i = 0; i < number_of_quarks; i++) {
        SpatialVec xvec = {0.0, xvec_q[i*3], xvec_q[i*3+1], xvec_q[i*3+2]}; 
        std::shared_ptr<Quark> quark_ptr(new Quark(xvec, xQuark[i]));
        push_back_quark(quark_ptr);
    }
    accelerate_quarks(ecm, direction);
}


void Nucleon::readd_soft_parton_ball(real ecm, int direction,
                                     std::vector<double> xvec_q,
                                     real BG, MomentumVec soft_pvec,
                                     std::vector<std::shared_ptr<Quark>> valence_quark_list) {
    for (const auto & q_i: valence_quark_list) {
        auto quark_pvec = q_i->get_p();
        for (int i = 0; i < 4; i++) {
            soft_pvec[i] -= quark_pvec[i];
        }
    }
    real mass = PhysConsts::MQuarkValence;
    if (soft_pvec[0] > mass) {
        real rapidity = direction*acosh(soft_pvec[0]/mass);
        soft_pvec[3] = mass*sinh(rapidity);
        SpatialVec xvec;
        if (xvec_q.size()<10) {
            xvec = resample_valence_quark_position(BG);
        } else {
            xvec = {0.0, xvec_q[9], xvec_q[10], xvec_q[11]};
        }
        std::shared_ptr<Quark> quark_ptr(new Quark(xvec, soft_pvec));
        quark_ptr->set_rapidity(rapidity);
        push_back_quark(quark_ptr);
    }
}


void Nucleon::lorentz_contraction(real gamma) {
    for (auto &it: quark_list) {
        auto xvec = it->get_x();
        xvec[3] /= gamma;
        it->set_x(xvec);
    }
}


std::vector<double> Nucleon::output_quark_pos() {
    std::vector<double> quark_xvec;
    for (auto &it: quark_list) {
        auto xvec = it->get_x();
        for (unsigned int i=1; i<4; i++) {
            quark_xvec.push_back(xvec[i]);
        }
    }
    return (quark_xvec);
}


std::shared_ptr<Quark> Nucleon::get_a_valence_quark() {
    // return the quark with the minimum number of connections
    int minimum_connections = 1000;
    for (auto &iq: quark_list) {
        if (minimum_connections > iq->get_number_of_connections())
            minimum_connections = iq->get_number_of_connections();
    }
    std::shuffle(quark_list.begin(), quark_list.end(),
                 *ran_gen_ptr_->getRanGenerator());
    for (auto &iq: quark_list) {
        if (minimum_connections == iq->get_number_of_connections()) {
            iq->add_a_connection();
            return(iq);
        }
    }
    return(quark_list[0]);
}


std::shared_ptr<Quark> Nucleon::get_a_valence_quark_sub_mom(real sub_E) {
    // return the quark to subtract the momentum of hard collision
    std::shuffle(quark_list.begin(), quark_list.end(),
                 *ran_gen_ptr_->getRanGenerator());
    for (auto &iq: quark_list) {
        auto p_q = iq->get_p();
        if (p_q[0] > sub_E) return(iq);
    }
    return(quark_list[0]);
}


void Nucleon::erase_one_quark() {
    for (unsigned int idx = 0; idx < quark_list.size(); idx++) {
        if (quark_list[idx]->quark_is_subtracted()) {
            quark_list.erase(quark_list.begin()+idx);
            //std::cout << " Successfully erase one quark." << std::endl;
            break;
        }
    }
}

}
