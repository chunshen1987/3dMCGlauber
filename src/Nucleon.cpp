// Copyright @ Chun Shen 2018

#include "Nucleon.h"

#include <cmath>

#include "PhysConsts.h"

namespace MCGlb {

Nucleon::Nucleon(
    SpatialVec x_in, MomentumVec p_in,
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr) {
    ran_gen_ptr_ = ran_gen_ptr;
    set_particle_variables(x_in, p_in);
}

Nucleon::Nucleon(
    SpatialVec x_in, MomentumVec p_in, real mass_in,
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr) {
    ran_gen_ptr_ = ran_gen_ptr;
    set_particle_variables(x_in, p_in, mass_in);
}

Nucleon::~Nucleon() { quark_list.clear(); }

bool Nucleon::is_connected_with(std::shared_ptr<Nucleon> targ) {
    bool connected = false;
    for (auto &it : connected_with) {
        if (*(it.lock()) == *targ) {
            connected = true;
            break;
        }
    }
    return (connected);
}

int Nucleon::get_number_of_connections(std::shared_ptr<Nucleon> targ) const {
    int n_connections = 0;
    for (unsigned int idx = 0; idx < connected_with.size(); idx++) {
        if (*(connected_with[idx].lock()) == *targ) {
            n_connections = connected_times_[idx];
            break;
        }
    }
    return (n_connections);
}

void Nucleon::accelerate_quarks(real ecm, int direction) {
    const real mq = PhysConsts::MQuarkValence;
    const real mN = get_mass();
    const real ybeam = acosh(ecm / (2. * mN));
    for (auto &it : quark_list) {
        // real rap_local = direction*asinh(it->get_pdf_x()
        //                                  *sqrt(ecm*ecm/(4.*mq*mq) - 1.));
        real rap_local =
            direction * asinh(it->get_pdf_x() * mN / mq * sinh(ybeam));
        it->set_rapidity(rap_local);
        MomentumVec p_in = {
            mq * cosh(rap_local), 0.0, 0.0, mq * sinh(rap_local)};
        it->set_p(p_in);
    }
}

void Nucleon::lorentz_contraction(real gamma) {
    for (auto &it : quark_list) {
        auto xvec = it->get_x();
        xvec[3] /= gamma;
        it->set_x(xvec);
    }
}

std::shared_ptr<Quark> Nucleon::get_a_valence_quark() {
    // return the quark with the minimum number of connections
    int minimum_connections = 1000;
    for (auto &iq : quark_list) {
        if (minimum_connections > iq->get_number_of_connections())
            minimum_connections = iq->get_number_of_connections();
    }
    std::shuffle(
        quark_list.begin(), quark_list.end(), *ran_gen_ptr_->getRanGenerator());
    for (auto &iq : quark_list) {
        if (minimum_connections == iq->get_number_of_connections()) {
            iq->add_a_connection();
            return (iq);
        }
    }
    return (quark_list[0]);
}

std::shared_ptr<Quark> Nucleon::get_a_close_valence_quark(real xq, real yq) {
    // return the quark with the minimum number of connections
    //, and choose the quark with closest distrance
    std::vector<std::pair<real, int>> vec;
    int minimum_connections = 1000;
    int idex = 0;
    for (auto &iq : quark_list) {
        if (minimum_connections > iq->get_number_of_connections())
            minimum_connections = iq->get_number_of_connections();
        auto q_xvec = iq->get_x();
        real dis2 = (q_xvec[1] - xq) * (q_xvec[1] - xq)
                    + (q_xvec[2] - yq) * (q_xvec[2] - yq);
        vec.push_back({dis2, idex});
        idex++;
    }
    std::sort(vec.begin(), vec.end());
    for (auto &qidex : vec) {
        auto iq = quark_list[qidex.second];
        if (minimum_connections == iq->get_number_of_connections()) {
            iq->add_a_connection();
            return (iq);
        }
    }
    return (quark_list[0]);
}

}  // namespace MCGlb
