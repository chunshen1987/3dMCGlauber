// Copyright @ Chun Shen 2018

#include "Nucleon.h"
#include "PhysConsts.h"
#include<cmath>

namespace MCGlb {

Nucleon::Nucleon(SpatialVec x_in, MomentumVec p_in) {
    set_particle_variables(x_in, p_in);
}

Nucleon::~Nucleon() {
    quark_list.clear();
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

void Nucleon::lorentz_contraction(real gamma) {
    for (auto &it: quark_list) {
        auto xvec = it->get_x();
        xvec[3] /= gamma;
        it->set_x(xvec);
    }
}


std::shared_ptr<Quark> Nucleon::get_a_valence_quark() {
    // return the quark with the minimum number of connections
    int minimum_connections = 1000;
    for (auto &iq: quark_list) {
        if (minimum_connections > iq->get_number_of_connections())
            minimum_connections = iq->get_number_of_connections();
    }
    std::random_shuffle(quark_list.begin(), quark_list.end());
    for (auto &iq: quark_list) {
        if (minimum_connections == iq->get_number_of_connections()) {
            iq->add_a_connection();
            return(iq);
        }
    }
    return(quark_list[0]);
}


}
