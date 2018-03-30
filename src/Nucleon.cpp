// Copyright @ Chun Shen 2018

#include "Nucleon.h"

namespace MCGlb {

Nucleon::Nucleon(SpatialVec x_in, MomentumVec p_in) {
    collided_times = 0;
    wounded = false;
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

}
