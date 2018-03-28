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

}
