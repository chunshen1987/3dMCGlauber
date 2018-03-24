// Copyright @ Chun Shen 2018

#ifndef SRC_NUCLEON_H_
#define SRC_NUCLEON_H_

#include "Particle.h"
#include "Quark.h"
#include <vector>

namespace MCGlb {

class Nucleon : public Particle {
 private:
    std::vector<Quark> quark_list;

 public:
    Nucleon() = default;
    Nucleon(SpatialVec x_in, MomentumVec p_in) {
        set_particle_variables(x_in, p_in);
    }

    Nucleon(SpatialVec x_in, MomentumVec p_in, real mass_in) {
        set_particle_variables(x_in, p_in, mass_in);
    }

    ~Nucleon();

    int get_number_of_quarks() {return(quark_list.size());}

};

}

#endif  // SRC_NUCLEON_H_
