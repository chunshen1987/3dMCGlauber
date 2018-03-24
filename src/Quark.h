// Copyright @ Chun Shen 2018

#ifndef SRC_QUARK_H_
#define SRC_QUARK_H_

#include "Particle.h"
#include <cassert>

namespace MCGlb {

class Quark : public Particle {
 private:
    real pdf_x;

 public:
    Quark() = default;
    Quark(SpatialVec x_in, MomentumVec p_in) {
        set_particle_variables(x_in, p_in);
    }

    Quark(SpatialVec x_in, MomentumVec p_in, real mass_in) {
        set_particle_variables(x_in, p_in, mass_in);
    }

    void set_pdf_x(real x_in) {
        assert(x_in > 0.);
        assert(x_in <= 1.);
        pdf_x = x_in;
    }
    real get_pdf_x() {return(pdf_x);}

};

}


#endif  // SRC_NUCLEON_H_
