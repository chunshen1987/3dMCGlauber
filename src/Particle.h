// Copyright @ Chun Shen 2018

#ifndef SRC_PARTICLE_H_
#define SRC_PARTICLE_H_

#include "data_structs.h"
#include <cassert>
#include <iostream>
#include <cmath>

namespace MCGlb {

class Particle {
 private:
    SpatialVec x;
    MomentumVec p;
    real mass;

 public:
    Particle() = default;

    Particle(SpatialVec x_in, MomentumVec p_in) {
        set_particle_variables(x_in, p_in);
    }
    bool operator == (const Particle &p1) const {
        return(x == p1.get_x() && p == p1.get_p() && mass == p1.get_mass());
    }

    Particle(SpatialVec x_in, MomentumVec p_in, real mass_in) {
        set_particle_variables(x_in, p_in, mass_in);
    }

    void set_particle_variables(SpatialVec x_in, MomentumVec p_in) {
        set_x(x_in); set_p(p_in);
        set_mass_with_momentum_vector();
    }

    void set_particle_variables(
                    SpatialVec x_in, MomentumVec p_in, real mass_in) {
        assert(std::abs(mass_in*mass_in
                        - ( p_in[0]*p_in[0] - p_in[1]*p_in[1]
                           - p_in[2]*p_in[2] - p_in[3]*p_in[3])) < 1e-15) ;
        set_x(x_in); set_p(p_in); set_mass(mass_in);
    }

    void set_x(SpatialVec x_in) {x = x_in;}
    SpatialVec get_x() const {return(x);}

    real get_rapidity() const {return(atanh(p[3]/p[0]));}

    void set_p(MomentumVec p_in) {
        assert(p_in[0] >= 0.);
        p = p_in;
        set_mass_with_momentum_vector();
    }
    void set_p_first(MomentumVec p_in) {
        assert(p_in[0] >= 0.);
        p = p_in;
    }
    MomentumVec get_p() const {return(p);}

    void set_mass(real mass_in) {
        assert(mass_in >= 0);
        mass = mass_in;
    }

    void set_mass_with_momentum_vector() {
        real m_inv = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
        assert(m_inv >= 0.);
        mass = sqrt(m_inv);
    }
    real get_mass() const {return(mass);}
};

}

#endif  // SRC_PARTICLE_H_
