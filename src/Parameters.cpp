// Copyright (C) 2018 Chun Shen

#include "Parameters.h"
#include <cassert>

namespace MCGlb {

void Parameters::set_b(real b_in) {
    assert(b_in >= 0.);
    set_parameter("b", b_in);
}

real Parameters::get_b() const {
    real b = static_cast<real>(get_param_double("b"));
    assert(b >= 0.);
    return(b);
}

void Parameters::set_b_max(real b_in) {
    assert(b_in >= 0.);
    set_parameter("b_max", b_in);
}

real Parameters::get_b_max() const {
    real b = static_cast<real>(get_param_double("b_max"));
    assert(b >= 0.);
    return(b);
}

void Parameters::set_b_min(real b_in) {
    assert(b_in >= 0.);
    set_parameter("b_min", b_in);
}


real Parameters::get_b_min() const {
    real b = static_cast<real>(get_param_double("b_min"));
    assert(b >= 0.);
    return(b);
}

    
int Parameters::get_use_energy_dependent_cross_section() const {
    int flag = get_param_int("useEnergyDependentCrossSection");
    assert(flag == 0 || flag == 1);
    return(flag);
}


int Parameters::get_use_quarks() const {
    int flag = get_param_int("useQuarks");
    assert(flag >= 0 && flag < 3);
    return(flag);
}


int Parameters::get_gaussian_wounding() const {
    int flag = get_param_int("gaussianWounding");
    assert(flag == 0 || flag == 1);
    return(flag);
}
    
real Parameters::get_roots() const {
    real roots = static_cast<real>(get_param_double("roots"));
    assert(roots > 0.);
    return(roots);
}
    
int Parameters::get_QCD_string_production_mode() const {
    int flag = get_param_int("QCD_string_production_mode");
    assert(flag >= 0 && flag < 5);
    return(flag);
}
    
real Parameters::get_string_tension() const {
    real sigma = static_cast<real>(get_param_double("string_tension"));
    assert(sigma >= 0.);
    return(sigma);
}

}
