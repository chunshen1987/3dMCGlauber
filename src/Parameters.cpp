// Copyright (C) 2018 Chun Shen

#include "Parameters.h"
#include <cassert>

namespace MCGlb {

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

    
int Parameters::get_use_quarks() const {
    int flag = get_param_int("useQuarks");
    assert(flag >= 0 && flag < 3);
    return(flag);
}
    

real Parameters::get_quarks_Q2() const {
    real Q2 = static_cast<real>(get_param_double("Q2"));
    assert(Q2 >= 0.);
    return(Q2);
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


int Parameters::get_QCD_string_evolution_mode() const {
    int flag = get_param_int("evolve_QCD_string_mode");
    assert(flag > 0 && flag < 5);
    return(flag);
}


bool Parameters::get_only_event_statistics() const {
    int flag = get_param_int("only_event_statistics");
    if (flag == 0) {
        return(false);
    } else {
        return(true);
    }
}
    
}
