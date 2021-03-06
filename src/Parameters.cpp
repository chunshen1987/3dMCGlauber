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

real Parameters::get_lambdaB() const {
    real lambdaB = static_cast<real>(get_param_double("lambdaB"));
    assert(lambdaB >= 0.);
    return(lambdaB);
}


bool Parameters::get_cached_tabels() const {
    int flag = get_param_int("cache_tables");
    if (flag == 1)
        return(true);
    else
        return(false);
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


int Parameters::get_rapidity_loss_method() const {
    int flag = get_param_int("rapidity_loss_method");
    assert(flag > 0 && flag < 4);
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

bool Parameters::get_baryon_junctions() const {
    int flag = get_param_int("baryon_junctions");
    if (flag == 0) {
        return(false);
    } else {
        return(true);
    }
}

real Parameters::get_shadowing_factor() const {
    real shadowing = static_cast<real>(get_param_double("shadowing_factor"));
    assert(shadowing >= 0.);
    assert(shadowing <= 1.);
    return(shadowing);
}


real Parameters::get_yloss_param_slope() const {
    real slope = static_cast<real>(get_param_double("yloss_param_slope"));
    assert(slope >= 0.);
    //assert(slope <= 1.);
    return(slope);
}


real Parameters::get_yloss_param_alpha1() const {
    real a = static_cast<real>(get_param_double("yloss_param_alpha1"));
    assert(a >= 1.);
    return(a);
}


real Parameters::get_yloss_param_alpha2() const {
    real a = static_cast<real>(get_param_double("yloss_param_alpha2"));
    assert(a >= 0.);
    assert(a <= 1.);
    return(a);
}


real Parameters::get_yloss_param_fluct_var_LHC() const {
    real a = static_cast<real>(get_param_double("yloss_param_fluct_var_LHC"));
    assert(a >= 0.);
    return(a);
}


real Parameters::get_yloss_param_fluct_var_RHIC() const {
    real a = static_cast<real>(get_param_double("yloss_param_fluct_var_RHIC"));
    assert(a >= 0.);
    return(a);
}


real Parameters::get_tau_form_mean() const {
    real tau_form_mean = static_cast<real>(get_param_double("tau_form_mean"));
    assert(tau_form_mean > 0.);
    return(tau_form_mean);
}


real Parameters::get_tau_form_fluct_gamma_beta() const {
    real tau_form_beta = static_cast<real>(
            get_param_double("tau_form_fluct_gamma_beta"));
    return(tau_form_beta);
}

}
