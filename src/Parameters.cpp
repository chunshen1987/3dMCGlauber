// Copyright (C) 2018 Chun Shen

#include "Parameters.h"
#include <cassert>

namespace MCGlb {


real Parameters::getParam(std::string paramName, real defaultValue) const {
    return(static_cast<real>(get_param_double(paramName, defaultValue)));
}


void Parameters::set_b_max(real b_in) {
    assert(b_in >= 0.);
    set_parameter("b_max", b_in);
}


real Parameters::get_b_max() const {
    real b = static_cast<real>(get_param_double("b_max", 25.));
    assert(b >= 0.);
    return(b);
}


void Parameters::set_b_min(real b_in) {
    assert(b_in >= 0.);
    set_parameter("b_min", b_in);
}


real Parameters::get_b_min() const {
    real b = static_cast<real>(get_param_double("b_min", 0.));
    assert(b >= 0.);
    return(b);
}


real Parameters::get_d_min() const {
    real b = static_cast<real>(get_param_double("d_min", 0.));
    assert(b >= 0.);
    return(b);
}


int Parameters::get_use_quarks() const {
    int flag = get_param_int("useQuarks", 1);
    assert(flag >= 0 && flag < 3);
    return(flag);
}


real Parameters::get_quarks_Q2() const {
    real Q2 = static_cast<real>(get_param_double("Q2", 1.));
    assert(Q2 >= 0.);
    return(Q2);
}


real Parameters::get_roots() const {
    real roots = static_cast<real>(get_param_double("roots", 200.));
    assert(roots > 0.);
    return(roots);
}


real Parameters::get_lambdaB() const {
    real lambdaB = static_cast<real>(get_param_double("lambdaB", 0.));
    assert(lambdaB >= 0.);
    return(lambdaB);
}


real Parameters::get_lambdaBs() const {
    real lambdaBs = static_cast<real>(get_param_double("lambdaBs", 1.));
    assert(lambdaBs >= 0.);
    return(lambdaBs);
}


real Parameters::get_baryon_in_string_prob() const {
    real prob = static_cast<real>(get_param_double("baryonInStringProb", 1.));
    assert(prob >= 0. && prob <= 1.);
    return(prob);
}


bool Parameters::get_cached_tabels() const {
    int flag = get_param_int("cache_tables", 1);
    if (flag == 1)
        return(true);
    else
        return(false);
}


bool Parameters::get_fluct_Nstrings_per_NN_collision() const {
    int flag = get_param_int("fluct_Nstrings_per_NN_collision", 1);
    if (flag == 1)
        return(true);
    else
        return(false);
}


double Parameters::get_remnant_energy_loss_fraction() const {
    real frac = static_cast<real>(
            get_param_double("remnant_energy_loss_fraction", 0.5));
    assert(frac >= 0. && frac <= 1.);
    return(frac);
}


int Parameters::get_QCD_string_production_mode() const {
    int flag = get_param_int("QCD_string_production_mode", 1);
    assert(flag >= 0 && flag < 5);
    return(flag);
}


int Parameters::get_QCD_string_evolution_mode() const {
    int flag = get_param_int("evolve_QCD_string_mode", 4);
    assert(flag > 0 && flag < 5);
    return(flag);
}


int Parameters::get_rapidity_loss_method() const {
    int flag = get_param_int("rapidity_loss_method", 3);
    assert(flag > 0 && flag < 5);
    return(flag);
}


bool Parameters::get_only_event_statistics() const {
    int flag = get_param_int("only_event_statistics", 0);
    if (flag == 0) {
        return(false);
    } else {
        return(true);
    }
}


bool Parameters::get_batch_density_output() const {
    int flag = get_param_int("batch_density_output", 0);
    if (flag == 0) {
        return(false);
    } else {
        return(true);
    }
}


bool Parameters::get_initialEst_output() const {
    int flag = get_param_int("outputInitialEst", 0);
    if (flag == 0) {
        return(false);
    } else {
        return(true);
    }
}


bool Parameters::get_batch_2Ddensity_output() const {
    int flag = get_param_int("batch_2Ddensity_output", 0);
    if (flag == 0) {
        return(false);
    } else {
        return(true);
    }
}


bool Parameters::get_batch_eccentricity_output() const {
    int flag = get_param_int("batch_eccentricity_output", 0);
    if (flag == 0) {
        return(false);
    } else {
        return(true);
    }
}


bool Parameters::get_baryon_junctions() const {
    int flag = get_param_int("baryon_junctions", 0);
    if (flag == 0) {
        return(false);
    } else {
        return(true);
    }
}


real Parameters::get_shadowing_factor() const {
    real shadowing = static_cast<real>(
                    get_param_double("shadowing_factor", 1.0));
    assert(shadowing >= 0.);
    assert(shadowing <= 1.);
    return(shadowing);
}


real Parameters::get_yloss_param_slope() const {
    real slope = static_cast<real>(
                    get_param_double("yloss_param_slope", 1.32));
    assert(slope >= 0.);
    //assert(slope <= 1.);
    return(slope);
}


real Parameters::get_yloss_param_alpha1() const {
    real a = static_cast<real>(
                    get_param_double("yloss_param_alpha1", 1.8));
    assert(a >= 1.);
    return(a);
}


real Parameters::get_yloss_param_alpha2() const {
    real a = static_cast<real>(
                    get_param_double("yloss_param_alpha2", 0.35));
    assert(a >= 0.);
    assert(a <= 1.);
    return(a);
}


real Parameters::get_yloss_param_fluct_var_LHC() const {
    real a = static_cast<real>(
                    get_param_double("yloss_param_fluct_var_LHC", 0.6));
    assert(a >= 0.);
    return(a);
}


real Parameters::get_yloss_param_fluct_var_RHIC() const {
    real a = static_cast<real>(
                    get_param_double("yloss_param_fluct_var_RHIC", 0.6));
    assert(a >= 0.);
    return(a);
}


real Parameters::get_tau_form_mean() const {
    real tau_form_mean = static_cast<real>(
                    get_param_double("tau_form_mean", 0.5));
    assert(tau_form_mean > 0.);
    return(tau_form_mean);
}


real Parameters::get_tau_form_fluct_gamma_beta() const {
    real tau_form_beta = static_cast<real>(
            get_param_double("tau_form_fluct_gamma_beta", 1.0));
    return(tau_form_beta);
}


bool Parameters::nucleon_configuration_from_file() const {
    int flag = get_param_int("nucleon_configuration_from_file", 0);
    if (flag == 0) {
        return(false);
    } else {
        return(true);
    }
}


real Parameters::get_BG() const {
    real BG = static_cast<real>(get_param_double("BG", 5.));
    assert(BG > 0.);
    return(BG);
}

}
