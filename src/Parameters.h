// Copyright (C) 2018 Chun Shen

#ifndef SRC_PARAMETERS_H_
#define SRC_PARAMETERS_H_

#include <string>
#include "ParametersMap.h"
#include "data_structs.h"

using std::string;

namespace MCGlb {

class Parameters : public ParametersMap {
 public:
    Parameters() = default;
    ~Parameters() {};

    int get_seed() const {return(get_param_int("seed"));}

    string get_projectle_nucleus_name() const {
        return(get_param_val("Projectile"));
    }
    string get_target_nucleus_name() const {return(get_param_val("Target"));}

    void set_b_max(real b_in);
    real get_b_max() const;
    void set_b_min(real b_in);
    real get_b_min() const;

    int get_use_quarks() const;
    real get_quarks_Q2() const;

    real get_roots() const;
    real get_UPC_root_low_cut() const;
    real get_UPC_root_up_cut() const;

    bool use_roots_distribution() const;
    bool use_roots_cut() const;

    bool use_E_dependent_LB() const;
    real get_lambdaB() const;
    real get_CB() const;
    real get_shadowing_factor() const;

    int get_QCD_string_production_mode() const;
    int get_QCD_string_evolution_mode() const;
    int get_rapidity_loss_method() const;

    // if False do assume baryon number at string ends
    // if True transport baryon number according to cosh(y*/2)
    bool get_baryon_junctions() const;

    bool get_only_event_statistics() const;
    bool get_cached_tabels() const;
    bool set_remnant_x_ori() const;

    bool get_fluct_Nstrings_per_NN_collision() const;
    real get_remnant_energy_loss_fraction() const;

    real get_yloss_param_slope() const;
    real get_yloss_param_alpha1() const;
    real get_yloss_param_alpha2() const;
    real get_yloss_param_fluct_var_RHIC() const;
    real get_yloss_param_fluct_var_LHC() const;

    real get_tau_form_mean() const;
    real get_tau_form_fluct_gamma_beta() const;

    real get_BG() const;

    real get_BG_proj() const;
    real get_BG_targ() const;

    bool nucleon_configuration_from_file() const;
    int lightNucleusOption() const;
};

}

#endif  // SRC_PARAMETERS_H_
