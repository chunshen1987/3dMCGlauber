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

    void set_b(real b_in) {set_parameter("b", b_in);}
    real get_b() const {
        return(static_cast<real>(get_param_double("b")));
    }
    void set_b_max(real b_in) {set_parameter("b_max", b_in);}
    real get_b_max() const {
        return(static_cast<real>(get_param_double("b_max")));
    }
    void set_b_min(real b_in) {set_parameter("b_min", b_in);}
    real get_b_min() const {
        return(static_cast<real>(get_param_double("b_min")));
    }

    int get_use_energy_dependent_cross_section() const {
        return(get_param_int("useEnergyDependentCrossSection"));
    }

    int get_use_quarks() const {return(get_param_int("useQuarks"));}
    int get_gaussian_wounding() const {
        return(get_param_int("gaussianWounding"));
    }

    real get_roots() const {
        return(static_cast<real>(get_param_double("roots")));
    }

    int get_QCD_string_production_mode() const {
        return(get_param_int("QCD_string_production_mode"));
    }
    
    real get_string_tension() const {
        return(static_cast<real>(get_param_double("string_tension")));
    }

};

}

#endif  // SRC_PARAMETERS_H_
