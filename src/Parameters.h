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

    void set_b(real b_in);
    real get_b() const;
    void set_b_max(real b_in);
    real get_b_max() const;
    void set_b_min(real b_in);
    real get_b_min() const;

    int get_use_energy_dependent_cross_section() const;

    int get_use_quarks() const;

    int get_gaussian_wounding() const;

    real get_roots() const;

    int get_QCD_string_production_mode() const;
    
    real get_string_tension() const;

};

}

#endif  // SRC_PARAMETERS_H_
