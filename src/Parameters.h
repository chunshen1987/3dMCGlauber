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

    real get_lambdaB() const;

    int get_QCD_string_production_mode() const;
    int get_QCD_string_evolution_mode() const;

    bool get_baryon_junctions() const; // if 0 do assume baryon number at string ends, if 1 transport baryon number according to cosh(y*/2)

    bool get_only_event_statistics() const;

};

}

#endif  // SRC_PARAMETERS_H_
