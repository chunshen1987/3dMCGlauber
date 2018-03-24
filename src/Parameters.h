// Copyright (C) 2018 Chun Shen

#ifndef SRC_PARAMETERS_H_
#define SRC_PARAMETERS_H_

#include <string>
#include "ParametersMap.h"
#include "data_structs.h"

using std::string;
using std::stod;
using std::stoi;

namespace MCGlb {

class Parameters : public ParametersMap {
 public:
    Parameters() = default;
    ~Parameters() {};

    int get_seed() const {return(stoi(get_param_val("seed")));}

    string get_projectle_nucleus_name() const {
        return(get_param_val("Projectile"));
    }
    string get_target_nucleus_name() const {return(get_param_val("Target"));}

    void set_b(real b_in);
    real get_b() const {return(static_cast<real>(stod(get_param_val("b"))));}
    int get_use_energy_dependent_cross_section() const {
        return(stoi(get_param_val("useEnergyDependentCrossSection")));
    }

    int get_use_quarks() const {return(stoi(get_param_val("useQuarks")));}
    int get_gaussian_wounding() const {
        return(stoi(get_param_val("gaussianWounding")));
    }

    real get_roots() const {
        return(static_cast<real>(stod(get_param_val("roots"))));
    }
};

}

#endif  // SRC_PARAMETERS_H_
