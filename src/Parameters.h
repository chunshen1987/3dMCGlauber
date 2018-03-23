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

    void set_b(real b_in);
    real get_b() {return(static_cast<real>(std::stod(get_param_val("b"))));}

    string get_projectle_nucleus_name() {return(get_param_val("Projectile"));}
    string get_target_nucleus_name() {return(get_param_val("Target"));}
};

}

#endif  // SRC_PARAMETERS_H_
