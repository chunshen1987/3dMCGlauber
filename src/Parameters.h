// Copyright (C) 2018 Chun Shen

#ifndef SRC_PARAMETERS_H_
#define SRC_PARAMETERS_H_

#include <string>
#include "ParametersMap.h"
#include "data_structs.h"

namespace MCGlb {

class Parameters : public ParametersMap {
 private:
    std::string parameter_filename;
    real b;

    std::string projectile_nucleus_name;
    std::string target_nucleus_name;

 public:
    Parameters() = default;
    ~Parameters() {};

    void set_parameters();

    void set_parameter_filename(std::string input_filename) {
        parameter_filename = input_filename;
    }

    void set_b(real b_in);
    real get_b();

    void set_projectile_nucleus(std::string nucleus_name) {
        projectile_nucleus_name = nucleus_name;
    }
    std::string get_projectle_nucleus_name() {return(projectile_nucleus_name);}
    
    void set_target_nucleus(std::string nucleus_name) {
        target_nucleus_name = nucleus_name;
    }
    std::string get_target_nucleus_name() {return(target_nucleus_name);}

};

}

#endif  // SRC_PARAMETERS_H_
