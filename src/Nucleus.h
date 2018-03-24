// Copyright @ Chun Shen 2018

#ifndef SRC_NUCLEUS_H_
#define SRC_NUCLEUS_H_

#include "Nucleon.h"
#include "Quark.h"
#include <vector>
#include <string>

namespace MCGlb {

class Nucleus {
 private:
    std::string name;
    int density_function_type;
    int A;
    int Z;
    real WS_w;
    real WS_a;
    real WS_R;
    real WS_rho;
    std::vector<Nucleon> Nucleon_list;

 public:
    Nucleus() = default;
};

}

#endif  // SRC_NUCLEUS_H_
