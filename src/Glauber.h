// Copyright @ Chun Shen 2018

#ifndef SRC_GLAUBER_H_
#define SRC_GLAUBER_H_

#include "data_structs.h"
#include "Parameters.h"
#include "Random.h"
#include "Nucleus.h"
#include <memory>

namespace MCGlb {

class Glauber {
 private:
    const Parameters &parameter_list;
    std::unique_ptr<Nucleus> projectile;
    std::unique_ptr<Nucleus> target;

 public:
    Glauber() = default;
    Glauber(const MCGlb::Parameters &param_in);
    ~Glauber() {};

    void make_nuclei();
};

}

#endif   // SRC_GLAUBER_H_
