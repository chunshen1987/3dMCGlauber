// Copyright @ Chun Shen 2018

#ifndef SRC_GLAUBER_H_
#define SRC_GLAUBER_H_

#include "data_structs.h"
#include "Parameters.h"

class Glauber {
 private:
    const MCGlb::Parameters &parameter_list;

 public:
    Glauber() = default;
    Glauber(const MCGlb::Parameters &param_in);
    ~Glauber() {};
};

#endif   // SRC_GLAUBER_H_
