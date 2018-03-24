// Copyright @ Chun Shen 2018

#ifndef SRC_GLAUBER_H_
#define SRC_GLAUBER_H_

#include "data_structs.h"
#include "Parameters.h"
#include "Random.h"
#include <memory>

class Glauber {
 private:
    const MCGlb::Parameters &parameter_list;
    std::unique_ptr<RandomUtil::Random> ran_gen_ptr;

 public:
    Glauber() = default;
    Glauber(const MCGlb::Parameters &param_in);
    ~Glauber() {};
};

#endif   // SRC_GLAUBER_H_
