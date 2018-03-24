// Copyright @ Chun Shen 2018

#include "Glauber.h"
#include "data_structs.h"

#include <iostream>

using MCGlb::Nucleon;
using MCGlb::Quark;

using std::cout;
using std::endl;

Glauber::Glauber(const MCGlb::Parameters &param_in) :
    parameter_list(param_in){
    parameter_list.print_parameter_list();
    int seed = parameter_list.get_seed();
    ran_gen_ptr = (
            std::unique_ptr<RandomUtil::Random>(new RandomUtil::Random(seed)));
    for (int i = 0; i < 10; i++) {
        cout << ran_gen_ptr->rand_uniform() << endl;
    }
}
