// Copyright @ Chun Shen 2018

#include "Glauber.h"
#include "data_structs.h"

using MCGlb::Nucleon;
using MCGlb::Quark;

Glauber::Glauber(const MCGlb::Parameters &param_in) :
    parameter_list(param_in) {
    parameter_list.print_parameter_list();
}
