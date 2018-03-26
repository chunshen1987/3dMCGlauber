// Copyright @ Chun Shen 2018

#include "Glauber.h"
#include "data_structs.h"

#include <iostream>

using std::cout;
using std::endl;

namespace MCGlb {

Glauber::Glauber(const MCGlb::Parameters &param_in) :
    parameter_list(param_in) {
    parameter_list.print_parameter_list();
    int seed = parameter_list.get_seed();
    projectile = std::unique_ptr<Nucleus>(
            new Nucleus(parameter_list.get_projectle_nucleus_name(), seed));
    target = std::unique_ptr<Nucleus>(
            new Nucleus(parameter_list.get_target_nucleus_name(), seed));
    make_nuclei();
}

void Glauber::make_nuclei() {
    projectile->generate_nucleus_3d_configuration();
    target->generate_nucleus_3d_configuration();
    projectile->accelerate_nucleus(parameter_list.get_roots(), 1);
    target->accelerate_nucleus(parameter_list.get_roots(), -1);
    auto impact_b = parameter_list.get_b();
    SpatialVec proj_shift = {0., impact_b/2., 0., -projectile->get_z_max()};
    projectile->shift_nucleus(proj_shift);
    SpatialVec targ_shift = {0., -impact_b/2., 0., -target->get_z_min()};
    target->shift_nucleus(targ_shift);
    projectile->output_nucleon_positions("projectile.dat");
    target->output_nucleon_positions("target.dat");
}

}
