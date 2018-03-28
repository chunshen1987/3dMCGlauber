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
    ran_gen_ptr = (
        std::unique_ptr<RandomUtil::Random>(new RandomUtil::Random(seed)));
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
    //projectile->output_nucleon_positions("projectile.dat");
    //target->output_nucleon_positions("target.dat");
}


int Glauber::make_collision_schedule() {
    collision_schedule.clear();
    auto d2 = (compute_NN_inelastic_cross_section(parameter_list.get_roots())
               /(M_PI*10.));  // in fm^2 
    auto proj_nucleon_list = projectile->get_nucleon_list();
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &iproj: (*proj_nucleon_list)) {
        auto proj_x = iproj.get_x();
        for (auto &itarg: (*targ_nucleon_list)) {
            auto targ_x = itarg.get_x();
            auto dij = (  (targ_x[1] - proj_x[1])*(targ_x[1] - proj_x[1])
                        + (targ_x[2] - proj_x[2])*(targ_x[2] - proj_x[2]));
            if (hit(dij, d2)) {
                create_a_collision_event(iproj, itarg);
            }
        }
    }
    auto Npart = (projectile->get_number_of_wounded_nucleons()
                  + target->get_number_of_wounded_nucleons());
    std::cout << "Npart = " << Npart << std::endl;
    return(collision_schedule.size());
}

bool Glauber::hit(real d2, real d2_in) {
    real G = 0.92;
    return(ran_gen_ptr->rand_uniform() < G*exp(-G*d2/d2_in));
}

void Glauber::create_a_collision_event(Nucleon &proj, Nucleon &targ) {
    real t_coll, z_coll;
    auto x1 = proj.get_x();
    auto p1 = proj.get_p();
    auto x2 = targ.get_x();
    auto p2 = targ.get_p();
    auto v1 = p1[3]/p1[0];
    auto v2 = p2[3]/p2[0];
    auto collided = get_collision_point(x1[0], x1[3], v1, x2[3], v2,
                                        t_coll, z_coll);
    if (collided) {
        SpatialVec x_coll = {t_coll, (x1[1] + x2[1])/2.,
                             (x1[2] + x2[2])/2., z_coll};
        CollisionEvent new_event(x_coll, proj, targ);
        collision_schedule.insert(new_event);
        proj.set_wounded(true);
        targ.set_wounded(true);
    }
}

bool Glauber::get_collision_point(real t, real z1, real v1, real z2, real v2,
                                  real &t_coll, real &z_coll) const {
    bool collided;
    if ((z2 - z1)*(v2 - v1) > 0.) {
        collided = false;
        t_coll = -1.;
        z_coll = -1.;
    } else {
        collided = true;
        real delta_t = std::abs((z2 - z1)/(v1 - v2));
        t_coll = t  + delta_t;
        z_coll = z1 + v1*delta_t;
    }
    return(collided);
}

real Glauber::compute_NN_inelastic_cross_section(real ecm) const {
    real s = ecm*ecm;
    real sigma_NN_total = 44.4 - 2.9*log(s) + 0.33*log(s)*log(s);
    real sigma_NN_inel  = (sigma_NN_total
                           - (11.4 - 1.52*log(s) + 0.13*log(s)*log(s)));
    return(sigma_NN_inel);
}

}
