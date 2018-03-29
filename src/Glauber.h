// Copyright @ Chun Shen 2018

#ifndef SRC_GLAUBER_H_
#define SRC_GLAUBER_H_

#include "data_structs.h"
#include "Parameters.h"
#include "Random.h"
#include "Nucleus.h"
#include "CollisionEvent.h"
#include <memory>
#include <set>

namespace MCGlb {

class Glauber {
 private:
    const Parameters &parameter_list;
    std::unique_ptr<Nucleus> projectile;
    std::unique_ptr<Nucleus> target;
    std::set<CollisionEvent> collision_schedule;
    std::weak_ptr<RandomUtil::Random> ran_gen_ptr;

 public:
    Glauber() = default;
    Glauber(const MCGlb::Parameters &param_in,
            std::shared_ptr<RandomUtil::Random> ran_gen);
    ~Glauber() {};

    void make_nuclei();

    int make_collision_schedule();
    bool hit(real d2, real d2_in);

    int get_Npart();

    //! This function creates a new collision event between two nucleons
    void create_a_collision_event(Nucleon &proj, Nucleon &targ);
    bool get_collision_point(real t, real z1, real v1, real z2, real v2,
                             real &t_coll, real &z_coll) const;

    real compute_NN_inelastic_cross_section(real ecm) const;
};

}

#endif   // SRC_GLAUBER_H_
