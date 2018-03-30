// Copyright @ Chun Shen 2018

#ifndef SRC_GLAUBER_H_
#define SRC_GLAUBER_H_

#include "data_structs.h"
#include "Parameters.h"
#include "Random.h"
#include "Nucleus.h"
#include "CollisionEvent.h"
#include "QCDString.h"
#include <memory>
#include <string>
#include <set>
#include <vector>

using std::shared_ptr;

namespace MCGlb {

class Glauber {
 private:
    const Parameters &parameter_list;
    std::unique_ptr<Nucleus> projectile;
    std::unique_ptr<Nucleus> target;
    std::set<shared_ptr<CollisionEvent>, compare_collision_time> collision_schedule;
    std::vector<shared_ptr<QCDString>> QCD_string_list;
    std::weak_ptr<RandomUtil::Random> ran_gen_ptr;

 public:
    Glauber() = default;
    Glauber(const MCGlb::Parameters &param_in,
            shared_ptr<RandomUtil::Random> ran_gen);
    ~Glauber() {};

    void make_nuclei();

    int make_collision_schedule();
    bool hit(real d2, real d2_in);

    int get_Npart();

    //! This function creates a new collision event between two nucleons
    void create_a_collision_event(shared_ptr<Nucleon> proj,
                                  shared_ptr<Nucleon> targ);
    bool get_collision_point(real t, real z1, real v1, real z2, real v2,
                             real &t_coll, real &z_coll) const;

    real compute_NN_inelastic_cross_section(real ecm) const;

    //! this function decides which of those binary collisions will produce
    //! QCD strings
    int decide_QCD_strings_production();

    //! this function determines whether a given binary collision event
    //! will produce a string
    bool decide_produce_string(shared_ptr<CollisionEvent> event_ptr) const;

    //! This function performs string production between each nucleon pair
    int perform_string_production();
    //! This function propagate individual nucleon inside the nucleus by dt
    void propagate_nuclei(real dt);
    void propagate_nucleon(shared_ptr<Nucleon> n_i, real dt);
    void update_momentum(shared_ptr<Nucleon> n_i, real y_shift);
    //! This function updates the collision schedule
    void update_collision_schedule(shared_ptr<CollisionEvent> event_happened);
    
    void output_QCD_strings(std::string filename);
};

}

#endif   // SRC_GLAUBER_H_
