// Copyright @ Chun Shen 2018

#ifndef SRC_NUCLEON_H_
#define SRC_NUCLEON_H_

#include "Particle.h"
#include "Quark.h"
#include <vector>
#include <memory>
#include <algorithm>

namespace MCGlb {

class Nucleon : public Particle {
 private:
    std::vector<std::shared_ptr<Quark>> quark_list;
    int collided_times = 0;
    bool wounded = false;
    bool baryon_used = false;
    bool remnant_set_ = false;
    std::vector<std::weak_ptr<Nucleon>> collide_with;
    std::vector<std::weak_ptr<Nucleon>> connected_with;
    MomentumVec remnant_p_ = {0.0, 0.0, 0.0, 0.0};
    SpatialVec remnant_x_frez_ = {0.0, 0.0, 0.0, 0.0};


 public:
    Nucleon() = default;
    Nucleon(SpatialVec x_in, MomentumVec p_in);

    Nucleon(SpatialVec x_in, MomentumVec p_in, real mass_in) {
        set_particle_variables(x_in, p_in, mass_in);
    }

    ~Nucleon();

    int get_number_of_quarks() const {return(quark_list.size());}
    void push_back_quark(std::shared_ptr<Quark> q) {quark_list.push_back(q);}
    std::vector<std::shared_ptr<Quark>> get_quark_list() {return(quark_list);}

    bool is_wounded() const {return(wounded);}
    bool baryon_was_used() const {return(baryon_used);}
    void set_wounded(bool hit) {wounded = hit;}
    void set_baryon_used(bool hit) {baryon_used = hit;}

    void increment_collided_times() {collided_times++;}
    int get_collided_times() const {return(collided_times);}
    void add_collide_nucleon(std::weak_ptr<Nucleon> collide_nucleon) {
        collide_with.push_back(collide_nucleon);
    }
    std::vector<std::weak_ptr<Nucleon>>* get_collide_nucleon_list() {
        return(&collide_with);
    }

    int get_number_of_connections() {return(connected_with.size());}
    void add_connected_nucleon(std::weak_ptr<Nucleon> connected_nucleon) {
        connected_with.push_back(connected_nucleon);
    }
    bool is_connected_with(std::shared_ptr<Nucleon> targ);
    void accelerate_quarks(real ecm, int direction);
    void lorentz_contraction(real gamma);

    std::shared_ptr<Quark> get_a_valence_quark();

    bool is_remnant_set() const {return(remnant_set_);}
    void set_remnant(bool remnant) {remnant_set_ = remnant;}

    void set_remnant_p(MomentumVec p_in) {remnant_p_ = p_in;}
    MomentumVec get_remnant_p() const {return(remnant_p_);}
    void substract_momentum_from_remnant(MomentumVec p_q) {
        for (int i = 0; i < 4; i++)
            remnant_p_[i] -= p_q[i];
    }

    void set_remnant_x_frez(SpatialVec x_in) {remnant_x_frez_ = x_in;}
    SpatialVec get_remnant_x_frez() const {return(remnant_x_frez_);}
};

}

#endif  // SRC_NUCLEON_H_
