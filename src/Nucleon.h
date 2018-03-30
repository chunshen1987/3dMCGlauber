// Copyright @ Chun Shen 2018

#ifndef SRC_NUCLEON_H_
#define SRC_NUCLEON_H_

#include "Particle.h"
#include "Quark.h"
#include <vector>
#include <memory>

namespace MCGlb {

class Nucleon : public Particle {
 private:
    std::vector<std::shared_ptr<Quark>> quark_list;
    int collided_times;
    bool wounded;
    std::vector<std::weak_ptr<Nucleon>> collide_with;
    std::vector<std::weak_ptr<Nucleon>> connected_with;

 public:
    Nucleon() = default;
    Nucleon(SpatialVec x_in, MomentumVec p_in);

    Nucleon(SpatialVec x_in, MomentumVec p_in, real mass_in) {
        set_particle_variables(x_in, p_in, mass_in);
    }

    ~Nucleon();

    int get_number_of_quarks() const {return(quark_list.size());}
    bool is_wounded() const {return(wounded);}
    void set_wounded(bool hit) {wounded = hit;}

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
};

}

#endif  // SRC_NUCLEON_H_
