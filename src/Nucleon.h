// Copyright @ Chun Shen 2018

#ifndef SRC_NUCLEON_H_
#define SRC_NUCLEON_H_

#include "Particle.h"
#include "Random.h"
#include "Quark.h"
#include <vector>
#include <memory>
#include <algorithm>

namespace MCGlb {

class Nucleon : public Particle {
 private:
    std::vector<std::shared_ptr<Quark>> quark_list;
    int collided_times = 0;
    int total_connected_times_ = 0;
    int electric_charge_ = 0;
    int baryon_number_=1;
    bool wounded_ = false;
    bool baryon_used = false;
    bool remnant_set_ = false;
    bool remnant_carry_baryon_number_ = false;
    std::vector<std::weak_ptr<Nucleon>> collide_with;
    std::vector<std::weak_ptr<Nucleon>> connected_with;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr_;
    std::vector<int> connected_times_;
    MomentumVec remnant_p_ = {0.0, 0.0, 0.0, 0.0};
    SpatialVec remnant_x_frez_ = {0.0, 0.0, 0.0, 0.0};
    MomentumVec fermi_momentum_ = {0.0, 0.0, 0.0, 0.0};

 public:
    Nucleon() = default;
    Nucleon(SpatialVec x_in, MomentumVec p_in,
            std::shared_ptr<RandomUtil::Random> ran_gen_ptr);

    Nucleon(SpatialVec x_in, MomentumVec p_in, real mass_in,
            std::shared_ptr<RandomUtil::Random> ran_gen_ptr);

    ~Nucleon();


    void set_baryon_number(int baryon_number) {baryon_number_ = baryon_number;}
    int get_baryon_number() const {return(baryon_number_);}

    void set_electric_charge(int charge) {electric_charge_ = charge;}
    int get_electric_charge() const {return(electric_charge_);}

    int get_number_of_quarks() const {return(quark_list.size());}
    void push_back_quark(std::shared_ptr<Quark> q) {quark_list.push_back(q);}
    std::vector<std::shared_ptr<Quark>> get_quark_list() {return(quark_list);}

    bool is_wounded() const {return(wounded_);}
    bool baryon_was_used() const {return(baryon_used);}
    void set_wounded(bool hit) {wounded_ = hit;}
    void set_baryon_used(bool hit) {baryon_used = hit;}

    void increment_collided_times() {collided_times++;}
    int get_collided_times() const {return(collided_times);}
    void add_collide_nucleon(std::weak_ptr<Nucleon> collide_nucleon) {
        collide_with.push_back(collide_nucleon);
    }
    std::vector<std::weak_ptr<Nucleon>>* get_collide_nucleon_list() {
        return(&collide_with);
    }

    void add_connected_nucleon(std::weak_ptr<Nucleon> connected_nucleon) {
        connected_with.push_back(connected_nucleon);
    }

    void add_num_connections(int N_connections) {
        connected_times_.push_back(N_connections);
        total_connected_times_ += N_connections;
    }

    int get_number_of_connections() const {
        return(total_connected_times_);
    }

    int get_number_of_connections(const int idx) const {
        return(connected_times_[idx]);
    }

    int get_number_of_connections(std::shared_ptr<Nucleon> targ) const;

    bool is_connected_with(std::shared_ptr<Nucleon> targ);
    void accelerate_quarks(real ecm, int direction);
    void accelerate_quarks_in_dipole(real ecm, int direction);
    void lorentz_contraction(real gamma);

    std::shared_ptr<Quark> get_a_valence_quark();
    std::shared_ptr<Quark> get_a_close_valence_quark(real xq, real yq);

    bool is_remnant_set() const {return(remnant_set_);}
    void set_remnant(bool remnant) {remnant_set_ = remnant;}

    bool is_remnant_carry_baryon_number() const {
        return(remnant_carry_baryon_number_);
    }
    void set_remnant_carry_baryon_number(bool remnant) {
        remnant_carry_baryon_number_ = remnant;
    }

    void set_remnant_p(MomentumVec p_in) {remnant_p_ = p_in;}
    MomentumVec get_remnant_p() const {return(remnant_p_);}
    void substract_momentum_from_remnant(MomentumVec p_q) {
        for (int i = 0; i < 4; i++)
            remnant_p_[i] -= p_q[i];
    }

    void set_fermi_momentum(real px, real py, real pz) {
        fermi_momentum_[1] = px;
        fermi_momentum_[2] = py;
        fermi_momentum_[3] = pz;
    }
    MomentumVec get_fermi_momentum() const {return(fermi_momentum_);}

    void set_remnant_x_frez(SpatialVec x_in) {remnant_x_frez_ = x_in;}
    SpatialVec get_remnant_x_frez() const {return(remnant_x_frez_);}
};

}

#endif  // SRC_NUCLEON_H_
