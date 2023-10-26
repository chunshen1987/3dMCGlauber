// Copyright @ Chun Shen 2018

#ifndef SRC_NUCLEON_H_
#define SRC_NUCLEON_H_

#include "Particle.h"
#include "Quark.h"
#include "RandomUlty.h"
#include <vector>
#include <memory>
#include <algorithm>
#include <string>
namespace MCGlb {

class Nucleon : public Particle {
 private:
    std::vector<std::shared_ptr<Quark>> quark_list;
    int collided_times = 0;
    int total_connected_times_ = 0;
    int electric_charge_ = 0;
    bool wounded_ = false;
    bool baryon_used = false;
    bool remnant_set_ = false;
    bool remnant_carry_baryon_number_ = false;
    bool collided_ = false;
    bool subtracted_ = false;
    std::vector<std::weak_ptr<Nucleon>> collide_with;
    std::vector<std::weak_ptr<Nucleon>> connected_with;
    std::vector<int> connected_times_;
    MomentumVec remnant_p_ = {0.0, 0.0, 0.0, 0.0};
    SpatialVec remnant_x_frez_ = {0.0, 0.0, 0.0, 0.0};
    std::vector< std::array<float, 3> > proton_resample_quark_x_;
    std::vector< std::array<float, 3> > neutron_resample_quark_x_;
    int number_of_valence_quark_resamples_;
    int nucleon_system_status_;
    static int random_value_;
    MomentumVec fermi_momentum_ = {0.0, 0.0, 0.0, 0.0};

 public:
    Nucleon() = default;
    Nucleon(SpatialVec x_in, MomentumVec p_in);

    Nucleon(SpatialVec x_in, MomentumVec p_in, real mass_in) {
        set_particle_variables(x_in, p_in, mass_in);
    }

    ~Nucleon();

    void set_electric_charge(int charge) {electric_charge_ = charge;}
    int get_electric_charge() const {return(electric_charge_);}
    int re_readin_valence_quark_samples();

    int get_number_of_quarks() const {return(quark_list.size());}
    void push_back_quark(std::shared_ptr<Quark> q) {quark_list.push_back(q);}
    void erase_quarks() {
        quark_list.erase(quark_list.begin(), quark_list.end());
    }
    void erase_one_quark();
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
    void lorentz_contraction(real gamma);

    std::shared_ptr<Quark> get_a_valence_quark(int ran_seed);
    std::shared_ptr<Quark> get_a_valence_quark_sub_mom(real sub_E,
                                                       int ran_seed);

    std::vector<double> output_quark_pos();

    bool is_remnant_set() const {return(remnant_set_);}
    void set_remnant(bool remnant) {remnant_set_ = remnant;}

    bool is_hard_collided() const {return(collided_);}
    void set_hard_collided(bool collided_index) {collided_ = collided_index;}

    bool nucleon_is_subtracted() const {return(subtracted_);}
    void set_hard_subtracted(bool subtracted_index) {
        subtracted_ = subtracted_index;
    }

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

    void resample_quark_momentum_fraction(
            std::vector<real> &xQuark, const int electric_charge,
            const real ecm,
            std::shared_ptr<RandomUtil::Random> nucleon_ran_gen_ptr) const;
    void resample_valence_quarks(
            real ecm, int direction, real charge, std::vector<double> xvec_q,
            std::shared_ptr<RandomUtil::Random> nucleon_ran_gen_ptr);
    void readd_soft_parton_ball(real ecm, int direction, std::vector<double> xvec_q,
                                real BG_, MomentumVec soft_pvec,
                                std::vector<std::shared_ptr<Quark>> valence_quark_list,
                                std::shared_ptr<RandomUtil::Random> nucleon_ran_gen_ptr);

    SpatialVec resample_valence_quark_position(real BG_, 
               std::shared_ptr<RandomUtil::Random> nucleon_ran_gen_ptr) const;

    void set_fermi_momentum(real px, real py, real pz) {
        fermi_momentum_[1] = px;
        fermi_momentum_[2] = py;
        fermi_momentum_[3] = pz;
    }
    MomentumVec get_fermi_momentum() const {return(fermi_momentum_);}

    void set_remnant_x_frez(SpatialVec x_in) {remnant_x_frez_ = x_in;}
    SpatialVec get_remnant_x_frez() const {return(remnant_x_frez_);}
    static int get_random_gen(int i) {return random_value_%i;}
    static void set_random_gen(int i) {random_value_ = i;}
};

}

#endif  // SRC_NUCLEON_H_
