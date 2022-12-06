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
    std::vector<QCDString> QCD_string_list;
    std::vector<QCDString> remnant_string_list_;
    std::vector<std::vector<real>> QCD_string_output_arr_;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr_;
    bool sample_valence_quark;
    bool fluct_Nstrings_per_NN_collision_;
    real remnant_energy_loss_fraction_;

    real impact_b;
    real yloss_param_slope;
    real yloss_param_a;
    real yloss_param_b;

    real ybeam;

    int system_status_;

    real sigma_eff_;
    real nucleon_width_;

 public:
    Glauber() = default;
    Glauber(const MCGlb::Parameters &param_in,
            shared_ptr<RandomUtil::Random> ran_gen);
    ~Glauber() {};

    void make_nuclei();
    real get_impact_parameter() const {return(impact_b);}

    int make_collision_schedule();
    bool hit(real d2) const;

    int get_Npart() const;

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
    int decide_produce_string_num(shared_ptr<CollisionEvent> event_ptr) const;

    real get_tau_form(const int string_evolution_mode) const;
    void get_tau_form_and_moversigma(const int string_evolution_mode,
                                     const real y_in_lrf,
                                     real &tau_form, real &m_over_sigma,
                                     real &y_loss);

    real sample_rapidity_loss_shell(real y_init) const;
    real sample_rapidity_loss_from_the_LEXUS_model(const real y_init) const;
    real sample_rapidity_loss_from_parametrization(const real y_init) const;
    real sample_rapidity_loss_from_piecewise_parametrization(
                                                const real y_init) const;
    real sample_rapidity_loss_from_parametrization_with_fluct(
                                                const real y_init) const;

    real sample_junction_rapidity_right(real y_left, real y_right) const;
    real sample_junction_rapidity_left(real y_left, real y_right) const;
    real sample_junction_rapidity_uniformed(real y_left, real y_right) const;

    //! This function performs string production between each nucleon pair
    int perform_string_production();
    void produce_remnant_strings();

    //! This function propagate individual nucleon inside the nucleus by dt
    void propagate_nuclei(real dt);
    void propagate_nucleon(shared_ptr<Nucleon> n_i, real dt);
    void update_momentum(shared_ptr<Nucleon> n_i, real y_shift);
    void update_momentum_quark(shared_ptr<Quark> q_i, real y_shift);
    //! This function updates the collision schedule
    void update_collision_schedule(shared_ptr<CollisionEvent> event_happened);

    void prepare_output_QCD_strings();
    void computeCenterOfMass(real &x_o, real &y_o);
    void output_QCD_strings(std::string filename, const real Npart,
                            const real Ncoll, const real Nstrings,
                            const real b, const unsigned int seed);
    void output_spectators(std::string filename);
    void outputParticipants(std::string filename);

    real get_sig_eff(const real siginNN);
};

}

#endif   // SRC_GLAUBER_H_
