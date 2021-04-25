// Copyright @ Chun Shen 2018

#include "Glauber.h"
#include "data_structs.h"
#include "PhysConsts.h"

#include <string>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <memory>

using std::cout;
using std::endl;
using std::shared_ptr;


namespace MCGlb {

Glauber::Glauber(const MCGlb::Parameters &param_in,
                 shared_ptr<RandomUtil::Random> ran_gen) :
    parameter_list(param_in) {
    sample_valence_quark = false;
    if (parameter_list.get_use_quarks() > 0) {
        sample_valence_quark = true;
        if (!parameter_list.get_cached_tabels()) {
            system_status_ = std::system(
                    "rm -fr tables/proton_valence_quark_samples*");
            system_status_ = std::system(
                    "rm -fr tables/neutron_valence_quark_samples*");
        }
    }
    projectile = std::unique_ptr<Nucleus>(
            new Nucleus(parameter_list.get_projectle_nucleus_name(), ran_gen,
                        sample_valence_quark));
    target = std::unique_ptr<Nucleus>(
            new Nucleus(parameter_list.get_target_nucleus_name(), ran_gen,
                        sample_valence_quark));
    if (sample_valence_quark) {
        projectile->set_valence_quark_Q2(parameter_list.get_quarks_Q2());
        target->set_valence_quark_Q2(parameter_list.get_quarks_Q2());
    }
    ran_gen_ptr_ = ran_gen;

    yloss_param_slope = parameter_list.get_yloss_param_slope();
    real alpha1 = parameter_list.get_yloss_param_alpha1();
    real alpha2 = parameter_list.get_yloss_param_alpha2();
    real alpha = alpha2/alpha1;
    yloss_param_a = alpha/(1. - alpha);
    yloss_param_b = alpha2/yloss_param_a;

    ybeam = acosh(parameter_list.get_roots()/(2.*PhysConsts::MProton));

    real siginNN = compute_NN_inelastic_cross_section(
                                            parameter_list.get_roots());
    sigma_eff_ = get_sig_eff(siginNN);
    cout << "sqrt{s} = " << parameter_list.get_roots() << " GeV, "
         << "siginNN = " << siginNN << " mb" << endl;
}

void Glauber::make_nuclei() {
    projectile->generate_nucleus_3d_configuration();
    target->generate_nucleus_3d_configuration();
    projectile->accelerate_nucleus(parameter_list.get_roots(), 1);
    target->accelerate_nucleus(parameter_list.get_roots(), -1);

    // sample impact parameters
    auto b_max = parameter_list.get_b_max();
    auto b_min = parameter_list.get_b_min();
    impact_b = sqrt(b_min*b_min +
            (b_max*b_max - b_min*b_min)*ran_gen_ptr_->rand_uniform());
    SpatialVec proj_shift = {0., impact_b/2., 0.,
                             -projectile->get_z_max() - 1e-15};
    projectile->shift_nucleus(proj_shift);
    SpatialVec targ_shift = {0., -impact_b/2., 0.,
                             -target->get_z_min() + 1e-15};
    target->shift_nucleus(targ_shift);
    //projectile->output_nucleon_positions("projectile.dat");
    //target->output_nucleon_positions("target.dat");
}


int Glauber::make_collision_schedule() {
    collision_schedule.clear();
    auto proj_nucleon_list = projectile->get_nucleon_list();
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &iproj: (*proj_nucleon_list)) {
        auto proj_x = iproj->get_x();
        for (auto &itarg: (*targ_nucleon_list)) {
            auto targ_x = itarg->get_x();
            auto dij = (  (targ_x[1] - proj_x[1])*(targ_x[1] - proj_x[1])
                        + (targ_x[2] - proj_x[2])*(targ_x[2] - proj_x[2]));
            if (hit(dij)) {
                create_a_collision_event(iproj, itarg);
                projectile->add_a_participant(iproj);
                target->add_a_participant(itarg);
                iproj->set_wounded(true);
                itarg->set_wounded(true);
                iproj->add_collide_nucleon(itarg);
                itarg->add_collide_nucleon(iproj);
            }
        }
    }
    return(collision_schedule.size());
}

int Glauber::get_Npart() const {
    int Npart = (projectile->get_number_of_wounded_nucleons()
                  + target->get_number_of_wounded_nucleons());
    return(Npart);
}

bool Glauber::hit(real d2) const {
    //real G = 0.92;  // from Glassando
    //return(ran_gen_ptr_->rand_uniform() < G*exp(-G*d2/d2_in));
    const real T_nn = (exp(-d2/(4.*nucleon_width_*nucleon_width_))
                       /(4.*M_PI*nucleon_width_*nucleon_width_));
    const real hit_treshold = 1. - exp(-sigma_eff_*T_nn);
    return (ran_gen_ptr_->rand_uniform() < hit_treshold);
}

void Glauber::create_a_collision_event(shared_ptr<Nucleon> proj,
                                       shared_ptr<Nucleon> targ) {
    real t_coll, z_coll;
    auto x1 = proj->get_x();
    auto p1 = proj->get_p();
    auto x2 = targ->get_x();
    auto p2 = targ->get_p();
    auto v1 = p1[3]/p1[0];
    auto v2 = p2[3]/p2[0];
    auto collided = get_collision_point(x1[0], x1[3], v1, x2[3], v2,
                                        t_coll, z_coll);
    if (collided) {
        SpatialVec x_coll = {t_coll, (x1[1] + x2[1])/2.,
                             (x1[2] + x2[2])/2., z_coll};
        shared_ptr<CollisionEvent> event_ptr(
                                    new CollisionEvent(x_coll, proj, targ));
        if (proj->is_connected_with(targ)) {
            auto form_N_strings = proj->get_number_of_connections(targ);
            event_ptr->set_produced_n_strings(form_N_strings);
        }
        collision_schedule.insert(event_ptr);
    }
}

bool Glauber::get_collision_point(real t, real z1, real v1, real z2, real v2,
                                  real &t_coll, real &z_coll) const {
    bool collided = false;
    real delta_t = (z2 - z1)/(v1 - v2);
    if (delta_t > 0.) {
        collided = true;
        t_coll = t  + delta_t;
        z_coll = z1 + v1*delta_t;
    } else {
        t_coll = -1.;
        z_coll = -1.;
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


int Glauber::decide_produce_string_num(
                                shared_ptr<CollisionEvent> event_ptr) const {
    int form_n_string = 0;
    auto proj = event_ptr->get_proj_nucleon_ptr().lock();
    auto targ = event_ptr->get_targ_nucleon_ptr().lock();
    int minimum_allowed_connections = 1;
    if (sample_valence_quark) {
        // assume P(N) = 1/2^N distribution
        auto rand = ran_gen_ptr_->rand_uniform();
        int N = std::min(proj->get_number_of_quarks(),
                         targ->get_number_of_quarks());
        minimum_allowed_connections = static_cast<int>(
                -log((1. - rand*(1. - pow(2., -N))))/log(2.)) + 1;
    }

    if (   proj->get_number_of_connections() < minimum_allowed_connections
        && targ->get_number_of_connections() < minimum_allowed_connections) {
        form_n_string = std::min(
            minimum_allowed_connections - proj->get_number_of_connections(),
            minimum_allowed_connections - targ->get_number_of_connections());
    } else if (
           proj->get_number_of_connections() < minimum_allowed_connections
        || targ->get_number_of_connections() < minimum_allowed_connections) {
        form_n_string = 1;
    } else {
        int n_connects = (  proj->get_number_of_connections()
                          + targ->get_number_of_connections());
        real shadowing_factor = parameter_list.get_shadowing_factor();
        real production_prob = (
            (1. - shadowing_factor)
            *exp(-shadowing_factor*(n_connects
                                    - 2.*minimum_allowed_connections)));
        if (ran_gen_ptr_->rand_uniform() < production_prob) {
            form_n_string = 1;
        }
    }
    return(form_n_string);
}


int Glauber::decide_QCD_strings_production() {
    if (sample_valence_quark) {
        projectile->sample_valence_quarks_inside_nucleons(
                                    parameter_list.get_roots(), 1);
        target->sample_valence_quarks_inside_nucleons(
                                    parameter_list.get_roots(), -1);
        projectile->add_soft_parton_ball(parameter_list.get_roots(), 1);
        target->add_soft_parton_ball(parameter_list.get_roots(), -1);
    }
    std::vector<shared_ptr<CollisionEvent>> collision_list;
    for (auto &it: collision_schedule)
        collision_list.push_back(it);  // collision list is time ordered

    const auto QCD_string_production_mode =
                            parameter_list.get_QCD_string_production_mode();
    if (QCD_string_production_mode == 1) {
        // randomly ordered strings
        std::random_shuffle(collision_list.begin(), collision_list.end());
    } else if (QCD_string_production_mode == 2) {
        // anti-time ordered strings
        std::reverse(collision_list.begin(), collision_list.end());
    }

    int number_of_strings = 0;
    bool finished = false;
    while (!finished) {
        finished = true;
        for (auto &ievent: collision_list) {
            auto form_N_strings = decide_produce_string_num(ievent);
            auto proj = ievent->get_proj_nucleon_ptr();
            auto targ = ievent->get_targ_nucleon_ptr();
            if (form_N_strings > 0) {
                number_of_strings += form_N_strings;
                proj.lock()->add_connected_nucleon(targ);
                proj.lock()->add_num_connections(form_N_strings);
                targ.lock()->add_connected_nucleon(proj);
                targ.lock()->add_num_connections(form_N_strings);
                ievent->set_produced_n_strings(form_N_strings);
            } else {
                if (   proj.lock()->get_number_of_connections() == 0
                    || targ.lock()->get_number_of_connections() == 0)
                        finished = false;
            }
        }
    }
    return(number_of_strings);
}


//! This function propagate individual nucleon inside the nucleus by dt
void Glauber::propagate_nuclei(real dt) {
    for (auto &n_i: (*projectile->get_nucleon_list()))
        propagate_nucleon(n_i, dt);
    for (auto &n_i: (*target->get_nucleon_list()))
        propagate_nucleon(n_i, dt);
}


void Glauber::propagate_nucleon(shared_ptr<Nucleon> n_i, real dt) {
    auto x_vec = n_i->get_x();
    auto p_vec = n_i->get_p();
    for (int i = 0; i < 4; i++) {
        real v_i = p_vec[i]/p_vec[0];
        x_vec[i] += dt*v_i;
    }
    n_i->set_x(x_vec);
}


void Glauber::update_momentum_quark(shared_ptr<Quark> q_i, real y_shift) {
    auto pvec = q_i->get_p();
    auto mass = q_i->get_mass();
    auto y_i = q_i->get_rapidity();
    auto y_f = y_i + y_shift;
    q_i->set_rapidity(y_f);
    pvec[0] = mass*cosh(y_f);
    pvec[3] = mass*sinh(y_f);
    q_i->set_p(pvec);
}


void Glauber::update_momentum(shared_ptr<Nucleon> n_i, real y_shift) {
    auto pvec = n_i->get_p();
    auto y_i = n_i->get_rapidity();
    auto y_f = y_i + y_shift;
    pvec[0] = PhysConsts::MProton*cosh(y_f);
    pvec[3] = PhysConsts::MProton*sinh(y_f);
    n_i->set_p(pvec);
}


real Glauber::get_tau_form(const int string_evolution_mode) const {
    auto tau_form_mean = parameter_list.get_tau_form_mean();
    real tau_form = tau_form_mean;      // [fm]
    if (string_evolution_mode == 2) {
        // tau_form fluctuates with a gamma distribution
        tau_form = tau_form_mean/ran_gen_ptr_->rand_gamma_dis();
        tau_form = std::min(std::max(0.1, tau_form), 6*tau_form_mean);
    }
    return(tau_form);
}


void Glauber::get_tau_form_and_moversigma(const int string_evolution_mode,
                                          const real y_in_lrf,
                                          real &tau_form, real &m_over_sigma,
                                          real &y_loss) {
    tau_form = get_tau_form(string_evolution_mode);      // [fm]
    m_over_sigma = 1.0;  // [fm]
    y_loss = 0.;
    if (string_evolution_mode == 1) {
        // fixed rapidity loss
        y_loss = (
            acosh(tau_form*tau_form/(2.*m_over_sigma*m_over_sigma) + 1.));
        // maximum 90% of y_in_lrf
        y_loss = std::min(y_loss, y_in_lrf*0.9);
        m_over_sigma = tau_form/sqrt(2.*(cosh(y_loss) - 1.));
    } else if (string_evolution_mode == 2) {
        // both tau_form and sigma fluctuate
        y_loss = sample_rapidity_loss_shell(y_in_lrf);
        m_over_sigma = tau_form/sqrt(2.*(cosh(y_loss) - 1.));
    } else if (string_evolution_mode == 3) {
        // only tau_form fluctuates
        y_loss = sample_rapidity_loss_shell(y_in_lrf);
        tau_form = m_over_sigma*sqrt(2.*(cosh(y_loss) - 1.));
    } else if (string_evolution_mode == 4) {
        // only m_over_sigma fluctuates
        y_loss = sample_rapidity_loss_shell(y_in_lrf);
        m_over_sigma = tau_form/sqrt(2.*(cosh(y_loss) - 1.));
    }
}

// add flag if baryon used to proj and target
// - check it - only put in baryon if not used yet, then set it to used
int Glauber::perform_string_production() {
    QCD_string_list.clear();
    remnant_string_list_.clear();
    const auto string_evolution_mode = (
                    parameter_list.get_QCD_string_evolution_mode());
    const auto baryon_junctions = parameter_list.get_baryon_junctions();

    // sqrt(parameter_list.get_roots()); // ~s^{-1/4} 
    real lambdaB = parameter_list.get_lambdaB();
    lambdaB = std::min(1., lambdaB);

    //cout << lambdaB <<endl;

    real t_current = 0.0;
    int number_of_collided_events = 0;
    while (collision_schedule.size() > 0) {
        auto first_event = (*collision_schedule.begin());
        if (!first_event->is_valid()) {
            collision_schedule.erase((*collision_schedule.begin()));
            continue;
        }
        number_of_collided_events++;
        real dt = first_event->get_collision_time() - t_current;
        assert(dt > 0.);
        t_current = first_event->get_collision_time();
        propagate_nuclei(dt);
        if (first_event->get_num_strings() == 0) {
            collision_schedule.erase((*collision_schedule.begin()));
            continue;
        }
        for (int istring = 0; istring < first_event->get_num_strings();
                istring++) {
            auto x_coll = first_event->get_collision_position();
            auto proj = first_event->get_proj_nucleon_ptr().lock();
            auto targ = first_event->get_targ_nucleon_ptr().lock();
            real y_in_lrf = std::abs(proj->get_rapidity()
                                     - targ->get_rapidity())/2.;
            std::shared_ptr<Quark> proj_q;
            std::shared_ptr<Quark> targ_q;
            if (sample_valence_quark) {
                proj_q = proj->get_a_valence_quark();
                if (proj_q->get_number_of_connections() == 1) {
                    // first time pick-up the valence quark
                    // we need to substract the valence quark energy-momentum
                    // from the nucleon remnant energy-momentum vector
                    auto p_q = proj_q->get_p();
                    proj->substract_momentum_from_remnant(p_q);
                }
                targ_q = targ->get_a_valence_quark();
                if (targ_q->get_number_of_connections() == 1) {
                    // first time pick-up the valence quark
                    // we need to substract the valence quark energy-momentum
                    // from the nucleon remnant energy-momentum vector
                    auto p_q = targ_q->get_p();
                    targ->substract_momentum_from_remnant(p_q);
                }
                y_in_lrf = std::abs(  proj_q->get_rapidity()
                                    - targ_q->get_rapidity())/2.;
            }
            real tau_form = 0.5;      // [fm]
            real m_over_sigma = 1.0;  // [fm]
            real y_loss = 0.;
            get_tau_form_and_moversigma(string_evolution_mode, y_in_lrf,
                                        tau_form, m_over_sigma, y_loss);
            // set variables in case of no baryon junction transport
            bool has_baryon_left = false;
            bool has_baryon_right = false;
            if (!sample_valence_quark) {
                QCDString qcd_string(x_coll, tau_form, proj, targ, m_over_sigma,
                                     has_baryon_right, has_baryon_left);
                QCD_string_list.push_back(qcd_string);
            } else {
                auto proj_q_xvec = proj_q->get_x();
                auto targ_q_xvec = targ_q->get_x();
                SpatialVec x_coll_q = {
                    x_coll[0],
                    x_coll[1] + (proj_q_xvec[1] + targ_q_xvec[1])/2.,
                    x_coll[2] + (proj_q_xvec[2] + targ_q_xvec[2])/2.,
                    x_coll[3]};
                QCDString qcd_string(x_coll_q, tau_form, proj, targ,
                                     proj_q, targ_q, m_over_sigma,
                                     has_baryon_right, has_baryon_left);
                QCD_string_list.push_back(qcd_string);
            }
            real y_shift = y_loss;
            if (!sample_valence_quark) {
                update_momentum(proj, -y_shift);
                update_momentum(targ,  y_shift);
            } else {
                // shift the nucleon rapidity to avoid identical collision when
                // update the collision schedule
                update_momentum(proj, -1e-3*proj->get_rapidity());
                update_momentum(targ,  1e-3*targ->get_rapidity());
                update_momentum_quark(proj_q, -y_shift);
                update_momentum_quark(targ_q,  y_shift);
            }
        }
        update_collision_schedule(first_event);
        collision_schedule.erase((*collision_schedule.begin()));
    }

    // randomize the QCD_string_list and assign the baryon charge to
    // the strings
    std::vector<unsigned int> random_idx;
    unsigned int Nstrings = QCD_string_list.size();
    unsigned int Npart_proj = projectile->get_number_of_wounded_nucleons();
    unsigned int Npart_targ = target->get_number_of_wounded_nucleons();
    unsigned int total_length = Nstrings + Npart_proj + Npart_targ;
    for (unsigned int idx = 0; idx < total_length; idx++)
        random_idx.push_back(idx);
    std::random_shuffle(random_idx.begin(), random_idx.end());
    for (auto &idx: random_idx) {
        if (idx < Nstrings) {
            // put baryon of the projectile in the selected string
            auto proj = QCD_string_list[idx].get_proj();
            if (!proj.lock()->baryon_was_used()) {
                proj.lock()->set_baryon_used(true);
                QCD_string_list[idx].set_has_baryon_right(true);
            }
        } else if (idx < Nstrings + Npart_proj) {
            // put baryon of the projectile in the projectile remnant
            auto proj = projectile->get_participant(idx - Nstrings);
            auto p_i = proj.lock()->get_remnant_p();
            auto mass = 0.;
            if (std::abs(p_i[3]) < p_i[0]) {
                // a time-like beam remnant
                mass = sqrt(p_i[0]*p_i[0] - p_i[3]*p_i[3]);
            }
            if (!proj.lock()->baryon_was_used() && mass > 0.1) {
                proj.lock()->set_baryon_used(true);
                proj.lock()->set_remnant_carry_baryon_number(true);
            }
        }
    }
    std::random_shuffle(random_idx.begin(), random_idx.end());
    for (auto &idx: random_idx) {
        if (idx < Nstrings) {
            // put baryon of the target in the selected string
            auto targ = QCD_string_list[idx].get_targ();
            if (!targ.lock()->baryon_was_used()) {
                targ.lock()->set_baryon_used(true);
                QCD_string_list[idx].set_has_baryon_left(true);
            }
        } else if (idx > Nstrings + Npart_proj - 1) {
            // put baryon of the target in the target remnant
            auto targ = target->get_participant(idx - Nstrings - Npart_proj);
            auto p_i = targ.lock()->get_remnant_p();
            auto mass = 0.;
            if (std::abs(p_i[3]) < p_i[0]) {
                // a time-like beam remnant
                mass = sqrt(p_i[0]*p_i[0] - p_i[3]*p_i[3]);
            }
            if (!targ.lock()->baryon_was_used() && mass > 0.1) {
                targ.lock()->set_baryon_used(true);
                targ.lock()->set_remnant_carry_baryon_number(true);
            }
        }
    }

    // set baryons' rapidities
    for (auto &it: QCD_string_list) {
        it.evolve_QCD_string();
        if (!baryon_junctions) {
            // set baryon rapidities to string endpoint rapidities
            // if no junction transport is used
            it.set_final_baryon_rapidities(it.get_y_f_left(),
                                           it.get_y_f_right());
        } else {
            // sample HERE if baryon should be moved
            real y_baryon_right = 0.;
            if (it.get_has_baryon_right()) {
                if (ran_gen_ptr_->rand_uniform() < lambdaB) {
                    // y_baryon_right = sample_junction_rapidity_right(
                    //              it->get_y_i_left(), it->get_y_i_right());
                    y_baryon_right = sample_junction_rapidity_right(
                                    it.get_y_f_left(), it.get_y_f_right());
                } else {
                    // One should use the very initial rapidities of the
                    // colliding nucleons to determine the junction rapidity
                    // this may, however, lead to junctions lying outside the
                    // rapidity range of the final string which may cause
                    // problems in hydro and could be seen as unphysical...
                    // Another thing to discuss is the energy dependence.
                    // The cross section for baryon junctions stopping has
                    // a different root-s dependence than the usual inelastic
                    // cross section
                    // - so in principle it needs to be treated separately
                    // altogether... not sure how yet
                    y_baryon_right = it.get_y_f_right();
                }
            }
            real y_baryon_left = 0.;
            if (it.get_has_baryon_left()) {
                if (ran_gen_ptr_->rand_uniform() < lambdaB) {
                    //y_baryon_left = sample_junction_rapidity_left(
                    //              it->get_y_i_left(), it->get_y_i_right());
                    y_baryon_left = sample_junction_rapidity_left(
                                    it.get_y_f_left(), it.get_y_f_right());
                } else {
                    y_baryon_left = it.get_y_f_left();
                }
            }
            it.set_final_baryon_rapidities(y_baryon_left, y_baryon_right);
        }
    }


    // collision remnant is assigned to the last string which connects the
    // colliding nucleons or quarks
    for (std::vector<QCDString>::reverse_iterator it = QCD_string_list.rbegin();
            it != QCD_string_list.rend(); ++it) {
        if (sample_valence_quark) {
            // record the freeze-out space-time position for the remnant of
            // the collding nucleons at their last produced strings
            auto proj_n = it->get_proj();
            if (!proj_n.lock()->is_remnant_set()) {
                proj_n.lock()->set_remnant(true);
                auto x_frez = proj_n.lock()->get_x();
                proj_n.lock()->set_remnant_x_frez(x_frez);
            }
            auto targ_n = it->get_targ();
            if (!targ_n.lock()->is_remnant_set()) {
                targ_n.lock()->set_remnant(true);
                auto x_frez = targ_n.lock()->get_x();
                targ_n.lock()->set_remnant_x_frez(x_frez);
            }

            // set flags for quark remnants at their last connected strings
            auto proj_q = it->get_proj_q();
            if (!proj_q.lock()->is_remnant_set()) {
                proj_q.lock()->set_remnant(true);
                it->set_has_remnant_right(true);
            }
            auto targ_q = it->get_targ_q();
            if (!targ_q.lock()->is_remnant_set()) {
                targ_q.lock()->set_remnant(true);
                it->set_has_remnant_left(true);
            }
        } else {
            auto proj_n = it->get_proj();
            if (!proj_n.lock()->is_remnant_set()) {
                proj_n.lock()->set_remnant(true);
                it->set_has_remnant_right(true);
            }
            auto targ_n = it->get_targ();
            if (!targ_n.lock()->is_remnant_set()) {
                targ_n.lock()->set_remnant(true);
                it->set_has_remnant_left(true);
            }
        }
    }

    produce_remnant_strings();
    return(number_of_collided_events);
}


void Glauber::produce_remnant_strings() {
    // create strings for the beam remnants
    const auto string_evolution_mode = (
                    parameter_list.get_QCD_string_evolution_mode());
    real tau_form = 0.5;
    real m_over_sigma = 1.0;  // [fm]
    real y_loss = 0.;
    auto proj_nucleon_list = projectile->get_nucleon_list();
    for (auto &iproj: (*proj_nucleon_list)) {
        if (iproj->is_wounded()) {
            auto x_i = iproj->get_remnant_x_frez();

            auto p_i = iproj->get_remnant_p();
            if (p_i[0] <= 0.) continue;
            auto y_rem = ybeam;
            if (std::abs(p_i[3]) < p_i[0]) {
                // a time-like beam remnant
                y_rem = 0.5*log((p_i[0] + p_i[3])/(p_i[0] - p_i[3]));
            }
            auto cosh_y_rem = cosh(y_rem);
            auto m_rem = p_i[0]/cosh_y_rem;
            p_i[1] = 0.;
            p_i[2] = 0.;
            p_i[3] = m_rem*sinh(y_rem);
            MomentumVec targ_p_vec = {p_i[0], p_i[1], p_i[2], -p_i[3]};

            get_tau_form_and_moversigma(string_evolution_mode, y_rem,
                                        tau_form, m_over_sigma, y_loss);
            bool has_baryon_left = false;
            bool has_baryon_right = iproj->is_remnant_carry_baryon_number();
            QCDString qcd_string(x_i, tau_form, iproj, iproj,
                                 p_i, targ_p_vec, m_over_sigma,
                                 has_baryon_right, has_baryon_left);
            qcd_string.set_has_remnant_right(true);
            qcd_string.evolve_QCD_string();
            qcd_string.set_final_baryon_rapidities(0., y_rem - y_loss);
            remnant_string_list_.push_back(qcd_string);
        }
    }
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &itarg: (*targ_nucleon_list)) {
        if (itarg->is_wounded()) {
            auto x_i = itarg->get_remnant_x_frez();

            auto p_i = itarg->get_remnant_p();
            if (p_i[0] <= 0.) continue;
            auto y_rem = -ybeam;
            if (std::abs(p_i[3]) < p_i[0]) {
                // a time-like beam remnant
                y_rem = 0.5*log((p_i[0] + p_i[3])/(p_i[0] - p_i[3]));
            }
            auto cosh_y_rem = cosh(y_rem);
            auto m_rem = p_i[0]/cosh_y_rem;
            p_i[1] = 0.;
            p_i[2] = 0.;
            p_i[3] = m_rem*sinh(y_rem);
            MomentumVec proj_p_vec = {p_i[0], p_i[1], p_i[2], -p_i[3]};

            get_tau_form_and_moversigma(string_evolution_mode, std::abs(y_rem),
                                        tau_form, m_over_sigma, y_loss);
            bool has_baryon_left = itarg->is_remnant_carry_baryon_number();
            bool has_baryon_right = false;
            QCDString qcd_string(x_i, tau_form, itarg, itarg,
                                 proj_p_vec, p_i, m_over_sigma,
                                 has_baryon_right, has_baryon_left);
            qcd_string.set_has_remnant_left(true);
            qcd_string.evolve_QCD_string();
            qcd_string.set_final_baryon_rapidities(y_rem + y_loss, 0.);
            remnant_string_list_.push_back(qcd_string);
        }
    }
}


void Glauber::update_collision_schedule(
                                shared_ptr<CollisionEvent> event_happened) {
    auto proj = event_happened->get_proj_nucleon_ptr().lock();
    proj->increment_collided_times();
    for (auto &it: (*proj->get_collide_nucleon_list()))
        create_a_collision_event(proj, it.lock());
    auto targ = event_happened->get_targ_nucleon_ptr().lock();
    targ->increment_collided_times();
    for (auto &it: (*targ->get_collide_nucleon_list()))
        create_a_collision_event(it.lock(), targ);
}

void Glauber::output_QCD_strings(std::string filename, const real Npart,
                                 const real Ncoll, const real Nstrings,
                                 const real b) {
    std::ofstream output(filename.c_str());
    real total_energy = Npart*parameter_list.get_roots()/2.;
    real net_Pz = ((projectile->get_number_of_wounded_nucleons()
                    - target->get_number_of_wounded_nucleons())
                   *parameter_list.get_roots()/2.);
    output << "# b = " << b << " fm " << "Npart = " << Npart
           << " Ncoll = " << Ncoll << " Nstrings = " << Nstrings
           << " total_energy = " << total_energy << " GeV, "
           << "net_Pz = " << net_Pz << " GeV" << endl;

    output << "# mass[GeV]  m_over_sigma[fm]  tau_form[fm]  tau_0[fm]  eta_s_0  "
           << "x_perp[fm]  y_perp[fm]  x_l[fm]  y_l[fm]  x_r[fm]  y_r[fm]  "
           << "eta_s_left  eta_s_right  y_l  y_r  remnant_l  remnant_r "
           << "y_l_i  y_r_i "
           << "eta_s_baryon_left  eta_s_baryon_right  y_l_baryon  y_r_baryon  "
           << "baryon_fraction_l  baryon_fraction_r"
           << endl;

    // output strings
    for (auto &it: QCD_string_list) {
        auto x_prod = it.get_x_production();
        auto x_left = it.get_targ().lock()->get_x();
        auto x_right = it.get_proj().lock()->get_x();
        if (sample_valence_quark) {
            auto xq_left = it.get_targ_q().lock()->get_x();
            x_left[1]  += xq_left[1];
            x_left[2]  += xq_left[2];

            auto xq_right = it.get_proj_q().lock()->get_x();
            x_right[1] += xq_right[1];
            x_right[2] += xq_right[2];
        }
        auto tau_0  = sqrt(x_prod[0]*x_prod[0] - x_prod[3]*x_prod[3]);
        auto etas_0 = 0.5*log((x_prod[0] + x_prod[3])/(x_prod[0] - x_prod[3]));

        real remnant_left  = 0.;
        if (it.get_has_remnant_left()) {
            remnant_left = 1.0;
        }

        real remnant_right = 0.;
        if (it.get_has_remnant_right()) {
            remnant_right = 1.0;
        }

        real baryon_fraction_left  = 0.;
        if (it.get_has_baryon_left()) {
            baryon_fraction_left = 1.0;
        }

        real baryon_fraction_right = 0.;
        if (it.get_has_baryon_right()) {
            baryon_fraction_right = 1.0;
        }

        auto mass = it.get_string_mass();
        std::vector<real> output_array = {
            mass, it.get_m_over_sigma(), it.get_tau_form(),
            tau_0, etas_0, x_prod[1], x_prod[2],
            x_left[1], x_left[2], x_right[1], x_right[2],
            it.get_eta_s_left(), it.get_eta_s_right(),
            it.get_y_f_left(), it.get_y_f_right(),
            remnant_left, remnant_right,
            it.get_y_i_left(), it.get_y_i_right(),
            it.get_eta_s_baryon_left(), it.get_eta_s_baryon_right(),
            it.get_y_f_baryon_left(), it.get_y_f_baryon_right(),
            baryon_fraction_left, baryon_fraction_right,
        };

        output << std::scientific << std::setprecision(8);
        for (auto &ival : output_array) {
            output << std::setw(15) << ival << "  ";
        }
        output << endl;
    }

    // output the beam remnant strings
    if (sample_valence_quark) {
        for (auto &it: remnant_string_list_) {
            auto x_prod = it.get_x_production();
            auto x_left = it.get_targ().lock()->get_x();
            auto x_right = it.get_proj().lock()->get_x();
            auto tau_0  = sqrt(x_prod[0]*x_prod[0] - x_prod[3]*x_prod[3]);
            auto etas_0 = 0.5*log((x_prod[0] + x_prod[3])/(x_prod[0] - x_prod[3]));

            real remnant_left  = 0.;
            if (it.get_has_remnant_left()) {
                remnant_left = 1.0;
            }

            real remnant_right = 0.;
            if (it.get_has_remnant_right()) {
                remnant_right = 1.0;
            }

            real baryon_fraction_left  = 0.;
            if (it.get_has_baryon_left()) {
                baryon_fraction_left = 1.0;
            }

            real baryon_fraction_right = 0.;
            if (it.get_has_baryon_right()) {
                baryon_fraction_right = 1.0;
            }

            auto mass = it.get_string_mass();
            auto eta_s_center = (it.get_eta_s_left() + it.get_eta_s_right())/2.;
            auto eta_s_left = (  remnant_left*it.get_eta_s_left()
                               + (1. - remnant_left)*eta_s_center);
            auto eta_s_right = (  remnant_right*it.get_eta_s_right()
                                + (1. - remnant_right)*eta_s_center);
            std::vector<real> output_array = {
                mass, it.get_m_over_sigma(), it.get_tau_form(),
                tau_0, etas_0, x_prod[1], x_prod[2],
                x_left[1], x_left[2], x_right[1], x_right[2],
                eta_s_left, eta_s_right,
                remnant_left*it.get_y_f_left(),
                remnant_right*it.get_y_f_right(),
                remnant_left, remnant_right,
                remnant_left*it.get_y_i_left(),
                remnant_right*it.get_y_i_right(),
                remnant_left*it.get_eta_s_baryon_left(),
                remnant_right*it.get_eta_s_baryon_right(),
                remnant_left*it.get_y_f_baryon_left(),
                remnant_right*it.get_y_f_baryon_right(),
                baryon_fraction_left, baryon_fraction_right,
            };

            output << std::scientific << std::setprecision(8);
            for (auto &ival : output_array) {
                output << std::setw(15) << ival << "  ";
            }
            output << endl;
        }
    }
    output.close();
}


real Glauber::sample_rapidity_loss_shell(real y_init) const {
    real y_loss = 0.0;
    if (parameter_list.get_rapidity_loss_method() == 1) {
        y_loss = sample_rapidity_loss_from_the_LEXUS_model(y_init);
    } else if (parameter_list.get_rapidity_loss_method() == 2) {
        y_loss = sample_rapidity_loss_from_parametrization(y_init);
    } else if (parameter_list.get_rapidity_loss_method() == 3) {
        y_loss = sample_rapidity_loss_from_parametrization_with_fluct(y_init);
    }
    return(y_loss);
}


real Glauber::sample_rapidity_loss_from_the_LEXUS_model(
                                            const real y_init) const {
    const real shape_coeff = 1.0;
    real sinh_y_lrf = sinh(shape_coeff*y_init);
    real arcsinh_factor = (ran_gen_ptr_->rand_uniform()
                           *(sinh(2.*shape_coeff*y_init) - sinh_y_lrf)
                           + sinh_y_lrf);
    real y_loss = 2.*y_init - asinh(arcsinh_factor)/shape_coeff;
    return(y_loss);
}


real Glauber::sample_rapidity_loss_from_parametrization(
                                                const real y_init) const {
    auto y_loss = (
        yloss_param_slope*pow(pow(y_init, yloss_param_a)*tanh(y_init),
                              yloss_param_b));
    return(y_loss);
}


real Glauber::sample_rapidity_loss_from_parametrization_with_fluct(
                                                const real y_init) const {
    auto y_mean = sample_rapidity_loss_from_parametrization(y_init);
    auto y_rhic = 5.5;
    auto y_lhc = 9.0;
    auto var_rhic = parameter_list.get_yloss_param_fluct_var_RHIC();
    auto var_lhc = parameter_list.get_yloss_param_fluct_var_LHC();
    auto var = var_rhic;
    if (y_init > y_lhc) {
        var = var_lhc;
    } else if (y_init > y_rhic) {
        var = (var_rhic
               + (y_init - y_rhic)/(y_lhc - y_rhic)*(var_lhc - var_rhic));
    }

    auto random_x = ran_gen_ptr_->rand_normal(0., var);

    real logit_rand = 1./(1. + exp(-random_x));
    auto y_loss = y_mean;
    if (std::abs(2.*y_mean - y_init) > 1e-15) {
        real aa = (2.*y_mean - y_init)/(2.*y_mean*y_init*(y_init - y_mean));
        real bb = ((y_init*y_init - 2.*y_mean*y_mean)
                   /(2.*y_init*y_mean*(y_init - y_mean)));
        real cc = - logit_rand;
        y_loss = (-bb + sqrt(bb*bb - 4.*aa*cc))/(2.*aa);
    } else {
        y_loss = logit_rand*y_init;
    }
    return(y_loss);
}


// sample y from exp[(y - (0.5*(yt + yp)))/2]/(4.*Sinh[0.25*yp - 0.25*yt]),
// the new rapidity of the baryon number from the right moving particle
// after the collision in the lab frame
real Glauber::sample_junction_rapidity_right(real y_left, real y_right) const {
    real y = -2.*(-0.25*y_right - 0.25*y_left
                  - 1.*log(2.*(ran_gen_ptr_->rand_uniform()
                               + 0.5*exp(-0.25*y_right + 0.25*y_left)
                                 /sinh(0.25*y_right - 0.25*y_left))
                           *sinh(0.25*y_right - 0.25*y_left)));
    return(y);
}

// sample y from exp[-(y - (0.5*(yt + yp)))/2]/(4.*Sinh[0.25*yp - 0.25*yt]),
// the new rapidity of the baryon number from the left moving particle
// after the collision in the lab frame
real Glauber::sample_junction_rapidity_left(real y_left, real y_right) const {
    real y = 2.*(0.25*y_right + 0.25*y_left
                 - log(2.*(ran_gen_ptr_->rand_uniform()
                           + 0.5*exp(-0.25*y_right + 0.25*y_left)
                             /sinh(0.25*y_right - 0.25*y_left))
                       *sinh( 0.25 * y_right - 0.25 * y_left)));
    return(y);
}


// This function computes the sigeff(s) from the formula
//  sigmaNN_in(s) = int d^2b [1 - exp(-sigeff(s)*Tpp(b))]
//  Reads sigmaNN, returns guassian width
real Glauber::get_sig_eff(const real siginNN) {
    // rms-radius of a gaussian = rms-radius of a disc with radius R,
    // where 2*PI*(2R)^2=sigmaNN
    const real width = sqrt(0.1*siginNN/M_PI)/sqrt(8.);
    nucleon_width_ = width;
    const real sigin = siginNN*0.1;   // sigma_in(s) [mb --> fm^2]

    const int Nint = 100;    // # of integration points
    std::vector<real> b(Nint, 0.);
    std::vector<real> Tnn(Nint, 0.);
    const real Bmax = 10.0*width;
    const real db = Bmax/Nint;
    for (int i = 0; i < Nint; i++) {
        b[i] = (i + 0.5)*db;
        Tnn[i] = exp(-b[i]*b[i]/(4.*width*width))/(M_PI*(4.*width*width));
    }
    const real prefactor = 2.*M_PI*db;

    real sum, dN;

    real sigeff = 10.0;        // starting point of iteration [fm^2]
    real sigeff0;     // holds value from previous iteration step
    int tol = 0;
    do {                                       // iterate ...
        sigeff0 = sigeff;
        sum = 0.0;
        dN  = 0.0;
        for (int ib = 0; ib < Nint; ib++) {  // integral d^2b from 0 to Bmax
            sum += prefactor*b[ib]*(1.0 - exp(-sigeff*Tnn[ib]));
            dN  += prefactor*b[ib]*Tnn[ib]*exp(-sigeff*Tnn[ib]);
        }
        sigeff -= (sum - sigin)/dN;
        tol++;
        //cout << "iter: " << tol << ": sigeff = " << sigeff
        //     << " fm^2, sum = " << sum
        //     << " fm^2, sigin = " << sigin << " fm^2" << endl;
    } while(std::abs(sigeff - sigeff0) > 1e-4 && tol < 1000);
    // until sigeff has converged

    return(sigeff);
}

}
