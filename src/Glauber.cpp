// Copyright @ Chun Shen 2018

#include "Glauber.h"
#include "data_structs.h"
#include "PhysConsts.h"

#include <string>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <memory>

using std::cout;
using std::endl;
using std::shared_ptr;

// put in use of baryon used flags

namespace MCGlb {

Glauber::Glauber(const MCGlb::Parameters &param_in,
                 shared_ptr<RandomUtil::Random> ran_gen) :
    parameter_list(param_in) {
    sample_valence_quark = false;
    if (parameter_list.get_use_quarks() > 0) {
        sample_valence_quark = true;
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
    ran_gen_ptr = ran_gen;
    string_production_mode = parameter_list.get_QCD_string_production_mode();

    yloss_param_slope = parameter_list.get_yloss_param_slope();
    real alpha1 = parameter_list.get_yloss_param_alpha1();
    real alpha2 = parameter_list.get_yloss_param_alpha2();
    real alpha = alpha2/alpha1;
    yloss_param_a = alpha/(1. - alpha);
    yloss_param_b = alpha2/yloss_param_a;
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
            (b_max*b_max - b_min*b_min)*ran_gen_ptr.lock()->rand_uniform());
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
    auto d2 = (compute_NN_inelastic_cross_section(parameter_list.get_roots())
               /(M_PI*10.));  // in fm^2 
    auto proj_nucleon_list = projectile->get_nucleon_list();
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &iproj: (*proj_nucleon_list)) {
        auto proj_x = iproj->get_x();
        for (auto &itarg: (*targ_nucleon_list)) {
            auto targ_x = itarg->get_x();
            auto dij = (  (targ_x[1] - proj_x[1])*(targ_x[1] - proj_x[1])
                        + (targ_x[2] - proj_x[2])*(targ_x[2] - proj_x[2]));
            if (hit(dij, d2)) {
                create_a_collision_event(iproj, itarg);
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

bool Glauber::hit(real d2, real d2_in) const {
    real G = 0.92;  // from Glassando
    return(ran_gen_ptr.lock()->rand_uniform() < G*exp(-G*d2/d2_in));
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
            event_ptr->set_produced_a_string(true);
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


bool Glauber::decide_produce_string(shared_ptr<CollisionEvent> event_ptr) const {
    bool flag_form_a_string = false;
    auto proj = event_ptr->get_proj_nucleon_ptr().lock();
    auto targ = event_ptr->get_targ_nucleon_ptr().lock();
    int minimum_allowed_connections = 1;
    if (sample_valence_quark)
        minimum_allowed_connections = PhysConsts::NumValenceQuark;
    if (   proj->get_number_of_connections() < minimum_allowed_connections
        || targ->get_number_of_connections() < minimum_allowed_connections) {
        // bps: string gets flag for whether it has left or right
        // or both baryon numbers.
        flag_form_a_string = true;
    } else {
        int n_connects = (  proj->get_number_of_connections()
                          + targ->get_number_of_connections());
        real shadowing_factor = parameter_list.get_shadowing_factor();
        real production_prob = (
            (1. - shadowing_factor)
            *exp(-shadowing_factor*(n_connects
                                    - 2.*minimum_allowed_connections)));
        if (ran_gen_ptr.lock()->rand_uniform() < production_prob) {
            flag_form_a_string = true;
        }
    }
    return(flag_form_a_string);
}


int Glauber::decide_QCD_strings_production() {
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
            auto flag_form_a_string = decide_produce_string(ievent);
            auto proj = ievent->get_proj_nucleon_ptr();
            auto targ = ievent->get_targ_nucleon_ptr();
            if (flag_form_a_string) {
                number_of_strings++;
                proj.lock()->add_connected_nucleon(targ);
                targ.lock()->add_connected_nucleon(proj);
                ievent->set_produced_a_string(true);
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
    auto y_i = q_i->get_rapidity();
    auto y_f = y_i + y_shift;
    q_i->set_rapidity(y_f);
    pvec[0] = PhysConsts::MProton*cosh(y_f);
    pvec[3] = PhysConsts::MProton*sinh(y_f);
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

// add flag if baryon used to proj and target
// - check it - only put in baryon if not used yet, then set it to used
int Glauber::perform_string_production() {
    if (sample_valence_quark) {
        projectile->sample_valence_quarks_inside_nucleons(
                                    parameter_list.get_roots(), 1);
        target->sample_valence_quarks_inside_nucleons(
                                    parameter_list.get_roots(), -1);
    }
    QCD_string_list.clear();
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
        if (!first_event->is_produced_a_string()) {
            collision_schedule.erase((*collision_schedule.begin()));
            continue;
        }
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
                MomentumVec p_q = {
                    PhysConsts::MQuarkValence*cosh(proj_q->get_rapidity()),
                    0.0,
                    0.0,
                    PhysConsts::MQuarkValence*sinh(proj_q->get_rapidity())};
                proj->substract_momentum_from_remnant(p_q);
            }
            targ_q = targ->get_a_valence_quark();
            if (targ_q->get_number_of_connections() == 1) {
                // first time pick-up the valence quark
                // we need to substract the valence quark energy-momentum
                // from the nucleon remnant energy-momentum vector
                MomentumVec p_q = {
                    PhysConsts::MQuarkValence*cosh(targ_q->get_rapidity()),
                    0.0,
                    0.0,
                    PhysConsts::MQuarkValence*sinh(targ_q->get_rapidity())};
                targ->substract_momentum_from_remnant(p_q);
            }
            y_in_lrf = std::abs(  proj_q->get_rapidity()
                                - targ_q->get_rapidity())/2.;
        }
        real tau_form = 0.5;      // [fm]
        real m_over_sigma = 1.0;  // [fm]
        real y_loss = 0.;
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
            tau_form = 0.5 + 1.*ran_gen_ptr.lock()->rand_uniform();
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
        // set variables in case of no baryon junction transport
        bool has_baryon_left = false;
        bool has_baryon_right = false;
        if (!sample_valence_quark) {
            QCDString qcd_string(x_coll, tau_form, proj, targ, m_over_sigma,
                                 has_baryon_right, has_baryon_left);
            QCD_string_list.push_back(qcd_string);
        } else {
            QCDString qcd_string(x_coll, tau_form, proj, targ,
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
            update_momentum(targ, -1e-3*targ->get_rapidity());
            update_momentum_quark(proj_q, -y_shift);
            update_momentum_quark(targ_q, y_shift);
        }
        update_collision_schedule(first_event);
        collision_schedule.erase((*collision_schedule.begin()));
    }

    // randomize the QCD_string_list and assign the baryon charge to
    // the strings
    std::vector<unsigned int> random_idx;
    for (unsigned int idx = 0; idx < QCD_string_list.size(); idx++)
        random_idx.push_back(idx);
    std::random_shuffle(random_idx.begin(), random_idx.end());
    for (auto &idx: random_idx) {
        auto proj = QCD_string_list[idx].get_proj();
        if (!proj.lock()->baryon_was_used()) {
            proj.lock()->set_baryon_used(true);
            QCD_string_list[idx].set_has_baryon_right(true);
        }
        auto targ = QCD_string_list[idx].get_targ();
        if (!targ.lock()->baryon_was_used()) {
            targ.lock()->set_baryon_used(true);
            QCD_string_list[idx].set_has_baryon_left(true);
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
                if (ran_gen_ptr.lock()->rand_uniform() < lambdaB) {
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
                if (ran_gen_ptr.lock()->rand_uniform() < lambdaB) {
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
                auto x_frez = it->get_x_production();
                proj_n.lock()->set_remnant_x_frez(x_frez);
            }
            auto targ_n = it->get_targ();
            if (!targ_n.lock()->is_remnant_set()) {
                targ_n.lock()->set_remnant(true);
                auto x_frez = it->get_x_production();
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
    return(number_of_collided_events);
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
    output << "# b = " << b << " fm " << "Npart = " << Npart
           << " Ncoll = " << Ncoll << " Nstrings = " << Nstrings
           << " total_energy = " << total_energy << " GeV" << endl;

    output << "# norm  m_over_sigma[fm]  tau_form[fm]  tau_0[fm]  eta_s_0  "
           << "x_perp[fm]  y_perp[fm]  "
           << "eta_s_left  eta_s_right  y_l  y_r  remnant_l  remnant_r "
           << "y_l_i  y_r_i "
           << "eta_s_baryon_left  eta_s_baryon_right  y_l_baryon  y_r_baryon  "
           << "baryon_fraction_l  baryon_fraction_r"
           << endl;

    real remnant_left  = 0.;
    real remnant_right = 0.;
    real baryon_fraction_left  = 0.;
    real baryon_fraction_right = 0.;

    for (auto &it: QCD_string_list) {
        auto x_prod = it.get_x_production();
        auto tau_0  = sqrt(x_prod[0]*x_prod[0] - x_prod[3]*x_prod[3]);
        auto etas_0 = 0.5*log((x_prod[0] + x_prod[3])/(x_prod[0] - x_prod[3]));

        if (it.get_has_remnant_left()) {
            remnant_left = 1.0;
        } else {
            remnant_left = 0.0;
        }

        if (it.get_has_remnant_right()) {
            remnant_right = 1.0;
        } else {
            remnant_right = 0.0;
        }

        if (it.get_has_baryon_left()) {
            baryon_fraction_left = 1.0;
        } else {
            baryon_fraction_left = 0.0;
        }

        if (it.get_has_baryon_right()) {
            baryon_fraction_right = 1.0;
        } else {
            baryon_fraction_right = 0.0;
        }

        std::vector<real> output_array = {
            1.0, it.get_m_over_sigma(), it.get_tau_form(),
            tau_0, etas_0, x_prod[1], x_prod[2],
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
    output.close();
}


void Glauber::output_remnants(std::string filename) const {
    std::ofstream output(filename.c_str());
    output << "# t[fm]  x[fm]  y[fm]  z[fm]  E[GeV]  px[GeV]  py[GeV]  pz[GeV]"
           << endl;
    auto proj_nucleon_list = projectile->get_nucleon_list();
    auto targ_nucleon_list = target->get_nucleon_list();
    output << std::scientific << std::setprecision(8);
    for (auto &iproj: (*proj_nucleon_list)) {
        if (iproj->is_wounded()) {
            auto x_i = iproj->get_remnant_x_frez();
            for (int i = 0; i < 4; i++)
                output << std::setw(15) << x_i[i] << "  ";
            auto p_i = iproj->get_remnant_p();
            for (int i = 0; i < 4; i++)
                output << std::setw(15) << p_i[i] << "  ";
            output << endl;
        }
    }
    for (auto &itarg: (*targ_nucleon_list)) {
        if (itarg->is_wounded()) {
            auto x_i = itarg->get_remnant_x_frez();
            for (int i = 0; i < 4; i++)
                output << std::setw(15) << x_i[i] << "  ";
            auto p_i = itarg->get_remnant_p();
            for (int i = 0; i < 4; i++)
                output << std::setw(15) << p_i[i] << "  ";
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
    }
    return(y_loss);
}

real Glauber::sample_rapidity_loss_from_the_LEXUS_model(real y_init) const {
    const real shape_coeff = 1.0;
    real sinh_y_lrf = sinh(shape_coeff*y_init);
    real arcsinh_factor = (ran_gen_ptr.lock()->rand_uniform()
                           *(sinh(2.*shape_coeff*y_init) - sinh_y_lrf)
                           + sinh_y_lrf);
    real y_loss = 2.*y_init - asinh(arcsinh_factor)/shape_coeff;
    return(y_loss);
}

real Glauber::sample_rapidity_loss_from_parametrization(real y_init) const {
    auto y_loss = (
        yloss_param_slope*pow(pow(y_init, yloss_param_a)*tanh(y_init),
                              yloss_param_b));
    return(y_loss);
}

// sample y from exp[(y - (0.5 (yt + yp)))/2]/(4.` Sinh[0.25` yp - 0.25` yt]),
// the new rapidity of the baryon number from the right moving particle
// after the collision in the lab frame
real Glauber::sample_junction_rapidity_right(real y_left, real y_right) const {
    real y = -2.*(-0.25*y_right - 0.25*y_left
                  - 1.*log(2.*(ran_gen_ptr.lock()->rand_uniform()
                               + 0.5*exp(-0.25*y_right + 0.25*y_left)
                                 /sinh(0.25*y_right - 0.25*y_left))
                           *sinh(0.25*y_right - 0.25*y_left)));
    return(y);
}

// sample y from exp[-(y - (0.5 (yt + yp)))/2]/(4.` Sinh[0.25` yp - 0.25` yt]),
// the new rapidity of the baryon number from the left moving particle
// after the collision in the lab frame
real Glauber::sample_junction_rapidity_left(real y_left, real y_right) const {
    real y = 2.*(0.25*y_right + 0.25*y_left
                 - log(2.*(ran_gen_ptr.lock()->rand_uniform()
                           + 0.5*exp(-0.25*y_right + 0.25*y_left)
                             /sinh(0.25*y_right - 0.25*y_left))
                       *sinh( 0.25 * y_right - 0.25 * y_left)));
    return(y);
}

}
