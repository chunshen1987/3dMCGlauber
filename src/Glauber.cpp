// Copyright @ Chun Shen 2018

#include "Glauber.h"

#include <math.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "PhysConsts.h"
#include "data_structs.h"

using std::shared_ptr;

namespace MCGlb {

Glauber::Glauber(
    const MCGlb::Parameters &param_in, shared_ptr<RandomUtil::Random> ran_gen)
    : parameter_list(param_in) {
    sample_valence_quark = false;
    fluct_Nstrings_per_NN_collision_ = false;
    remnant_energy_loss_fraction_ =
        (parameter_list.get_remnant_energy_loss_fraction());
    if (parameter_list.get_use_quarks() > 0) {
        sample_valence_quark = true;
        fluct_Nstrings_per_NN_collision_ =
            (parameter_list.get_fluct_Nstrings_per_NN_collision());
        if (!parameter_list.get_cached_tabels()) {
            system_status_ =
                std::system("rm -fr tables/proton_valence_quark_samples*");
            system_status_ =
                std::system("rm -fr tables/neutron_valence_quark_samples*");
        }
    }

    real d_min = parameter_list.get_d_min();

    int N_sea_partons = parameter_list.get_N_sea_partons();

    bool deformed = true;
    bool nucleonConfFromFile = parameter_list.nucleon_configuration_from_file();
    projectile = std::unique_ptr<Nucleus>(new Nucleus(
        parameter_list.get_projectle_nucleus_name(), ran_gen,
        sample_valence_quark, parameter_list.get_BG(), d_min, deformed,
        nucleonConfFromFile, N_sea_partons));
    int resetProjWS =
        static_cast<int>(parameter_list.getParam("resetProjWS", 0.0));
    if (resetProjWS != 0) {
        auto defaultWSParam = projectile->get_woods_saxon_parameters();
        real WS_rho = parameter_list.getParam("ProjWS_rho0", defaultWSParam[0]);
        real WS_w = parameter_list.getParam("ProjWS_w", defaultWSParam[1]);
        real WS_R = parameter_list.getParam("ProjWS_R", defaultWSParam[2]);
        real WS_a = parameter_list.getParam("ProjWS_a", defaultWSParam[3]);
        real WS_beta2 =
            parameter_list.getParam("ProjWS_beta2", defaultWSParam[4]);
        real WS_beta3 =
            parameter_list.getParam("ProjWS_beta3", defaultWSParam[5]);
        real WS_beta4 =
            parameter_list.getParam("ProjWS_beta4", defaultWSParam[6]);
        real WS_gamma =
            parameter_list.getParam("ProjWS_gamma", defaultWSParam[7]);
        projectile->setWoodsSaxonParameters(
            WS_rho, WS_w, WS_R, WS_a, WS_beta2, WS_beta3, WS_beta4, WS_gamma);
    }

    target = std::unique_ptr<Nucleus>(new Nucleus(
        parameter_list.get_target_nucleus_name(), ran_gen, sample_valence_quark,
        parameter_list.get_BG(), d_min, deformed, nucleonConfFromFile,
        N_sea_partons));
    int resetTargWS =
        static_cast<int>(parameter_list.getParam("resetTargWS", 0.0));
    if (resetTargWS != 0) {
        auto defaultWSParam = target->get_woods_saxon_parameters();
        real WS_rho, WS_w, WS_R, WS_a, WS_beta2, WS_beta3, WS_beta4, WS_gamma;
        if (parameter_list.get_target_nucleus_name()
            == parameter_list.get_projectle_nucleus_name()) {
            // if it is a symmetric collision, then use the projectile WS param
            WS_rho = parameter_list.getParam("ProjWS_rho0", defaultWSParam[0]);
            WS_w = parameter_list.getParam("ProjWS_w", defaultWSParam[1]);
            WS_R = parameter_list.getParam("ProjWS_R", defaultWSParam[2]);
            WS_a = parameter_list.getParam("ProjWS_a", defaultWSParam[3]);
            WS_beta2 =
                parameter_list.getParam("ProjWS_beta2", defaultWSParam[4]);
            WS_beta3 =
                parameter_list.getParam("ProjWS_beta3", defaultWSParam[5]);
            WS_beta4 =
                parameter_list.getParam("ProjWS_beta4", defaultWSParam[6]);
            WS_gamma =
                parameter_list.getParam("ProjWS_gamma", defaultWSParam[7]);
        } else {
            // for asymmetric collisions, set with target WS params
            WS_rho = parameter_list.getParam("TargWS_rho0", defaultWSParam[0]);
            WS_w = parameter_list.getParam("TargWS_w", defaultWSParam[1]);
            WS_R = parameter_list.getParam("TargWS_R", defaultWSParam[2]);
            WS_a = parameter_list.getParam("TargWS_a", defaultWSParam[3]);
            WS_beta2 =
                parameter_list.getParam("TargWS_beta2", defaultWSParam[4]);
            WS_beta3 =
                parameter_list.getParam("TargWS_beta3", defaultWSParam[5]);
            WS_beta4 =
                parameter_list.getParam("TargWS_beta4", defaultWSParam[6]);
            WS_gamma =
                parameter_list.getParam("TargWS_gamma", defaultWSParam[7]);
        }
        target->setWoodsSaxonParameters(
            WS_rho, WS_w, WS_R, WS_a, WS_beta2, WS_beta3, WS_beta4, WS_gamma);
    }
    if (nucleonConfFromFile) {
        projectile->setLightNucleusOption(
            parameter_list.getLightNucleusOption());
        target->setLightNucleusOption(parameter_list.getLightNucleusOption());
    }

    if (sample_valence_quark) {
        projectile->set_valence_quark_Q2(parameter_list.get_quarks_Q2());
        target->set_valence_quark_Q2(parameter_list.get_quarks_Q2());
    }
    ran_gen_ptr_ = ran_gen;

    yloss_param_slope = parameter_list.get_yloss_param_slope();
    real alpha1 = parameter_list.get_yloss_param_alpha1();
    real alpha2 = parameter_list.get_yloss_param_alpha2();
    real alpha = alpha2 / alpha1;
    yloss_param_a = alpha / (1. - alpha);
    yloss_param_b = alpha2 / yloss_param_a;
    collision_energy_ = parameter_list.get_roots();
    real ybeam_AA = acosh(collision_energy_ / (2. * PhysConsts::MProton));
    real rootgammaN_low = 2.0;
    real rootgammaN_up = parameter_list.get_roots() / 2.0;

    if (parameter_list.use_roots_distribution()) {
        if (parameter_list.use_roots_cut()) {
            rootgammaN_low = parameter_list.get_UPC_root_low_cut();
            rootgammaN_up = parameter_list.get_UPC_root_up_cut();
        }
        collision_energy_ = get_roots_from_distribution(
            parameter_list.get_roots(), rootgammaN_low, rootgammaN_up);
    }
    ybeam = acosh(collision_energy_ / (2. * PhysConsts::MProton));
    std::ofstream output_rapidity_shift("rapidity_shift");
    output_rapidity_shift << parameter_list.use_roots_distribution() << "  "
                          << collision_energy_ << "  " << ybeam - ybeam_AA
                          << std::endl;
    output_rapidity_shift.close();
    real siginNN = compute_NN_inelastic_cross_section(collision_energy_);
    sigma_eff_ = get_sig_eff(siginNN);
    std::cout << "sqrt{s} = " << collision_energy_ << " GeV, "
              << "siginNN = " << siginNN << " mb" << std::endl;
}

real Glauber::get_roots_from_distribution(
    real roots, real rootgammaN_low_cut, real rootgammaN_up_cut) {
    auto defaultWSParam = target->get_woods_saxon_parameters();
    real radius = defaultWSParam[2];
    int Z_in = target->get_nucleus_Z();

    double gammaL = roots / (2. * PhysConsts::MProton);
    double probability_sum = 0.0;
    real rootgammaN_sampled = roots;
    std::vector<double> dNdrootgammaN;
    for (double rootgammaN = rootgammaN_low_cut; rootgammaN < rootgammaN_up_cut;
         rootgammaN++) {
        double k = rootgammaN * rootgammaN / 2. / roots;
        double omegaAA = 2. * k * radius / gammaL / PhysConsts::HBARC;
        // The dN/dk get from the Eq. (6) in arXiv: 0706.3356
        double temp = rootgammaN / roots * 2. * Z_in * Z_in / M_PI / k
                      * (omegaAA * std::cyl_bessel_k(0, omegaAA)
                             * std::cyl_bessel_k(1, omegaAA)
                         - omegaAA * omegaAA / 2.
                               * (std::cyl_bessel_k(1, omegaAA)
                                      * std::cyl_bessel_k(1, omegaAA)
                                  - std::cyl_bessel_k(0, omegaAA)
                                        * std::cyl_bessel_k(0, omegaAA)));
        dNdrootgammaN.push_back(temp);
        probability_sum = probability_sum + temp;
    }
    real MCprobability = ran_gen_ptr_->rand_uniform();
    int index = 0;
    for (double rootgammaN = rootgammaN_low_cut; rootgammaN < rootgammaN_up_cut;
         rootgammaN++) {
        MCprobability = MCprobability - dNdrootgammaN[index] / probability_sum;
        if (MCprobability <= 0.) {
            rootgammaN_sampled = rootgammaN + 1.;
            break;
        }
        index++;
    }
    return (rootgammaN_sampled);
}

void Glauber::make_nuclei() {
    projectile->generate_nucleus_3d_configuration();
    target->generate_nucleus_3d_configuration();
    // projectile->output_nucleon_positions("projectile_before_accelaration.dat");
    // target->output_nucleon_positions("target_before_accelaration.dat");
    projectile->accelerate_nucleus(collision_energy_, 1);
    target->accelerate_nucleus(collision_energy_, -1);

    // sample impact parameters
    auto b_max = parameter_list.get_b_max();
    auto b_min = parameter_list.get_b_min();
    impact_b = sqrt(
        b_min * b_min
        + (b_max * b_max - b_min * b_min) * ran_gen_ptr_->rand_uniform());
    SpatialVec proj_shift = {
        0., impact_b / 2., 0., -projectile->get_z_max() - 1e-15};
    projectile->shift_nucleus(proj_shift);
    SpatialVec targ_shift = {
        0., -impact_b / 2., 0., -target->get_z_min() + 1e-15};
    target->shift_nucleus(targ_shift);
    // projectile->output_nucleon_positions("projectile.dat");
    // target->output_nucleon_positions("target.dat");
}

int Glauber::make_collision_schedule() {
    collision_schedule.clear();
    auto proj_nucleon_list = projectile->get_nucleon_list();
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &iproj : (*proj_nucleon_list)) {
        auto proj_x = iproj->get_x();
        for (auto &itarg : (*targ_nucleon_list)) {
            auto targ_x = itarg->get_x();
            auto dij =
                ((targ_x[1] - proj_x[1]) * (targ_x[1] - proj_x[1])
                 + (targ_x[2] - proj_x[2]) * (targ_x[2] - proj_x[2]));
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

    // copy the collision_schedule information to outside
    collision_schedule_list_.clear();
    int pos = 0;
    for (auto &it : collision_schedule) {
        collision_schedule_list_.push_back(
            *it);  // collision list is time ordered
        // auto xvec = collision_schedule_list_[pos].get_collision_position();
        // std::cout << xvec[0]<<" "<< xvec[1]<<" "<< xvec[2]<<" "<< xvec[3]<<"
        // " << std::endl;
        pos++;
    }
    get_collision_information();
    return (collision_schedule.size());
}

int Glauber::get_Npart() const {
    int Npart =
        (projectile->get_number_of_wounded_nucleons()
         + target->get_number_of_wounded_nucleons());
    return (Npart);
}

bool Glauber::hit(real d2) const {
    // real G = 0.92;  // from Glassando
    // return(ran_gen_ptr_->rand_uniform() < G*exp(-G*d2/d2_in));
    const real T_nn =
        (exp(-d2 / (4. * nucleon_width_ * nucleon_width_))
         / (4. * M_PI * nucleon_width_ * nucleon_width_));
    const real hit_treshold = 1. - exp(-sigma_eff_ * T_nn);
    return (ran_gen_ptr_->rand_uniform() < hit_treshold);
}

void Glauber::create_a_collision_event(
    shared_ptr<Nucleon> proj, shared_ptr<Nucleon> targ) {
    real t_coll, z_coll;
    auto x1 = proj->get_x();
    auto p1 = proj->get_p();
    auto x2 = targ->get_x();
    auto p2 = targ->get_p();
    auto v1 = p1[3] / p1[0];
    auto v2 = p2[3] / p2[0];
    auto collided =
        get_collision_point(x1[0], x1[3], v1, x2[3], v2, t_coll, z_coll);
    if (collided) {
        SpatialVec x_coll = {
            t_coll, (x1[1] + x2[1]) / 2., (x1[2] + x2[2]) / 2., z_coll};
        shared_ptr<CollisionEvent> event_ptr(
            new CollisionEvent(x_coll, proj, targ));
        if (proj->is_connected_with(targ)) {
            auto form_N_strings = proj->get_number_of_connections(targ);
            event_ptr->set_produced_n_strings(form_N_strings);
        }
        collision_schedule.insert(event_ptr);
    }
}

bool Glauber::get_collision_point(
    real t, real z1, real v1, real z2, real v2, real &t_coll,
    real &z_coll) const {
    bool collided = false;
    real delta_t = (z2 - z1) / (v1 - v2);
    if (delta_t > 0.) {
        collided = true;
        t_coll = t + delta_t;
        z_coll = z1 + v1 * delta_t;
    } else {
        t_coll = -1.;
        z_coll = -1.;
    }
    return (collided);
}

real Glauber::compute_NN_inelastic_cross_section(real ecm) const {
    real s = ecm * ecm;
    real sigma_NN_total = 44.4 - 2.9 * log(s) + 0.33 * log(s) * log(s);
    real sigma_NN_inel =
        (sigma_NN_total - (11.4 - 1.52 * log(s) + 0.13 * log(s) * log(s)));
    return (sigma_NN_inel);
}

int Glauber::decide_produce_string_num(
    shared_ptr<CollisionEvent> event_ptr) const {
    int form_n_string = 0;
    auto proj = event_ptr->get_proj_nucleon_ptr().lock();
    auto targ = event_ptr->get_targ_nucleon_ptr().lock();
    int minimum_allowed_connections = 1;
    if (sample_valence_quark && fluct_Nstrings_per_NN_collision_) {
        // assume P(N) = 1/e^N distribution
        int N = std::min(
            proj->get_number_of_quarks(), targ->get_number_of_quarks());
        do {
            auto rand = ran_gen_ptr_->rand_uniform();
            minimum_allowed_connections = static_cast<int>(-log(1 - rand)) + 1;
        } while (minimum_allowed_connections > N);
    }

    if (proj->get_number_of_connections() < minimum_allowed_connections
        && targ->get_number_of_connections() < minimum_allowed_connections) {
        form_n_string = std::min(
            minimum_allowed_connections - proj->get_number_of_connections(),
            minimum_allowed_connections - targ->get_number_of_connections());
    } else if (
        proj->get_number_of_connections() < minimum_allowed_connections
        || targ->get_number_of_connections() < minimum_allowed_connections) {
        form_n_string = 1;
    } else {
        int n_connects =
            (proj->get_number_of_connections()
             + targ->get_number_of_connections());
        real shadowing_factor = parameter_list.get_shadowing_factor();
        real production_prob =
            ((1. - shadowing_factor)
             * exp(
                 -shadowing_factor
                 * (n_connects - 2. * minimum_allowed_connections)));
        if (ran_gen_ptr_->rand_uniform() < production_prob) {
            form_n_string = 1;
        }
    }
    return (form_n_string);
}

int Glauber::decide_QCD_strings_production() {
    if (sample_valence_quark) {
        projectile->sample_valence_quarks_inside_nucleons(collision_energy_, 1);
        target->sample_valence_quarks_inside_nucleons(collision_energy_, -1);
        projectile->add_soft_parton_ball(collision_energy_, 1);
        target->add_soft_parton_ball(collision_energy_, -1);
    }
    std::vector<shared_ptr<CollisionEvent>> collision_list;
    for (auto &it : collision_schedule)
        collision_list.push_back(it);  // collision list is time ordered

    const auto QCD_string_production_mode =
        parameter_list.get_QCD_string_production_mode();
    if (QCD_string_production_mode == 1) {
        // randomly ordered strings
        std::shuffle(
            collision_list.begin(), collision_list.end(),
            *ran_gen_ptr_->getRanGenerator());
    } else if (QCD_string_production_mode == 2) {
        // anti-time ordered strings
        std::reverse(collision_list.begin(), collision_list.end());
    }

    int number_of_strings = 0;
    bool finished = false;
    while (!finished) {
        finished = true;
        for (auto &ievent : collision_list) {
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
                if (proj.lock()->get_number_of_connections() == 0
                    || targ.lock()->get_number_of_connections() == 0)
                    finished = false;
            }
        }
    }
    return (number_of_strings);
}

//! This function propagate individual nucleon inside the nucleus by dt
void Glauber::propagate_nuclei(real dt) {
    for (auto &n_i : (*projectile->get_nucleon_list()))
        propagate_nucleon(n_i, dt);
    for (auto &n_i : (*target->get_nucleon_list())) propagate_nucleon(n_i, dt);
}

void Glauber::propagate_nucleon(shared_ptr<Nucleon> n_i, real dt) {
    auto x_vec = n_i->get_x();
    auto p_vec = n_i->get_p();
    for (int i = 0; i < 4; i++) {
        real v_i = p_vec[i] / p_vec[0];
        x_vec[i] += dt * v_i;
    }
    n_i->set_x(x_vec);
}

void Glauber::update_momentum_quark(shared_ptr<Quark> q_i, real y_shift) {
    auto pvec = q_i->get_p();
    auto mass = q_i->get_mass();
    auto y_i = q_i->get_rapidity();
    auto y_f = y_i + y_shift;
    q_i->set_rapidity(y_f);
    pvec[0] = mass * cosh(y_f);
    pvec[3] = mass * sinh(y_f);
    q_i->set_p(pvec);
}

void Glauber::update_momentum(shared_ptr<Nucleon> n_i, real y_shift) {
    auto pvec = n_i->get_p();
    auto y_i = n_i->get_rapidity();
    auto y_f = y_i + y_shift;
    pvec[0] = PhysConsts::MProton * cosh(y_f);
    pvec[3] = PhysConsts::MProton * sinh(y_f);
    n_i->set_p(pvec);
}

real Glauber::get_tau_form(const int string_evolution_mode) const {
    auto tau_form_mean = parameter_list.get_tau_form_mean();
    real tau_form = tau_form_mean;  // [fm]
    if (string_evolution_mode == 2) {
        // tau_form fluctuates with a gamma distribution
        tau_form = tau_form_mean / ran_gen_ptr_->rand_gamma_dis();
        tau_form = std::min(std::max(0.1, tau_form), 6 * tau_form_mean);
    }
    return (tau_form);
}

void Glauber::get_tau_form_and_moversigma(
    const int string_evolution_mode, const real y_in_lrf, real &tau_form,
    real &m_over_sigma, real &y_loss) {
    tau_form = get_tau_form(string_evolution_mode);  // [fm]
    m_over_sigma = 1.0;                              // [fm]
    y_loss = 0.;
    real eps = 1e-16;
    if (string_evolution_mode == 1) {
        // fixed rapidity loss
        y_loss = (acosh(
            tau_form * tau_form / (2. * m_over_sigma * m_over_sigma) + 1.));
        // maximum 90% of y_in_lrf
        y_loss = std::min(y_loss, y_in_lrf * 0.9);
        m_over_sigma = tau_form / std::max(eps, sqrt(2. * (cosh(y_loss) - 1.)));
    } else if (string_evolution_mode == 2) {
        // both tau_form and sigma fluctuate
        y_loss = sample_rapidity_loss_shell(y_in_lrf);
        m_over_sigma = tau_form / std::max(eps, sqrt(2. * (cosh(y_loss) - 1.)));
    } else if (string_evolution_mode == 3) {
        // only tau_form fluctuates
        y_loss = sample_rapidity_loss_shell(y_in_lrf);
        tau_form = m_over_sigma * sqrt(2. * (cosh(y_loss) - 1.));
    } else if (string_evolution_mode == 4) {
        // only m_over_sigma fluctuates
        y_loss = sample_rapidity_loss_shell(y_in_lrf);
        m_over_sigma = tau_form / std::max(eps, sqrt(2. * (cosh(y_loss) - 1.)));
    } else if (string_evolution_mode == -4) {
        // only m_over_sigma fluctuates for Beam Remnants
        y_loss =
            (sample_rapidity_loss_shell(y_in_lrf)
             * remnant_energy_loss_fraction_);
        m_over_sigma = tau_form / std::max(eps, sqrt(2. * (cosh(y_loss) - 1.)));
    }
}

// add flag if baryon used to proj and target
// - check it - only put in baryon if not used yet, then set it to used
int Glauber::perform_string_production() {
    QCD_string_list.clear();
    remnant_string_list_.clear();
    const auto string_evolution_mode =
        (parameter_list.get_QCD_string_evolution_mode());
    const auto baryon_junctions = parameter_list.get_baryon_junctions();

    real lambdaB = parameter_list.get_lambdaB();
    lambdaB = std::min(1., lambdaB);
    real lambdaBs = parameter_list.get_lambdaBs();
    lambdaBs = std::min(1., lambdaBs);

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
            real y_in_lrf =
                std::abs(proj->get_rapidity() - targ->get_rapidity()) / 2.;
            std::shared_ptr<Quark> proj_q;
            std::shared_ptr<Quark> targ_q;
            if (sample_valence_quark) {
                auto proj_xvec = proj->get_x();
                auto targ_xvec = targ->get_x();
                proj_q = proj->get_a_valence_quark();
                auto proj_q_xvec = proj_q->get_x();
                if (proj_q->get_number_of_connections() == 1) {
                    // first time pick-up the valence quark
                    // we need to substract the valence quark energy-momentum
                    // from the nucleon remnant energy-momentum vector
                    auto p_q = proj_q->get_p();
                    proj->substract_momentum_from_remnant(p_q);
                }
                // targ_q = targ->get_a_valence_quark();
                targ_q = targ->get_a_close_valence_quark(
                    proj_q_xvec[1] + proj_xvec[1] - targ_xvec[1],
                    proj_q_xvec[2] + proj_xvec[2] - targ_xvec[2]);
                if (targ_q->get_number_of_connections() == 1) {
                    // first time pick-up the valence quark
                    // we need to substract the valence quark energy-momentum
                    // from the nucleon remnant energy-momentum vector
                    auto p_q = targ_q->get_p();
                    targ->substract_momentum_from_remnant(p_q);
                }
                y_in_lrf =
                    std::abs(proj_q->get_rapidity() - targ_q->get_rapidity())
                    / 2.;
            }
            real tau_form = 0.5;      // [fm]
            real m_over_sigma = 1.0;  // [fm]
            real y_loss = 0.;
            get_tau_form_and_moversigma(
                string_evolution_mode, y_in_lrf, tau_form, m_over_sigma,
                y_loss);
            // set variables in case of no baryon junction transport
            bool has_baryon_left = false;
            bool has_baryon_right = false;
            if (!sample_valence_quark) {
                QCDString qcd_string(
                    x_coll, tau_form, proj, targ, m_over_sigma,
                    has_baryon_right, has_baryon_left);
                QCD_string_list.push_back(qcd_string);
            } else {
                // Connect the closest quark pairs into string
                auto proj_q_xvec = proj_q->get_x();
                auto targ_q_xvec = targ_q->get_x();
                SpatialVec x_coll_q = {
                    x_coll[0],
                    x_coll[1] + (proj_q_xvec[1] + targ_q_xvec[1]) / 2.,
                    x_coll[2] + (proj_q_xvec[2] + targ_q_xvec[2]) / 2.,
                    x_coll[3]};
                QCDString qcd_string(
                    x_coll_q, tau_form, proj, targ, proj_q, targ_q,
                    m_over_sigma, has_baryon_right, has_baryon_left);
                QCD_string_list.push_back(qcd_string);
            }
            real y_shift = y_loss;
            if (!sample_valence_quark) {
                update_momentum(proj, -y_shift);
                update_momentum(targ, y_shift);
            } else {
                // shift the nucleon rapidity to avoid identical collision when
                // update the collision schedule
                update_momentum(proj, -1e-3 * proj->get_rapidity());
                update_momentum(targ, 1e-3 * targ->get_rapidity());
                update_momentum_quark(proj_q, -y_shift);
                update_momentum_quark(targ_q, y_shift);
            }
        }
        update_collision_schedule(first_event);
        collision_schedule.erase((*collision_schedule.begin()));
    }

    // randomize the QCD_string_list and assign the baryon charge to
    // the strings
    std::vector<unsigned int> random_idx;
    const real baryonInStringProb = parameter_list.get_baryon_in_string_prob();
    unsigned int Nstrings = QCD_string_list.size();
    unsigned int Npart_proj = projectile->get_number_of_wounded_nucleons();
    unsigned int Npart_targ = target->get_number_of_wounded_nucleons();
    unsigned int total_length = Nstrings + Npart_proj + Npart_targ;
    for (unsigned int idx = 0; idx < total_length; idx++)
        random_idx.push_back(idx);
    std::shuffle(
        random_idx.begin(), random_idx.end(), *ran_gen_ptr_->getRanGenerator());
    for (auto &idx : random_idx) {
        if (idx < Nstrings) {
            // put baryon of the projectile in the selected string
            auto proj = QCD_string_list[idx].get_proj();
            if (proj->get_baryon_number() == 0) proj->set_baryon_used(true);
            if (!proj->baryon_was_used()) {
                if (ran_gen_ptr_->rand_uniform() < baryonInStringProb) {
                    proj->set_baryon_used(true);
                    QCD_string_list[idx].set_has_baryon_right(true);
                }
            }
        } else if (idx < Nstrings + Npart_proj) {
            // put baryon of the projectile in the projectile remnant
            auto proj = projectile->get_participant(idx - Nstrings);
            auto p_i = proj->get_remnant_p();
            if (p_i[0] <= 0) continue;
            // auto mass = 0.;
            // if (std::abs(p_i[3]) < p_i[0]) {
            //     // a time-like beam remnant
            //     mass = sqrt(p_i[0]*p_i[0] - p_i[3]*p_i[3]);
            // }
            // if (!proj->baryon_was_used() && mass > 0.1) {}
            if (proj->get_baryon_number() == 0) proj->set_baryon_used(true);
            if (!proj->baryon_was_used()) {
                proj->set_baryon_used(true);
                proj->set_remnant_carry_baryon_number(true);
            }
        }
    }
    std::shuffle(
        random_idx.begin(), random_idx.end(), *ran_gen_ptr_->getRanGenerator());
    for (auto &idx : random_idx) {
        if (idx < Nstrings) {
            // put baryon of the target in the selected string
            auto targ = QCD_string_list[idx].get_targ();
            if (targ->get_baryon_number() == 0) targ->set_baryon_used(true);
            if (!targ->baryon_was_used()) {
                if (ran_gen_ptr_->rand_uniform() < baryonInStringProb) {
                    targ->set_baryon_used(true);
                    QCD_string_list[idx].set_has_baryon_left(true);
                }
            }
        } else if (idx > Nstrings + Npart_proj - 1) {
            // put baryon of the target in the target remnant
            auto targ = target->get_participant(idx - Nstrings - Npart_proj);
            auto p_i = targ->get_remnant_p();
            if (p_i[0] <= 0) continue;
            // auto mass = 0.;
            // if (std::abs(p_i[3]) < p_i[0]) {
            //     // a time-like beam remnant
            //     mass = sqrt(p_i[0]*p_i[0] - p_i[3]*p_i[3]);
            // }
            // if (!targ->baryon_was_used() && mass > 0.1) {}
            if (targ->get_baryon_number() == 0) targ->set_baryon_used(true);
            if (!targ->baryon_was_used()) {
                targ->set_baryon_used(true);
                targ->set_remnant_carry_baryon_number(true);
            }
        }
    }

    // set baryons' rapidities
    for (auto &it : QCD_string_list) {
        it.evolve_QCD_string();
        if (!baryon_junctions) {
            // set baryon rapidities to string endpoint rapidities
            // if no junction transport is used
            it.set_final_baryon_rapidities(
                it.get_y_f_left(), it.get_y_f_right());
        } else {
            // sample HERE if baryon should be moved
            real y_baryon_right = 0.;
            if (it.get_has_baryon_right()) {
                if (ran_gen_ptr_->rand_uniform() < lambdaB) {
                    if (ran_gen_ptr_->rand_uniform() < lambdaBs) {
                        // y_baryon_right = sample_junction_rapidity_right(
                        //          it->get_y_i_left(), it->get_y_i_right());
                        y_baryon_right = sample_junction_rapidity_right(
                            it.get_y_f_left(), it.get_y_f_right());
                    } else {
                        y_baryon_right = sample_junction_rapidity_uniformed(
                            it.get_y_f_left(), it.get_y_f_right());
                    }
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
                    if (ran_gen_ptr_->rand_uniform() < lambdaBs) {
                        // y_baryon_left = sample_junction_rapidity_left(
                        //           it->get_y_i_left(), it->get_y_i_right());
                        y_baryon_left = sample_junction_rapidity_left(
                            it.get_y_f_left(), it.get_y_f_right());
                    } else {
                        y_baryon_left = sample_junction_rapidity_uniformed(
                            it.get_y_f_left(), it.get_y_f_right());
                    }
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
            SpatialVec x_frez_proj = {0., 0., 0., 0.};
            SpatialVec x_frez_targ = {0., 0., 0., 0.};
            if (!proj_n->is_remnant_set()) {
                proj_n->set_remnant(true);
                auto x_frez = proj_n->get_x();
                x_frez_proj = x_frez;
                proj_n->set_remnant_x_frez(x_frez);
            }
            auto targ_n = it->get_targ();
            if (!targ_n->is_remnant_set()) {
                targ_n->set_remnant(true);
                auto x_frez = targ_n->get_x();
                x_frez_targ = x_frez;
                targ_n->set_remnant_x_frez(x_frez);
            }

            // set flags for quark remnants at their last connected strings
            auto proj_q = it->get_proj_q();
            if (!proj_q->is_remnant_set()) {
                proj_q->set_remnant(true);
                it->set_has_remnant_right(true);
            }
            auto targ_q = it->get_targ_q();
            if (!targ_q->is_remnant_set()) {
                targ_q->set_remnant(true);
                it->set_has_remnant_left(true);
            }
        } else {
            auto proj_n = it->get_proj();
            if (!proj_n->is_remnant_set()) {
                proj_n->set_remnant(true);
                it->set_has_remnant_right(true);
            }
            auto targ_n = it->get_targ();
            if (!targ_n->is_remnant_set()) {
                targ_n->set_remnant(true);
                it->set_has_remnant_left(true);
            }
        }
    }

    produce_remnant_strings();
    return (number_of_collided_events);
}

void Glauber::produce_remnant_strings() {
    // create strings for the beam remnants
    const auto string_evolution_mode = -4;
    real tau_form = 0.5;
    real m_over_sigma = 1.0;  // [fm]
    real y_loss = 0.;
    auto proj_nucleon_list = projectile->get_nucleon_list();
    for (auto &iproj : (*proj_nucleon_list)) {
        if (iproj->is_wounded()) {
            auto x_i = iproj->get_remnant_x_frez();

            auto p_i = iproj->get_remnant_p();
            if (p_i[0] <= 0.) continue;
            auto y_rem = ybeam;
            if (std::abs(p_i[3]) < p_i[0]) {
                // a time-like beam remnant
                y_rem = 0.5 * log((p_i[0] + p_i[3]) / (p_i[0] - p_i[3]));
            }
            auto cosh_y_rem = cosh(y_rem);
            auto m_rem = p_i[0] / cosh_y_rem;
            p_i[1] = 0.;
            p_i[2] = 0.;
            p_i[3] = m_rem * sinh(y_rem);
            MomentumVec targ_p_vec = {p_i[0], p_i[1], p_i[2], -p_i[3]};

            get_tau_form_and_moversigma(
                string_evolution_mode, y_rem, tau_form, m_over_sigma, y_loss);
            bool has_baryon_left = false;
            bool has_baryon_right = iproj->is_remnant_carry_baryon_number();
            QCDString qcd_string(
                x_i, tau_form, iproj, iproj, p_i, targ_p_vec, m_over_sigma,
                has_baryon_right, has_baryon_left);
            qcd_string.set_has_remnant_right(true);
            qcd_string.evolve_QCD_string();
            qcd_string.set_final_baryon_rapidities(0., y_rem - y_loss);
            remnant_string_list_.push_back(qcd_string);
        }
    }
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &itarg : (*targ_nucleon_list)) {
        if (itarg->is_wounded()) {
            auto x_i = itarg->get_remnant_x_frez();

            auto p_i = itarg->get_remnant_p();
            if (p_i[0] <= 0.) continue;
            auto y_rem = -ybeam;
            if (std::abs(p_i[3]) < p_i[0]) {
                // a time-like beam remnant
                y_rem = 0.5 * log((p_i[0] + p_i[3]) / (p_i[0] - p_i[3]));
            }
            auto cosh_y_rem = cosh(y_rem);
            auto m_rem = p_i[0] / cosh_y_rem;
            p_i[1] = 0.;
            p_i[2] = 0.;
            p_i[3] = m_rem * sinh(y_rem);
            MomentumVec proj_p_vec = {p_i[0], p_i[1], p_i[2], -p_i[3]};

            get_tau_form_and_moversigma(
                string_evolution_mode, std::abs(y_rem), tau_form, m_over_sigma,
                y_loss);
            bool has_baryon_left = itarg->is_remnant_carry_baryon_number();
            bool has_baryon_right = false;
            QCDString qcd_string(
                x_i, tau_form, itarg, itarg, proj_p_vec, p_i, m_over_sigma,
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
    for (auto &it : (*proj->get_collide_nucleon_list()))
        create_a_collision_event(proj, it.lock());
    auto targ = event_happened->get_targ_nucleon_ptr().lock();
    targ->increment_collided_times();
    for (auto &it : (*targ->get_collide_nucleon_list()))
        create_a_collision_event(it.lock(), targ);
}

void Glauber::computeCenterOfMass(real &x_o, real &y_o) {
    // compute the center of mass
    x_o = 0.;
    y_o = 0.;
    int n_strings = QCD_string_list.size();
    for (auto &it : QCD_string_list) {
        auto x_prod = it.get_x_production();
        x_o += x_prod[1];
        y_o += x_prod[2];
    }
    if (sample_valence_quark) {
        for (auto &it : remnant_string_list_) {
            auto x_prod = it.get_x_production();
            x_o += x_prod[1];
            y_o += x_prod[2];
        }
        n_strings += remnant_string_list_.size();
    }
    x_o /= n_strings;
    y_o /= n_strings;
}

void Glauber::prepare_output_QCD_strings() {
    QCD_string_output_arr_.clear();
    // compute the center of mass
    real x_o = 0.;
    real y_o = 0.;
    computeCenterOfMass(x_o, y_o);

    // prepare output strings
    for (auto &it : QCD_string_list) {
        auto x_prod = it.get_x_production();
        auto x_left = it.get_targ()->get_x();
        auto x_right = it.get_proj()->get_x();
        if (sample_valence_quark) {
            auto xq_left = it.get_targ_q()->get_x();
            x_left[1] += xq_left[1];
            x_left[2] += xq_left[2];

            auto xq_right = it.get_proj_q()->get_x();
            x_right[1] += xq_right[1];
            x_right[2] += xq_right[2];
        }
        auto tau_0 = sqrt(x_prod[0] * x_prod[0] - x_prod[3] * x_prod[3]);
        auto etas_0 =
            0.5 * log((x_prod[0] + x_prod[3]) / (x_prod[0] - x_prod[3]));

        real remnant_left = 0.;
        if (it.get_has_remnant_left()) {
            remnant_left = 1.0;
        }

        real remnant_right = 0.;
        if (it.get_has_remnant_right()) {
            remnant_right = 1.0;
        }

        real baryon_fraction_left = 0.;
        if (it.get_has_baryon_left()) {
            baryon_fraction_left = 1.0;
        }

        real baryon_fraction_right = 0.;
        if (it.get_has_baryon_right()) {
            baryon_fraction_right = 1.0;
        }

        auto mass = it.get_string_mass();
        std::vector<real> output_array = {
            mass,
            it.get_m_over_sigma(),
            it.get_tau_form(),
            tau_0,
            etas_0,
            x_prod[1] - x_o,
            x_prod[2] - y_o,
            x_left[1] - x_o,
            x_left[2] - y_o,
            x_right[1] - x_o,
            x_right[2] - y_o,
            it.get_eta_s_left(),
            it.get_eta_s_right(),
            it.get_y_f_left(),
            it.get_y_f_right(),
            remnant_left,
            remnant_right,
            it.get_y_i_left(),
            it.get_y_i_right(),
            it.get_eta_s_baryon_left(),
            it.get_eta_s_baryon_right(),
            it.get_y_f_baryon_left(),
            it.get_y_f_baryon_right(),
            baryon_fraction_left,
            baryon_fraction_right,
        };
        QCD_string_output_arr_.push_back(output_array);
    }

    // output the beam remnant strings
    if (sample_valence_quark) {
        for (auto &it : remnant_string_list_) {
            auto x_prod = it.get_x_production();
            auto x_left = it.get_targ()->get_x();
            auto x_right = it.get_proj()->get_x();
            auto tau_0 = sqrt(x_prod[0] * x_prod[0] - x_prod[3] * x_prod[3]);
            auto etas_0 =
                0.5 * log((x_prod[0] + x_prod[3]) / (x_prod[0] - x_prod[3]));

            real remnant_left = 0.;
            if (it.get_has_remnant_left()) {
                remnant_left = 1.0;
            }

            real remnant_right = 0.;
            if (it.get_has_remnant_right()) {
                remnant_right = 1.0;
            }

            real baryon_fraction_left = 0.;
            if (it.get_has_baryon_left()) {
                baryon_fraction_left = 1.0;
            }

            real baryon_fraction_right = 0.;
            if (it.get_has_baryon_right()) {
                baryon_fraction_right = 1.0;
            }

            auto mass = it.get_string_mass();
            auto eta_s_center =
                (it.get_eta_s_left() + it.get_eta_s_right()) / 2.;
            auto eta_s_left =
                (remnant_left * it.get_eta_s_left()
                 + (1. - remnant_left)
                       * (std::max(
                           (eta_s_center + it.get_eta_s_right()) / 2.,
                           it.get_eta_s_right() - 1.0)));
            auto eta_s_right =
                (remnant_right * it.get_eta_s_right()
                 + (1. - remnant_right)
                       * (std::min(
                           (eta_s_center + it.get_eta_s_left()) / 2.,
                           it.get_eta_s_left() + 1.0)));
            std::vector<real> output_array = {
                mass,
                it.get_m_over_sigma(),
                it.get_tau_form(),
                tau_0,
                etas_0,
                x_prod[1] - x_o,
                x_prod[2] - y_o,
                x_left[1] - x_o,
                x_left[2] - y_o,
                x_right[1] - x_o,
                x_right[2] - y_o,
                eta_s_left,
                eta_s_right,
                remnant_left * it.get_y_f_left()
                    + (1. - remnant_left) * eta_s_left,
                remnant_right * it.get_y_f_right()
                    + (1. - remnant_right) * eta_s_right,
                remnant_left,
                remnant_right,
                remnant_left * it.get_y_i_left()
                    + (1. - remnant_left) * eta_s_left,
                remnant_right * it.get_y_i_right()
                    + (1. - remnant_right) * eta_s_right,
                remnant_left * it.get_eta_s_baryon_left(),
                remnant_right * it.get_eta_s_baryon_right(),
                remnant_left * it.get_y_f_baryon_left(),
                remnant_right * it.get_y_f_baryon_right(),
                baryon_fraction_left,
                baryon_fraction_right,
            };
            QCD_string_output_arr_.push_back(output_array);
        }
    }
}

void Glauber::output_QCD_strings(
    std::string filename, const real Npart, const real Ncoll,
    const real Nstrings, const real b, const unsigned int seed) {
    if (QCD_string_output_arr_.size() == 0) prepare_output_QCD_strings();

    std::ofstream output(filename.c_str());
    real total_energy = Npart * parameter_list.get_roots() / 2.;
    real net_Pz =
        ((projectile->get_number_of_wounded_nucleons()
          - target->get_number_of_wounded_nucleons())
         * parameter_list.get_roots() / 2.);
    output << "# b = " << b << " fm " << "Npart = " << Npart
           << " Ncoll = " << Ncoll << " Nstrings = " << Nstrings
           << " total_energy = " << total_energy << " GeV, "
           << "net_Pz = " << net_Pz << " GeV, "
           << "seed = " << seed << std::endl;

    output
        << "# mass[GeV]  m_over_sigma[fm]  tau_form[fm]  tau_0[fm]  eta_s_0  "
        << "x_perp[fm]  y_perp[fm]  x_l[fm]  y_l[fm]  x_r[fm]  y_r[fm]  "
        << "eta_s_left  eta_s_right  y_l  y_r  remnant_l  remnant_r "
        << "y_l_i  y_r_i "
        << "eta_s_baryon_left  eta_s_baryon_right  y_l_baryon  y_r_baryon  "
        << "baryon_fraction_l  baryon_fraction_r" << std::endl;

    for (auto &string_i : QCD_string_output_arr_) {
        output << std::scientific << std::setprecision(8);
        for (auto &ival : string_i) {
            output << std::setw(15) << ival << "  ";
        }
        output << std::endl;
    }
    output.close();
}

void Glauber::prepareParticipantList() {
    participantList_.clear();
    // compute the center of mass
    real x_o = 0.;
    real y_o = 0.;
    computeCenterOfMass(x_o, y_o);
    auto proj_nucleon_list = projectile->get_nucleon_list();
    int dir = 1;
    for (auto &iproj : (*proj_nucleon_list)) {
        if (iproj->is_wounded()) {
            auto proj_x = iproj->get_x();
            proj_x[1] -= x_o;
            proj_x[2] -= y_o;
            std::vector<real> part_i;
            part_i.push_back(proj_x[0]);
            part_i.push_back(proj_x[1]);
            part_i.push_back(proj_x[2]);
            part_i.push_back(proj_x[3]);
            part_i.push_back(dir);
            part_i.push_back(iproj->get_electric_charge());
            participantList_.push_back(part_i);
        }
    }
    dir = -1;
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &itarg : (*targ_nucleon_list)) {
        if (itarg->is_wounded()) {
            auto targ_x = itarg->get_x();
            targ_x[1] -= x_o;
            targ_x[2] -= y_o;
            std::vector<real> part_i;
            part_i.push_back(targ_x[0]);
            part_i.push_back(targ_x[1]);
            part_i.push_back(targ_x[2]);
            part_i.push_back(targ_x[3]);
            part_i.push_back(dir);
            part_i.push_back(itarg->get_electric_charge());
            participantList_.push_back(part_i);
        }
    }
}

void Glauber::outputParticipants(std::string filename) {
    prepareParticipantList();

    std::ofstream output(filename.c_str());
    output << "# t[fm]  x[fm]  y[fm]  z[fm]  dir  e" << std::endl;
    for (auto &ipart : participantList_) {
        output << std::scientific << std::setprecision(6);
        for (auto &x_i : ipart) {
            output << std::setw(10) << x_i << "  ";
        }
        output << std::endl;
    }
    output.close();
}

void Glauber::output_spectators(std::string filename) {
    // compute the center of mass
    real x_o = 0.;
    real y_o = 0.;
    computeCenterOfMass(x_o, y_o);

    std::ofstream output(filename.c_str());
    output << "# t[fm]  x[fm]  y[fm]  z[fm]  m[GeV]  px[GeV]  py[GeV]  y  e"
           << std::endl;
    auto proj_nucleon_list = projectile->get_nucleon_list();
    for (auto &iproj : (*proj_nucleon_list)) {
        if (!iproj->is_wounded()) {
            output << std::scientific << std::setprecision(6);
            auto proj_x = iproj->get_x();
            auto proj_p = iproj->get_p();
            auto mass = iproj->get_mass();

            // compute the when and where the spectator nucleon enters the
            // light cone
            real vz = proj_p[3] / proj_p[0];
            proj_x[3] = vz / (1. + vz) * (proj_x[3] / vz - proj_x[0]);
            proj_x[0] = -proj_x[3];
            proj_x[1] -= x_o;
            proj_x[2] -= y_o;

            // output spectator's position and momentum
            for (const auto &x_i : proj_x)
                output << std::setw(10) << x_i << "  ";

            auto fermiMomentum = iproj->get_fermi_momentum();
            for (int i = 1; i < 4; i++) proj_p[i] += fermiMomentum[i];
            auto rap = asinh(
                proj_p[3]
                / (sqrt(
                    mass * mass + proj_p[1] * proj_p[1]
                    + proj_p[2] * proj_p[2])));
            output << std::setw(10) << mass << "  " << std::setw(10)
                   << proj_p[1] << "  " << std::setw(10) << proj_p[2] << "  "
                   << std::setw(10) << rap << "  "
                   << iproj->get_electric_charge() << std::endl;
        }
    }
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &itarg : (*targ_nucleon_list)) {
        if (!itarg->is_wounded()) {
            output << std::scientific << std::setprecision(6);
            auto targ_x = itarg->get_x();
            auto targ_p = itarg->get_p();
            auto mass = itarg->get_mass();

            // compute the when and where the spectator nucleon enters the
            // light cone
            real vz = targ_p[3] / targ_p[0];
            targ_x[3] = vz / (1. - vz) * (targ_x[3] / vz - targ_x[0]);
            targ_x[0] = targ_x[3];
            targ_x[1] -= x_o;
            targ_x[2] -= y_o;

            // output spectator's position and momentum
            for (const auto &x_i : targ_x)
                output << std::setw(10) << x_i << "  ";

            auto fermiMomentum = itarg->get_fermi_momentum();
            for (int i = 1; i < 4; i++) targ_p[i] += fermiMomentum[i];
            auto rap = asinh(
                targ_p[3]
                / (sqrt(
                    mass * mass + targ_p[1] * targ_p[1]
                    + targ_p[2] * targ_p[2])));
            output << std::setw(10) << mass << "  " << std::setw(10)
                   << targ_p[1] << "  " << std::setw(10) << targ_p[2] << "  "
                   << std::setw(10) << rap << "  "
                   << itarg->get_electric_charge() << std::endl;
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
    } else if (parameter_list.get_rapidity_loss_method() > 2) {
        y_loss = sample_rapidity_loss_from_parametrization_with_fluct(y_init);
    }
    return (y_loss);
}

real Glauber::sample_rapidity_loss_from_the_LEXUS_model(
    const real y_init) const {
    const real shape_coeff = 1.0;
    real sinh_y_lrf = sinh(shape_coeff * y_init);
    real arcsinh_factor =
        (ran_gen_ptr_->rand_uniform()
             * (sinh(2. * shape_coeff * y_init) - sinh_y_lrf)
         + sinh_y_lrf);
    real y_loss = 2. * y_init - asinh(arcsinh_factor) / shape_coeff;
    return (y_loss);
}

real Glauber::sample_rapidity_loss_from_parametrization(
    const real y_init) const {
    auto y_loss =
        (yloss_param_slope
         * pow(pow(y_init, yloss_param_a) * tanh(y_init), yloss_param_b));
    return (y_loss);
}

real Glauber::sample_rapidity_loss_from_piecewise_parametrization(
    const real y_init) const {
    auto y_loss1 = parameter_list.getParam("ylossParam4At2", 1.60);
    auto y_loss2 = parameter_list.getParam("ylossParam4At4", 2.15);
    auto y_loss3 = parameter_list.getParam("ylossParam4At6", 2.45);
    auto y_loss4 = parameter_list.getParam("ylossParam4At10", 2.95);

    real y_loss = 0.;
    if (y_init < 2) {
        y_loss = y_loss1 / 2. * y_init;
    } else if (y_init < 4) {
        y_loss = (y_loss2 - y_loss1) / 2. * y_init + (2. * y_loss1 - y_loss2);
    } else if (y_init < 6) {
        y_loss =
            (y_loss3 - y_loss2) / 2. * y_init + (3. * y_loss2 - 2. * y_loss3);
    } else {
        y_loss =
            (y_loss4 - y_loss3) / 4. * y_init + (2.5 * y_loss3 - 1.5 * y_loss4);
    }
    return (y_loss);
}

real Glauber::sample_rapidity_loss_from_parametrization_with_fluct(
    const real y_init) const {
    real y_mean = y_init / 2.;
    real var = 0.5;
    if (parameter_list.get_rapidity_loss_method() == 3) {
        y_mean = sample_rapidity_loss_from_parametrization(y_init);
        auto y_rhic = 5.5;
        auto y_lhc = 9.0;
        auto var_rhic = parameter_list.get_yloss_param_fluct_var_RHIC();
        auto var_lhc = parameter_list.get_yloss_param_fluct_var_LHC();
        var = var_rhic;
        if (y_init > y_lhc) {
            var = var_lhc;
        } else if (y_init > y_rhic) {
            var =
                (var_rhic
                 + (y_init - y_rhic) / (y_lhc - y_rhic) * (var_lhc - var_rhic));
        }
    } else if (parameter_list.get_rapidity_loss_method() == 4) {
        y_mean = sample_rapidity_loss_from_piecewise_parametrization(y_init);
        var = parameter_list.getParam("ylossParam4var", 0.5);
    }

    // sample logit distribution for y_loss given the mean and variance
    auto random_x = ran_gen_ptr_->rand_normal(0., var);
    real logit_rand = 1. / (1. + exp(-random_x));
    auto y_loss = y_mean;
    if (std::abs(2. * y_mean - y_init) > 1e-15) {
        real aa =
            (2. * y_mean - y_init) / (2. * y_mean * y_init * (y_init - y_mean));
        real bb =
            ((y_init * y_init - 2. * y_mean * y_mean)
             / (2. * y_init * y_mean * (y_init - y_mean)));
        real cc = -logit_rand;
        y_loss = (-bb + sqrt(bb * bb - 4. * aa * cc)) / (2. * aa);
    } else {
        y_loss = logit_rand * y_init;
    }
    return (y_loss);
}

// sample y from exp[(y - (0.5*(yt + yp)))/2]/(4.*Sinh[0.25*yp - 0.25*yt]),
// the new rapidity of the baryon number from the right moving particle
// after the collision in the lab frame
real Glauber::sample_junction_rapidity_right(
    const real y_left, const real y_right) const {
    real y = -2.
             * (-0.25 * y_right - 0.25 * y_left
                - 1.
                      * log(
                          2.
                          * (ran_gen_ptr_->rand_uniform()
                             + 0.5 * exp(-0.25 * y_right + 0.25 * y_left)
                                   / sinh(0.25 * y_right - 0.25 * y_left))
                          * sinh(0.25 * y_right - 0.25 * y_left)));
    return (y);
}

// sample y from exp[-(y - (0.5*(yt + yp)))/2]/(4.*Sinh[0.25*yp - 0.25*yt]),
// the new rapidity of the baryon number from the left moving particle
// after the collision in the lab frame
real Glauber::sample_junction_rapidity_left(
    const real y_left, const real y_right) const {
    real y = 2.
             * (0.25 * y_right + 0.25 * y_left
                - log(
                    2.
                    * (ran_gen_ptr_->rand_uniform()
                       + 0.5 * exp(-0.25 * y_right + 0.25 * y_left)
                             / sinh(0.25 * y_right - 0.25 * y_left))
                    * sinh(0.25 * y_right - 0.25 * y_left)));
    return (y);
}

// sample y from a uniformed 1/(yp - yt) distribution
real Glauber::sample_junction_rapidity_uniformed(
    const real y_left, const real y_right) const {
    real y = y_left + ran_gen_ptr_->rand_uniform() * (y_right - y_left);
    return (y);
}

// This function computes the sigeff(s) from the formula
//  sigmaNN_in(s) = int d^2b [1 - exp(-sigeff(s)*Tpp(b))]
//  Reads sigmaNN, returns guassian width
real Glauber::get_sig_eff(const real siginNN) {
    // rms-radius of a gaussian = rms-radius of a disc with radius R,
    // where 2*PI*(2R)^2=sigmaNN
    const real width = sqrt(0.1 * siginNN / M_PI) / sqrt(8.);
    nucleon_width_ = width;
    const real sigin = siginNN * 0.1;  // sigma_in(s) [mb --> fm^2]

    const int Nint = 100;  // # of integration points
    std::vector<real> b(Nint, 0.);
    std::vector<real> Tnn(Nint, 0.);
    const real Bmax = 10.0 * width;
    const real db = Bmax / Nint;
    for (int i = 0; i < Nint; i++) {
        b[i] = (i + 0.5) * db;
        Tnn[i] = exp(-b[i] * b[i] / (4. * width * width))
                 / (M_PI * (4. * width * width));
    }
    const real prefactor = 2. * M_PI * db;

    real sum, dN;

    real sigeff = 10.0;  // starting point of iteration [fm^2]
    real sigeff0;        // holds value from previous iteration step
    int tol = 0;
    do {  // iterate ...
        sigeff0 = sigeff;
        sum = 0.0;
        dN = 0.0;
        for (int ib = 0; ib < Nint; ib++) {  // integral d^2b from 0 to Bmax
            sum += prefactor * b[ib] * (1.0 - exp(-sigeff * Tnn[ib]));
            dN += prefactor * b[ib] * Tnn[ib] * exp(-sigeff * Tnn[ib]);
        }
        sigeff -= (sum - sigin) / dN;
        tol++;
        // cout << "iter: " << tol << ": sigeff = " << sigeff
        //      << " fm^2, sum = " << sum
        //      << " fm^2, sigin = " << sigin << " fm^2" << endl;
    } while (std::abs(sigeff - sigeff0) > 1e-4 && tol < 1000);
    // until sigeff has converged

    return (sigeff);
}

}  // namespace MCGlb
