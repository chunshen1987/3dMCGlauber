// Copyright @ Chun Shen 2018

#include "EventGenerator.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

namespace MCGlb {

std::vector<CollisionEvent> EventGenerator::get_CollisionEventvector() {
    return mc_glauber_ptr_->get_collision_information();
}


EventGenerator::EventGenerator(std::string input_filename,
                               int argc, char* argv[], int seed) {
    parameter_list_.read_in_parameters_from_file(input_filename);
    parameter_list_.read_in_parameters_from_arguments(argc, argv, "=", 4);
    parameter_list_.print_parameter_list();
    int ran_seed = parameter_list_.get_seed();
    if (seed != 0) ran_seed = seed;
    auto gamma_beta = parameter_list_.get_tau_form_fluct_gamma_beta();
    ran_gen_ptr_ = std::shared_ptr<RandomUtil::Random>(
                    new RandomUtil::Random(ran_seed, 0.0, 1.0, gamma_beta));
    mc_glauber_ptr_ = std::unique_ptr<Glauber>(
                    new Glauber(parameter_list_, ran_gen_ptr_));

    statistics_only_ = parameter_list_.get_only_event_statistics();
    ecm_ = parameter_list_.get_roots();

    batchDensityOutput_ = parameter_list_.get_batch_density_output();
    batchEccOutput_ = parameter_list_.get_batch_eccentricity_output();
    initialEstOutput_ = parameter_list_.get_initialEst_output();
    density_maker_ptr_ = std::unique_ptr<MakeDensity>(new MakeDensity());
    // set grid information
    density_maker_ptr_->set_1D_grid_info_eta(72, 0.2);
    density_maker_ptr_->set_2D_grid_info_etax(72, 0.2, 100, 0.2);
    density_maker_ptr_->set_3D_grid_info(100, 0.2, 100, 0.2, 72, 0.2);

    // set Gaussian widths for sigma_x, sigma_eta
    density_maker_ptr_->setGaussianWidths(0.4, 0.5);
    density_maker_ptr_->setStringTransShiftFrac(0.0);
    if (batchDensityOutput_) {
        batchDensity2DOutput_ = parameter_list_.get_batch_2Ddensity_output();
    } else {
        batchDensity2DOutput_ = false;
    }
    cenEstMax_ = 1e16;
    cenEstMin_ = 0;
}


void EventGenerator::set_parameter(std::string parameter_name, double value) {
    parameter_list_.set_parameter(parameter_name, value);
}



void EventGenerator::generate_pre_events() {
    //messager << "Random seed = " << ran_gen_ptr_->get_seed();
    //messager.flush("info");
    //messager << "Generating 1 events ... ";
    //messager.flush("info");
    int iev =0;
    while (iev < 1) {
        mc_glauber_ptr_->make_nuclei();
        Ncoll_ = mc_glauber_ptr_->make_collision_schedule();
        Npart_ = mc_glauber_ptr_->get_Npart();
        auto Nstrings = mc_glauber_ptr_->decide_QCD_strings_production();
        if (event_of_interest_trigger(Npart_, Ncoll_, Nstrings))iev++;
    }
}


double EventGenerator::MCGlb_nucleon_density(double t, double x,
                                             double y, double z) {
    return(mc_glauber_ptr_->get_nucleon_density(t, x, y, z));
}


double EventGenerator::MCGlb_target_nucleon_density(double t, double x,
                                                    double y, double z) {
    return(mc_glauber_ptr_->get_targ_nucleon_density(t, x, y, z));
}


double EventGenerator::MCGlb_projectile_nucleon_density(double t, double x,
                                                        double y, double z) {
    return(mc_glauber_ptr_->get_proj_nucleon_density(t, x, y, z));
}

std::vector<double> EventGenerator::MCGlb_projectile_nucleon_z() {
    return(mc_glauber_ptr_->get_all_proj_nucleon_z());
}

std::vector<double> EventGenerator::MCGlb_Proj_hot_spot_x() {
    auto Hot_spot_pos = mc_glauber_ptr_->get_Proj_hot_spot_x();
    std::vector<double> hot_spot_pos = {Hot_spot_pos[0], Hot_spot_pos[1], 
                                        Hot_spot_pos[2], Hot_spot_pos[3]};
    return(hot_spot_pos);
}

std::vector<double> EventGenerator::MCGlb_Targ_hot_spot_x() {
    auto Hot_spot_pos = mc_glauber_ptr_->get_Targ_hot_spot_x();
    std::vector<double> hot_spot_pos = {Hot_spot_pos[0], Hot_spot_pos[1], 
                                        Hot_spot_pos[2], Hot_spot_pos[3]};
    return(hot_spot_pos);
}

std::vector<double> EventGenerator::MCGlb_target_nucleon_z() {
    return(mc_glauber_ptr_->get_all_targ_nucleon_z());
}

void EventGenerator::GetMomandPos_Proj(double &t, double &x, double &y, double &z,
                                       double &E, double &px, double &py, double &pz) {
    proj_t = t;
    proj_x = x;
    proj_y = y;
    proj_z = z;
    proj_E = E;
    proj_px = px;
    proj_py = py;
    proj_pz = pz;
}

void EventGenerator::GetMomandPos_Targ(double &t, double &x, double &y, double &z,
                                       double &E, double &px, double &py, double &pz) {
    targ_t = t;
    targ_x = x;
    targ_y = y;
    targ_z = z;
    targ_E = E;
    targ_px = px;
    targ_py = py;
    targ_pz = pz;
}

void EventGenerator::GetMom_Proj(double &E, double &px, double &py, double &pz) {
    proj_E  = E;
    proj_px = px;
    proj_py = py;
    proj_pz = pz;
}

void EventGenerator::GetMom_Targ(double &E, double &px, double &py, double &pz) {
    targ_E  = E;
    targ_px = px;
    targ_py = py;
    targ_pz = pz;
}

void EventGenerator::GetHardPos(double &t, double &x, double &y, double &z) {
    proj_t = t;
    proj_x = x;
    proj_y = y;
    proj_z = z;
}

std::vector<double> EventGenerator::GetQuarkPosProj() {
    std::vector<double> HardPartonPos;
    HardPartonPos.push_back(proj_t);
    HardPartonPos.push_back(proj_x);
    HardPartonPos.push_back(proj_y);
    HardPartonPos.push_back(proj_z);

    mc_glauber_ptr_->Set_hard_collisions_Pos(HardPartonPos);
    std::vector<double> quarkxvec = mc_glauber_ptr_->OutputquarkPosProj();
    return(quarkxvec);
}

std::vector<double> EventGenerator::GetQuarkPosTarg() {
    std::vector<double> HardPartonPos;
    HardPartonPos.push_back(proj_t);
    HardPartonPos.push_back(proj_x);
    HardPartonPos.push_back(proj_y);
    HardPartonPos.push_back(proj_z);

    mc_glauber_ptr_->Set_hard_collisions_Pos(HardPartonPos);
    std::vector<double> quarkxvec = mc_glauber_ptr_->OutputquarkPosTarg();
    return(quarkxvec);
}

std::vector<double> EventGenerator::GetRemMom_Proj() {
    std::vector<double> RemMomProj = mc_glauber_ptr_->GetRemMomProj();
    return(RemMomProj);
}

std::vector<double> EventGenerator::GetRemMom_Targ() {
    std::vector<double> RemMomTarg = mc_glauber_ptr_->GetRemMomTarg();
    return(RemMomTarg);
}

void EventGenerator::generate_strings() {
    //messager << "Random seed = " << ran_gen_ptr_->get_seed();
    //messager.flush("info");
    //messager << "Generating 1 events after subtracted four momentum of hard partons ... ";
    //messager.flush("info");
    // this file records all the essential information for the generated events
    std::ofstream record_file("events_summary.dat", std::ios::out);
    record_file << "# event_id  Npart  Ncoll  Nstrings  b(fm)" << std::endl;

    int iev = 0;
    //real mean_Npart = 0;
    HardPartonPosAndMomProj.clear();
    HardPartonPosAndMomProj.push_back(proj_t);
    HardPartonPosAndMomProj.push_back(proj_x);
    HardPartonPosAndMomProj.push_back(proj_y);
    HardPartonPosAndMomProj.push_back(proj_z);
    HardPartonPosAndMomProj.push_back(proj_E);
    HardPartonPosAndMomProj.push_back(proj_px);
    HardPartonPosAndMomProj.push_back(proj_py);
    HardPartonPosAndMomProj.push_back(proj_pz);

    HardPartonPosAndMomTarg.clear();
    HardPartonPosAndMomTarg.push_back(proj_t);
    HardPartonPosAndMomTarg.push_back(proj_x);
    HardPartonPosAndMomTarg.push_back(proj_y);
    HardPartonPosAndMomTarg.push_back(proj_z);
    HardPartonPosAndMomTarg.push_back(targ_E);
    HardPartonPosAndMomTarg.push_back(targ_px);
    HardPartonPosAndMomTarg.push_back(targ_py);
    HardPartonPosAndMomTarg.push_back(targ_pz);

    int event_id = iev; 
    //mean_Npart += Npart_;
    mc_glauber_ptr_->Set_hard_parton_momentum(
                     HardPartonPosAndMomProj, HardPartonPosAndMomTarg);
    mc_glauber_ptr_->Pick_and_subtract_hard_parton_momentum();
    Npart_ = mc_glauber_ptr_->get_Npart();
    auto Nstrings = mc_glauber_ptr_->decide_QCD_strings_production_second_stage();
    Ncoll_ = mc_glauber_ptr_->perform_string_production();
    auto b = mc_glauber_ptr_->get_impact_parameter();
    if (!statistics_only_) {
        std::ostringstream filename;
        filename << "strings_event_" << event_id << ".dat";
        mc_glauber_ptr_->output_QCD_strings(
                        filename.str(), Npart_, Ncoll_, Nstrings, b,
                        ran_gen_ptr_->get_seed());
    }

    // write event information to the record file
    record_file << event_id << "  " << Npart_ << "  " << Ncoll_ << "  "
                << Nstrings << "  " << b << std::endl;
    record_file.close();
    //mean_Npart = static_cast<real>(mean_Npart);
    //messager << "Completed. <Npart> = " << mean_Npart;
    //messager.flush("info");
    //auto b_max = parameter_list_.get_b_max();
    //auto b_min = parameter_list_.get_b_min();
    //auto total_cross_section = (
    //    M_PI*(b_max*b_max - b_min*b_min)/100.);
    //messager << "Total cross section sig_tot = " << total_cross_section
    //         << " b";
    //messager.flush("info");
}


float EventGenerator::computeCenEstimator(const int Npart, const int Ncoll,
                                          const int Nstrings) const {
    return(static_cast<float>(Nstrings));
}


void EventGenerator::generateMinBiasEventList() {
    real cenMin = parameter_list_.getParam("cenMin", 0.);
    real cenMax = parameter_list_.getParam("cenMax", 100.);
    if (cenMin < 1e-8 && 100 - cenMax < 1e-8) return;

    auto b_max_tmp = parameter_list_.get_b_max();
    auto b_min_tmp = parameter_list_.get_b_min();
    messager << "Random seed = " << ran_gen_ptr_->get_seed();
    messager.flush("info");
    const int nev = 50000;
    messager << "Generating minimum bias with " << nev << " events ... ";
    messager.flush("info");
    parameter_list_.set_b_min(0.);
    parameter_list_.set_b_max(25.);
    auto b_max_local = 0.;

    int iev = 0;
    int nev_progress = std::max(1, nev/10);
    while (iev < nev) {
        mc_glauber_ptr_->make_nuclei();
        auto Ncoll = mc_glauber_ptr_->make_collision_schedule();
        auto Npart = mc_glauber_ptr_->get_Npart();

        if (Npart < 2) continue;

        auto b = mc_glauber_ptr_->get_impact_parameter();
        if (b > b_max_local) b_max_local = b;

        if (iev%nev_progress == 0) {
            messager << "Progress: " << iev << " out of " << nev
                      << " is done.";
            messager.flush("info");
            if (iev > 0) {
                parameter_list_.set_b_max(b_max_local + 1.);
            }
        }

        auto Nstrings = mc_glauber_ptr_->decide_QCD_strings_production();
        cenEstMinBiasList_.push_back(computeCenEstimator(Npart, Ncoll,
                                                         Nstrings));
        iev++;
    }
    std::sort(cenEstMinBiasList_.begin(), cenEstMinBiasList_.end(),
              std::greater<int>());

    int idx = std::min(nev - 1, static_cast<int>(nev*cenMax/100.));
    cenEstMin_ = cenEstMinBiasList_[idx];
    idx = std::max(0, static_cast<int>(nev*cenMin/100.));
    cenEstMax_ = cenEstMinBiasList_[idx];
    messager << "centrality cut [" << cenMin << ", " << cenMax
             << "]: cenEstMin = " << cenEstMin_ << ", cenEstMax = "
             << cenEstMax_;
    messager.flush("info");

    parameter_list_.set_b_min(b_min_tmp);
    parameter_list_.set_b_max(b_max_tmp);
}


void EventGenerator::generate_events(int nev, int event_id_offset) {
    //messager << "Random seed = " << ran_gen_ptr_->get_seed();
    //messager.flush("info");
    //messager << "Generating " << nev << " events ... ";
    //messager.flush("info");
    // this file records all the essential information for the generated events
    std::ofstream record_file("events_summary.dat", std::ios::out);
    record_file << "# event_id  Npart  Ncoll  Nstrings  b(fm)" << std::endl;

    int iev = 0;
    int icollisions = 0;
    int nev_progress = std::max(1, nev/10);
    int mean_Npart = 0;
    while (iev < nev) {
        mc_glauber_ptr_->make_nuclei();
        auto Ncoll = mc_glauber_ptr_->make_collision_schedule();
        auto Npart = mc_glauber_ptr_->get_Npart();
        auto Nstrings = mc_glauber_ptr_->decide_QCD_strings_production();
        if (event_of_interest_trigger(Npart, Ncoll, Nstrings))  {
            int event_id = iev + event_id_offset;
            mean_Npart += Npart;
            if (iev%nev_progress == 0) {
                messager << "Progress: " << iev << " out of " << nev
                          << " is done.";
                messager.flush("info");
            }

            Ncoll = mc_glauber_ptr_->perform_string_production();
            auto b = mc_glauber_ptr_->get_impact_parameter();
            // write event information to the record file
            record_file << event_id << "  " << Npart << "  " << Ncoll << "  "
                        << Nstrings << "  " << b << std::endl;
            iev++;
            if (statistics_only_) continue;

            density_maker_ptr_->set_QCD_string_output_arr(
                        mc_glauber_ptr_->get_QCD_strings_output_list());
            if (initialEstOutput_) {
                density_maker_ptr_->output_eccentricity("ecc_ed_n",
                                                        event_id, 1);
                density_maker_ptr_->output_netBaryon_eta_distribution(
                        "nB_etas_distribution", event_id, 1);
                density_maker_ptr_->output_energyDensity_eta_distribution(
                        "ed_etas_distribution", event_id, 1);
            }
            if (batchDensityOutput_) {
                density_maker_ptr_->output_netBaryon_eta_distribution(
                        "nB_etas_distribution", event_id);
                density_maker_ptr_->output_energyDensity_eta_distribution(
                        "ed_etas_distribution", event_id);
                if (batchDensity2DOutput_) {
                    density_maker_ptr_->output_energyDensity_xeta_distribution(
                                        "ed2D_xetas_distribution", event_id);
                    density_maker_ptr_->output_energyDensity_3d(
                                        "ed3D", event_id);
                }
                if (batchEccOutput_) {
                    density_maker_ptr_->output_eccentricity("ecc_ed_n",
                                                            event_id);
                    density_maker_ptr_->setParticipantList(
                                    mc_glauber_ptr_->getParticipantList());
                    density_maker_ptr_->outputTATBEccentricity("ecc_ed_n",
                                                               event_id);
                }
            } else {
                std::ostringstream filename;
                filename << "strings_event_" << event_id << ".dat";
                mc_glauber_ptr_->output_QCD_strings(
                        filename.str(), Npart, Ncoll, Nstrings, b,
                        ran_gen_ptr_->get_seed());
                std::ostringstream specFileName;
                specFileName << "spectators_event_" << event_id << ".dat";
                mc_glauber_ptr_->output_spectators(specFileName.str());

                std::ostringstream partFileName;
                partFileName << "participants_event_" << event_id << ".dat";
                mc_glauber_ptr_->outputParticipants(partFileName.str());
            }
        }
        icollisions++;
    }
    record_file.close();
    mean_Npart = static_cast<real>(mean_Npart)/static_cast<real>(nev);
    messager << "Completed. <Npart> = " << mean_Npart;
    messager.flush("info");

    if (nev > 1000) {
        auto b_max = parameter_list_.get_b_max();
        auto b_min = parameter_list_.get_b_min();
        auto total_cross_section = (
            M_PI*(b_max*b_max - b_min*b_min)*static_cast<real>(nev)
            /static_cast<real>(icollisions)/100.);
        messager << "Total cross section sig_tot = " << total_cross_section
                 << " b";
        messager.flush("info");
    }
}


bool EventGenerator::event_of_interest_trigger(const int Npart,
                                               const int Ncoll,
                                               const int Nstrings) const {
    bool pick = false;
    float cenEst = computeCenEstimator(Npart, Ncoll, Nstrings);
    pick = (Npart > 1) && (cenEst >= cenEstMin_) && (cenEst <= cenEstMax_);
    return(pick);
}


};
