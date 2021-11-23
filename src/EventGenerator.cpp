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

EventGenerator::EventGenerator(std::string input_filename, int seed) {
    parameter_list_.read_in_parameters_from_file(input_filename);
    parameter_list_.print_parameter_list();
    int ran_seed = parameter_list_.get_seed();
    if (seed != 0) ran_seed = seed;
    auto gamma_beta = parameter_list_.get_tau_form_fluct_gamma_beta();
    ran_gen_ptr_ = std::shared_ptr<RandomUtil::Random>(
                    new RandomUtil::Random(ran_seed, 0.0, 1.0, gamma_beta));
    mc_glauber_ptr_ = std::unique_ptr<Glauber>(
                    new Glauber(parameter_list_, ran_gen_ptr_));

    statistics_only_ = parameter_list_.get_only_event_statistics();
}

/* get the collisions information for the JETSCAPE framework*/
void EventGenerator::generate_pre_events(int nev, int event_id_offset) {
    messager << "Random seed = " << ran_gen_ptr_->get_seed();
    messager.flush("info");
    messager << "Generating " << nev << " events ... ";
    messager.flush("info");
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
            iev++;
        }
        icollisions++;
    }
}

void EventGenerator::generate_events(int nev, int event_id_offset) {
    messager << "Random seed = " << ran_gen_ptr_->get_seed();
    messager.flush("info");
    messager << "Generating " << nev << " events ... ";
    messager.flush("info");
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
            if (!statistics_only_) {
                std::ostringstream filename;
                filename << "strings_event_" << event_id << ".dat";
                mc_glauber_ptr_->output_QCD_strings(filename.str(), Npart,
                                                    Ncoll, Nstrings, b);
            }

            // write event information to the record file
            record_file << event_id << "  " << Npart << "  " << Ncoll << "  "
                        << Nstrings << "  " << b << std::endl;
            iev++;
        }
        icollisions++;
    }
    record_file.close();
    mean_Npart = static_cast<real>(mean_Npart)/static_cast<real>(nev);
    messager << "Completed. <Npart> = " << mean_Npart;
    messager.flush("info");
    auto b_max = parameter_list_.get_b_max();
    auto b_min = parameter_list_.get_b_min();
    auto total_cross_section = (
        M_PI*(b_max*b_max - b_min*b_min)*static_cast<real>(nev)
        /static_cast<real>(icollisions)/100.);
    messager << "Total cross section sig_tot = " << total_cross_section
             << " b";
    messager.flush("info");
}



bool EventGenerator::event_of_interest_trigger(int Npart, int Ncoll,
                                               int Nstrings) {
    bool pick = false;
    if (Npart > 1) {
        pick = true;
    }
    return(pick);
}

};
