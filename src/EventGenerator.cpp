// Copyright @ Chun Shen 2018

#include "EventGenerator.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

namespace MCGlb {

EventGenerator::EventGenerator(std::string input_filename, int seed_add) {
    parameter_list.read_in_parameters_from_file(input_filename);
    parameter_list.print_parameter_list();
    int seed = parameter_list.get_seed();
    ran_gen_ptr = std::shared_ptr<RandomUtil::Random>(
                            new RandomUtil::Random(seed, 0.0, 1.0, seed_add));
    mc_glauber_ptr = std::unique_ptr<Glauber>(
                                new Glauber(parameter_list, ran_gen_ptr));

    statistics_only = parameter_list.get_only_event_statistics();
}


void EventGenerator::generate_events(int nev, int event_id_offset) {
    messager << "Random seed = " << ran_gen_ptr->get_seed();
    messager.flush("info");
    messager << "Generating " << nev << " events ... ";
    messager.flush("info");
    // this file records all the essential information for the generated events
    std::ofstream record_file("events_summary.dat", std::ios::out);
    record_file << "# event_id  Npart  Ncoll  Nstrings  b(fm)" << std::endl;

    int iev = 0;
    int nev_progress = std::max(1, nev/10);
    int mean_Npart = 0;
    while (iev < nev) {
        mc_glauber_ptr->make_nuclei();
        auto Ncoll = mc_glauber_ptr->make_collision_schedule();
        auto Npart = mc_glauber_ptr->get_Npart();
        auto Nstrings = mc_glauber_ptr->decide_QCD_strings_production();
        if (event_of_interest_trigger(Npart, Ncoll, Nstrings))  {
            int event_id = iev + event_id_offset;
            mean_Npart += Npart;
            if (iev%nev_progress == 0) {
                messager << "Progress: " << iev << " out of " << nev
                          << " is done.";
                messager.flush("info");
            }

            Ncoll = mc_glauber_ptr->perform_string_production();
            auto b = mc_glauber_ptr->get_impact_parameter();
            if (!statistics_only) {
                std::ostringstream filename;
                filename << "strings_event_" << event_id << ".dat";
                mc_glauber_ptr->output_QCD_strings(filename.str(), Npart,
                                                   Ncoll, Nstrings, b);
            }

            // write event information to the record file
            record_file << event_id << "  " << Npart << "  " << Ncoll << "  "
                        << Nstrings << "  " << b << std::endl;
            iev++;
        }
    }
    record_file.close();
    mean_Npart = static_cast<real>(mean_Npart)/static_cast<real>(nev);
    messager << "Completed. <Npart> = " << mean_Npart;
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
