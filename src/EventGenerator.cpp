// Copyright @ Chun Shen 2018

#include "EventGenerator.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

namespace MCGlb {

EventGenerator::EventGenerator(std::string input_filename, int nev_in) {
    nev = nev_in;
    parameter_list.read_in_parameters_from_file(input_filename);
    int seed = parameter_list.get_seed();
    ran_gen_ptr = std::shared_ptr<RandomUtil::Random>(
                                            new RandomUtil::Random(seed));
    mc_glauber_ptr = std::unique_ptr<Glauber>(
                                new Glauber(parameter_list, ran_gen_ptr));
    std::cout << "Generating " << nev << " events ... " << std::endl;
}


void EventGenerator::generate_events() {
    // this file records all the essential information for the generated events
    std::ofstream record_file("events_summary.dat", std::ios::out);
    record_file << "# event_id  Npart  Ncoll  Nstrings  b(fm)" << std::endl;

    int iev = 0;
    while (iev < nev) {
        mc_glauber_ptr->make_nuclei();
        auto Ncoll = mc_glauber_ptr->make_collision_schedule();
        auto Npart = mc_glauber_ptr->get_Npart();
        auto Nstrings = mc_glauber_ptr->decide_QCD_strings_production();
        if (event_of_interest_trigger(Npart, Ncoll, Nstrings)) {
            iev++;
            Ncoll = mc_glauber_ptr->perform_string_production();

            std::ostringstream filename;
            filename << "strings_event_" << iev << ".dat";
            mc_glauber_ptr->output_QCD_strings(filename.str());
            
            // write event information to the record file
            auto b = mc_glauber_ptr->get_impact_parameter();
            record_file << iev << "  " << Npart << "  " << Ncoll << "  "
                        << Nstrings << "  " << b << std::endl;
        }
    }
    record_file.close();
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
