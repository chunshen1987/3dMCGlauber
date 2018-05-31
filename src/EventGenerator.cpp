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
}

void EventGenerator::generate_events() {
    std::ofstream of("events_summary.dat", std::ios::out);
    of << "# event_id  Npart  Ncoll  Nstrings" << std::endl;
    int iev = 0;
    while (iev < nev) {
        mc_glauber_ptr->make_nuclei();
        auto Ncoll = mc_glauber_ptr->make_collision_schedule();
        auto Npart = mc_glauber_ptr->get_Npart();
        if (Npart > 1) {
            iev++;
            auto Nstrings = mc_glauber_ptr->decide_QCD_strings_production();
            Ncoll = mc_glauber_ptr->perform_string_production();
            std::ostringstream filename;
            filename << "strings_event_" << iev << ".dat";
            mc_glauber_ptr->output_QCD_strings(filename.str());
            of << iev << "  " << Npart << "  " << Ncoll << "  "
               << Nstrings << std::endl;
        }
    }
    of.close();
}

};
