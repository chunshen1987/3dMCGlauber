// Copyright @ Chun Shen 2018

#include <iostream>
#include <string>
#include <memory>
#include "Glauber.h"
#include "Parameters.h"
#include "Random.h"

int main(int argc, char* argv[]) {
    std::string input_filename = "input";
    if (argc > 1) {
        input_filename = *(argv + 1);
    }
    MCGlb::Parameters parameter_list;
    parameter_list.read_in_parameters_from_file(input_filename);
    int seed = parameter_list.get_seed();
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr(
                                            new RandomUtil::Random(seed));
    MCGlb::Glauber testGlauber(parameter_list, ran_gen_ptr);
    int iev = 0;
    while (iev < 100) {
        testGlauber.make_nuclei();
        auto Ncoll = testGlauber.make_collision_schedule();
        auto Npart = testGlauber.get_Npart();
        if (Npart > 1) {
            iev++;
            auto Nstrings = testGlauber.decide_QCD_strings_production();
            std::cout << "Npart = " << Npart << ", Ncoll = " << Ncoll
                      << ", Nstring = " << Nstrings
                      << std::endl;
            std::cout << "averaged connected rate: " << Nstrings/(Npart/2.)
                      << std::endl;
            auto Ncollided = testGlauber.perform_string_production();
            std::cout << "Ncollisions = " << Ncollided << std::endl;
        }
    }
    return(0);
}
