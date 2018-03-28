// Copyright @ Chun Shen 2018

#include <iostream>
#include <string>
#include "Glauber.h"
#include "Parameters.h"

int main(int argc, char* argv[]) {
    std::string input_filename = "input";
    if (argc > 1) {
        input_filename = *(argv + 1);
    }
    MCGlb::Parameters parameter_list;
    parameter_list.read_in_parameters_from_file(input_filename);
    MCGlb::Glauber testGlauber(parameter_list);
    for (int iev = 0; iev < 100; iev++) {
        testGlauber.make_nuclei();
        auto Ncoll = testGlauber.make_collision_schedule();
        auto Npart = testGlauber.get_Npart();
        std::cout << "Npart = " << Npart << ", Ncoll = " << Ncoll << std::endl;
    }
    return(0);
}
