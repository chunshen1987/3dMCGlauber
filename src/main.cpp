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
    return(0);
}
