// Copyright @ Chun Shen 2018

#include "EventGenerator.h"
#include <string>
#include <iostream>

void print_help() {
    std::cout << "Usage: ./3dMCGlb.e nev [input_file] [seed]" << std::endl;
    std::cout << "Default input_file: input" << std::endl;
    std::cout << "Default seed: -1 (use device random seed)" << std::endl;
}

int main(int argc, char* argv[]) {
    std::string input_filename = "input";
    int nev = 100;
    long long int seed = 0;
    if (argc == 1) {
        print_help();
        exit(1);
    }

    if (argc > 1) {
        nev = std::stoi(*(argv + 1));
    }
    if (argc > 2) {
        input_filename = *(argv + 2);
    }
    if (argc > 3) {
        seed = std::stoll(*(argv + 3));
    }

    MCGlb::EventGenerator mc_gen(input_filename, seed);
    mc_gen.generateMinBiasEventList();
    mc_gen.generate_events(nev);
    return(0);
}
