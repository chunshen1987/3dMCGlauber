// Copyright @ Chun Shen 2018

#include "EventGenerator.h"
#include <string>

int main(int argc, char* argv[]) {
    std::string input_filename = "input";
    int nev = 256;
    int seed_add = 0;
    if (argc > 1) {
        nev = std::stoi(*(argv + 1));
    }
    if (argc > 2) {
        input_filename = *(argv + 2);
    }
    if (argc > 3) {
        seed_add = std::stoi(*(argv + 3));
    }

    MCGlb::EventGenerator mc_gen(input_filename, seed_add);
    mc_gen.generate_events(nev);
    return(0);
}
