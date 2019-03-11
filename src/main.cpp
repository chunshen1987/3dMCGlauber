// Copyright @ Chun Shen 2018

#include "EventGenerator.h"
#include <string>

int main(int argc, char* argv[]) {
    std::string input_filename = "input";
    int nev = 200;
    if (argc > 1) {
        nev = std::stoi(*(argv + 1));
    }
    if (argc > 2) {
        input_filename = *(argv + 2);
    }

    MCGlb::EventGenerator mc_gen(input_filename, nev);
    mc_gen.generate_events();
    return(0);
}
