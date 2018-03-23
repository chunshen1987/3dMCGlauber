// Copyright @ Chun Shen 2018

#include <iostream>
#include "Glauber.h"
#include "Parameters.h"
#include "Util.h"

int main() {
    auto test = StringUtility::parse_a_line("test a line ", " ");
    for (int i = 0; i < test.size(); i++) {
        std::cout << test[i] << std::endl;
    }
    MCGlb::Parameters Partest;
    Partest.set_parameters();
    std::cout << "b = " << Partest.get_b() << std::endl;
    Glauber testGlauber;
    return(0);
}
