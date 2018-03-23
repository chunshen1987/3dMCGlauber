// Copyright @ Chun Shen 2018

#include <iostream>
#include "Glauber.h"
#include "Parameters.h"
#include "ParametersMap.h"
#include "Util.h"

int main() {
    MCGlb::Parameters Partest;
    Partest.set_parameters();
    std::cout << "b = " << Partest.get_b() << std::endl;
    Glauber testGlauber;
    return(0);
}
