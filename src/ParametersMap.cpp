// Copyright (C) 2018 Chun Shen
#include "ParametersMap.h"
#include "Util.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <sstream>

using std::string;

ParametersMap::~ParametersMap() {
    parameter_map.clear();
}

void ParametersMap::set_parameter(string par, double val) {
    std::ostringstream param_val;
    param_val << val;
    set_parameter(par, param_val.str());
}

void ParametersMap::set_parameter(string par, int val) {
    std::ostringstream param_val;
    param_val << val;
    set_parameter(par, param_val.str());
}

void ParametersMap::set_parameter(string par, float val) {
    std::ostringstream param_val;
    param_val << val;
    set_parameter(par, param_val.str());
}

void ParametersMap::read_in_parameters_from_file(string filename) {
    std::ifstream input_file(filename);
    if (!input_file.is_open()) {
        std::cout << "[ParametersMap]: can not open the file: " << filename 
                  << std::endl;
        exit(1);
    }
    string line;
    while (getline(input_file,line)) {
        auto param = StringUtility::parse_a_line(line, " ", "#");
        if (param.size() > 0) {
            set_parameter(param[0], param[1]);
        }
    }
}

void ParametersMap::print_parameter_list() const {
    std::cout << "==============================================" << std::endl;
    std::cout << "Input Parameter list:" << std::endl;
    std::cout << "==============================================" << std::endl;
    for (auto& it: parameter_map) {
        std::cout << it.first << " = " << it.second << "\n";
    }
    std::cout << "==============================================" << std::endl;
}
