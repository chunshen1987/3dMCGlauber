// Copyright (C) 2018 Chun Shen
#ifndef SRC_PARAMETERSMAP_H_
#define SRC_PARAMETERSMAP_H_

#include <string>
#include <map>

using std::string;
using std::map;

class ParametersMap {
 private:
     map<string, string> parameter_map;

 public:
     ParametersMap() = default;
     ~ParametersMap();

     void set_parameter(string par, string val) {parameter_map[par] = val;}
     string get_param_val(string par) const {return(parameter_map.at(par));}

     void read_in_parameters_from_file(string filename);
     int get_parameter_list_size() const {return(parameter_map.size());}
     void print_parameter_list() const;
};

#endif  // SRC_PARAMETERSMAP_H_
