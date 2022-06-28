// Copyright (C) 2018 Chun Shen
#ifndef SRC_PARAMETERSMAP_H_
#define SRC_PARAMETERSMAP_H_

#include <string>
#include <map>

using std::string;
using std::map;
using std::stod;
using std::stoi;

class ParametersMap {
 private:
    map<string, string> parameter_map;

 public:
    ParametersMap() = default;
    ~ParametersMap();

    void set_parameter(string par, string val) {parameter_map[par] = val;}
    void set_parameter(string par, double val);
    void set_parameter(string par, float val);
    void set_parameter(string par, int val);

    string get_param_val(string par) const {return(parameter_map.at(par));}
    int get_param_int(string par, int defaultVal=0) const {
        if (checkParamIsDefined(par)) {
            return(stoi(get_param_val(par)));
        } else {
            return(defaultVal);
        }
    }

    double get_param_double(string par, double defaultVal=0.) const {
        if (checkParamIsDefined(par)) {
            return(stod(get_param_val(par)));
        } else {
            return(defaultVal);
        }
    }

    bool checkParamIsDefined(string paramName) const;
    void read_in_parameters_from_file(string filename);
    int get_parameter_list_size() const {return(parameter_map.size());}
    void print_parameter_list() const;
};

#endif  // SRC_PARAMETERSMAP_H_
