// Copyright (C) 2018 Chun Shen

#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_

#include <vector>
#include <string>

using std::vector;
using std::string;

namespace StringUtility {

vector<string> parse_a_line(string line, string delimiter);
vector<string> parse_a_line(string line, string delimiter, string ignore);

}


#endif  // SRC_UTIL_H_
