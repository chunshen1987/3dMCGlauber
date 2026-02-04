// Copyright (C) 2018 Chun Shen

#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_

#include <string>
#include <vector>

using std::string;
using std::vector;

namespace StringUtility {

vector<string> parse_a_line(string line, string delimiter);
vector<string> parse_a_line(string line, string delimiter, string ignore);

}  // namespace StringUtility

#endif  // SRC_UTIL_H_
