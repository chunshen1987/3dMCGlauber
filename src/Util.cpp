// Copyright (C) 2018 Chun Shen

#include <vector>
#include <string>
#include "Util.h"

namespace StringUtility {

vector<string> parse_a_line(string line, string delimiter, string ignore) {
    string trimed_line = line;
    if (ignore.length() != 0) {
        trimed_line = line.substr(0, line.find(ignore));
    }
    return(parse_a_line(trimed_line, delimiter));
}

vector<string> parse_a_line(string line, string delimiter) {
    string word;
    vector<string> word_list;
    size_t pos = 0;
    while ((pos = line.find(delimiter)) != string::npos) {
        word = line.substr(0, pos);
        if (word.length() > 0) {
            word_list.push_back(word);
        }
        line.erase(0, pos + delimiter.length());
    }
    if (line.length() != 0) {
        word_list.push_back(line);
    }
    return(word_list);
}

}
