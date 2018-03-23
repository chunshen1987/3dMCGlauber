// Copyright (C) 2018 Chun Shen

#include <vector>
#include <string>
#include "Util.h"
#include "doctest.h"

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

TEST_CASE("Test parse_a_line level 1"){
    auto word_list1 = StringUtility::parse_a_line("parse a line", " ");
    CHECK(word_list1.size() == 3);
    CHECK(word_list1[0] == "parse");
    CHECK(word_list1[1] == "a");
    CHECK(word_list1[2] == "line");

    auto word_list2 = StringUtility::parse_a_line("one<two<three", "<");
    CHECK(word_list2.size() == 3);
    CHECK(word_list2[0] == "one");
    CHECK(word_list2[1] == "two");
    CHECK(word_list2[2] == "three");
    
    auto word_list3 = StringUtility::parse_a_line("a+=bc+=def+=++", "+=");
    CHECK(word_list3.size() == 4);
    CHECK(word_list3[0] == "a");
    CHECK(word_list3[1] == "bc");
    CHECK(word_list3[2] == "def");
    CHECK(word_list3[3] == "++");
}
    
TEST_CASE("Test parse_a_line level 2"){
    auto word_list1 = StringUtility::parse_a_line("one  two three   ", " ");
    CHECK(word_list1.size() == 3);
    CHECK(word_list1[0] == "one");
    CHECK(word_list1[1] == "two");
    CHECK(word_list1[2] == "three");
    
    word_list1 = StringUtility::parse_a_line("  one  two three   ", " ");
    CHECK(word_list1.size() == 3);
    CHECK(word_list1[0] == "one");
    CHECK(word_list1[1] == "two");
    CHECK(word_list1[2] == "three");

    word_list1 = StringUtility::parse_a_line("one  two three   ", "  ");
    CHECK(word_list1.size() == 3);
    CHECK(word_list1[0] == "one");
    CHECK(word_list1[1] == "two three");
    CHECK(word_list1[2] == " ");
    
    auto word_list2 = StringUtility::parse_a_line("one,,two,three,", ",");
    CHECK(word_list2.size() == 3);
    CHECK(word_list2[0] == "one");
    CHECK(word_list2[1] == "two");
    CHECK(word_list2[2] == "three");
    
    word_list2 = StringUtility::parse_a_line("one,,two,three,", ",,");
    CHECK(word_list2.size() == 2);
    CHECK(word_list2[0] == "one");
    CHECK(word_list2[1] == "two,three,");
}

TEST_CASE("Test parse_a_line_with_ignore"){
    auto word_list1 = (
                StringUtility::parse_a_line("one  two #three   ", " ", "#"));
    CHECK(word_list1.size() == 2);
    CHECK(word_list1[0] == "one");
    CHECK(word_list1[1] == "two");
    
    word_list1 = StringUtility::parse_a_line("#one  two #three   ", " ", "#");
    CHECK(word_list1.size() == 0);
    
    word_list1 = StringUtility::parse_a_line("one  two #three   ", " ", "");
    CHECK(word_list1.size() == 3);
    CHECK(word_list1[0] == "one");
    CHECK(word_list1[1] == "two");
    CHECK(word_list1[2] == "#three");
}

