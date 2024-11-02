// Copyright (C) 2018 Chun Shen
#include "Util.h"

#include "doctest.h"

TEST_CASE("Test parse_a_line level 1") {
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

TEST_CASE("Test parse_a_line level 2") {
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

TEST_CASE("Test parse_a_line_with_ignore") {
    auto word_list1 =
        (StringUtility::parse_a_line("one  two #three   ", " ", "#"));
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

    word_list1 = StringUtility::parse_a_line("      #three   ", " ", "");
    CHECK(word_list1.size() == 1);
    CHECK(word_list1[0] == "#three");

    word_list1 = StringUtility::parse_a_line("      #three   ", " ", "#");
    CHECK(word_list1.size() == 0);
}
