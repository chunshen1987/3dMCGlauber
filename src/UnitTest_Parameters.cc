// Copyright (C) 2018 Chun Shen
#include "doctest.h"
#include "Parameters.h"

TEST_CASE("Test set and get parameters") {
    MCGlb::Parameters parameter_list;
    parameter_list.set_b(2.0);
    CHECK(parameter_list.get_b() == 2.0);
}
