// Copyright (C) 2018 Chun Shen
#include "doctest.h"
#include "Parameters.h"

TEST_CASE("Test set and get parameters") {
    MCGlb::Parameters parameter_list;
    parameter_list.set_b(2.0);
    CHECK(parameter_list.get_b() == 2.0);
}

TEST_CASE("Test read in parameters from a file and get functions") {
    MCGlb::Parameters parameter_list;
    parameter_list.read_in_parameters_from_file("test/test_input");
    CHECK(parameter_list.get_parameter_list_size() == 31);
    CHECK(parameter_list.get_b() == 1);
    CHECK(parameter_list.get_use_quarks() == 2);
    CHECK(parameter_list.get_use_energy_dependent_cross_section() == 1);
    CHECK(parameter_list.get_roots() == 19);
    CHECK(parameter_list.get_time_for_seed() == 1);
    CHECK(parameter_list.get_projectle_nucleus_name() == "Au");
    CHECK(parameter_list.get_target_nucleus_name() == "Pb");
}
