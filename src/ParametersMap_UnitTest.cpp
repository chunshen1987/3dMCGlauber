// Copyright (C) 2018 Chun Shen
#include "doctest.h"
#include "ParametersMap.h"

TEST_CASE("Test set and get parameters in ParametersMap") {
    ParametersMap parameter_list;
    parameter_list.set_parameter("test", "1.0");
    CHECK(parameter_list.get_param_val("test") == "1.0");
}

TEST_CASE("Test read in parameters from a file") {
    ParametersMap parameter_list;
    parameter_list.read_in_parameters_from_file("test/test_input");
    CHECK(parameter_list.get_parameter_list_size() == 31);
    CHECK(parameter_list.get_param_val("Projectile") == "Au");
    CHECK(parameter_list.get_param_val("Target") == "Pb");
    CHECK(parameter_list.get_param_val("string_tension") == "1.0");
}
