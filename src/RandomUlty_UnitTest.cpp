// Copyright (C) 2018 Chun Shen
#include "doctest.h"
#include "RandomUlty.h"

TEST_CASE("Test fixed random seed") {
    MCGlb::RandomUtil::Random ran_test1(1, 0., 1.);
    MCGlb::RandomUtil::Random ran_test2(1, 0., 1.);
    CHECK(ran_test1.get_seed() == 1);
    for (int i = 0; i < 10; i++)
        CHECK(ran_test1.rand_uniform() == ran_test2.rand_uniform());
}

TEST_CASE("Test device random seed") {
    MCGlb::RandomUtil::Random ran_test1(-1, 0., 1.);
    MCGlb::RandomUtil::Random ran_test2(-1, 0., 1.);
    CHECK(ran_test1.get_seed() != ran_test2.get_seed());
    for (int i = 0; i < 10; i++)
        CHECK(ran_test1.rand_uniform() != ran_test2.rand_uniform());
}

TEST_CASE("Test random number range") {
    MCGlb::RandomUtil::Random ran_test1(-1, 1., 2.);
    for (int i = 0; i < 10; i++)
        CHECK(ran_test1.rand_uniform() >= 1.);
}
