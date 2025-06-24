// Copyright (C) 2018 Chun Shen
#include "Quark.h"

#include "doctest.h"

using MCGlb::real;

TEST_CASE("Test set and get functions") {
    MCGlb::Quark test_quark;
    test_quark.set_pdf_x(0.2);
    CHECK(test_quark.get_pdf_x() == 0.2);
}

TEST_CASE("Test constructors") {
    MCGlb::SpatialVec x = {1.0, 0.0, -2.0, 3.0};
    MCGlb::MomentumVec p = {5.0, 0.0, -4.0, 3.0};
    MCGlb::Quark testParticle1(x, p);
    CHECK(testParticle1.get_x() == x);
    CHECK(testParticle1.get_p() == p);
    CHECK(testParticle1.get_mass() == 0.0);
}

TEST_CASE("Test copy") {
    MCGlb::SpatialVec x = {1.0, 0.0, -2.0, 3.0};
    MCGlb::MomentumVec p = {5.0, 0.0, -4.0, 3.0};
    MCGlb::Quark testParticle1(x, p);
    auto testParticle2 = testParticle1;
    CHECK(testParticle2.get_x() == x);
    CHECK(testParticle2.get_p() == p);
    CHECK(testParticle2.get_mass() == 0.0);
}

TEST_CASE("Test new variables") {
    MCGlb::SpatialVec x = {1.0, 0.0, -2.0, 3.0};
    MCGlb::MomentumVec p = {5.0, 0.0, -4.0, 3.0};
    real b1 = 1;
    real c1 = 0;
    real s1 = -1;
    real b2 = 0;
    real c2 = 1;
    real s2 = 0;
    MCGlb::Quark testParticle1(x, p, b1, c1, s1);
    MCGlb::Quark testParticle2(x, p, b2, c2, s2);
    CHECK(testParticle1.get_baryon() == b1);
    CHECK(testParticle1.get_charge() == c1);
    CHECK(testParticle1.get_strange() == s1);
    CHECK(testParticle2.get_baryon() == b2);
    CHECK(testParticle2.get_charge() == c2);
    CHECK(testParticle2.get_strange() == s2);
}
