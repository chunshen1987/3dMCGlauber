// Copyright (C) 2018 Chun Shen
#include "Quark.h"

#include "doctest.h"

using MCGlb::real;

TEST_CASE("Test set and get functions") {
    MCGlb::SpatialVec x = {1.0, 0.0, -2.0, 3.0};
    MCGlb::MomentumVec p = {5.0, 0.0, -4.0, 3.0};
    MCGlb::Quark test_quark(x, p);
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

TEST_CASE("Test new variables in constructors") {
    MCGlb::SpatialVec x = {1.0, 0.0, -2.0, 3.0};
    MCGlb::MomentumVec p = {5.0, 0.0, -4.0, 3.0};
    real mass = 1.234;
    real pdf_x = 0.2;
    real b1 = 1;
    real c1 = 0;
    real s1 = -1;
    real b2 = 0;
    real c2 = 1;
    real s2 = 0;
    real b3 = 1;
    real c3 = -1;
    real s3 = 1;
    MCGlb::Quark testParticle1(x, p, b1, c1, s1);
    MCGlb::Quark testParticle2(x, p, mass, b2, c2, s2);
    MCGlb::Quark testParticle3(x, pdf_x, b3, c3, s3);
    CHECK(testParticle1.get_baryon() == b1);
    CHECK(testParticle1.get_charge() == c1);
    CHECK(testParticle1.get_strange() == s1);
    CHECK(testParticle2.get_baryon() == b2);
    CHECK(testParticle2.get_charge() == c2);
    CHECK(testParticle2.get_strange() == s2);
    CHECK(testParticle3.get_baryon() == b3);
    CHECK(testParticle3.get_charge() == c3);
    CHECK(testParticle3.get_strange() == s3);
}

TEST_CASE("Test new variables setters") {
    MCGlb::SpatialVec x = {1.0, 0.0, -2.0, 3.0};
    MCGlb::MomentumVec p = {5.0, 0.0, -4.0, 3.0};
    real b1 = 1;
    real c1 = 0;
    real s1 = -1;
    real b2 = 0;
    real c2 = 1;
    real s2 = 0;
    MCGlb::Quark testParticle1(x, p, b1, c1, s1);
    testParticle1.set_baryon(b2);
    testParticle1.set_charge(c2);
    testParticle1.set_strange(s2);

    MCGlb::Quark testParticle2(x, p, b2, c2, s2);
    testParticle2.set_charges(b1, c1, s1);

    CHECK(testParticle1.get_baryon() == b2);
    CHECK(testParticle1.get_charge() == c2);
    CHECK(testParticle1.get_strange() == s2);

    CHECK(testParticle2.get_baryon() == b1);
    CHECK(testParticle2.get_charge() == c1);
    CHECK(testParticle2.get_strange() == s1);
}
