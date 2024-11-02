// Copyright (C) 2018 Chun Shen
#include "Quark.h"

#include "doctest.h"

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
