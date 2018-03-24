// Copyright (C) 2018 Chun Shen
#include "doctest.h"
#include "data_structs.h"
#include "Particle.h"

using MCGlb::Particle;
using MCGlb::real;
using MCGlb::SpatialVec;
using MCGlb::MomentumVec;

TEST_CASE("Test set and get functions") {
    Particle testParticle;
    SpatialVec x = {1.0, 0.0, -2.0, 3.0};
    testParticle.set_x(x);
    CHECK(testParticle.get_x()[0] == 1.0);
    CHECK(testParticle.get_x()    == x);

    MomentumVec p = {3.0, 0.0, -2.0, 1.0};
    testParticle.set_p(p);
    CHECK(testParticle.get_p()[0] == 3.0);
    CHECK(testParticle.get_p()    == p);
    
    p = {0.0};
    testParticle.set_p(p);
    CHECK(testParticle.get_p()[0] == 0.0);

    testParticle.set_mass(1.0);
    CHECK(testParticle.get_mass() == 1.0);
}

TEST_CASE("Test constructors") {
    SpatialVec x = {1.0, 0.0, -2.0, 3.0};
    MomentumVec p = {5.0, 0.0, -4.0, 3.0};
    Particle testParticle1(x, p);
    CHECK(testParticle1.get_x()    == x);
    CHECK(testParticle1.get_p()    == p);
    CHECK(testParticle1.get_mass() == 0.0);
    
    x = {1.0, 0.0, -2.0, 3.0};
    p = {5.0, 0.0, -4.0, 3.0};
    Particle testParticle2(x, p, 0.0);
    CHECK(testParticle2.get_x()    == x);
    CHECK(testParticle2.get_p()    == p);
    CHECK(testParticle2.get_mass() == 0.0);
}

TEST_CASE("Test copy") {
    SpatialVec x = {1.0, 0.0, -2.0, 3.0};
    MomentumVec p = {5.0, 0.0, -4.0, 3.0};
    Particle testParticle1(x, p);
    auto testParticle2 = testParticle1;
    CHECK(testParticle2.get_x()    == x);
    CHECK(testParticle2.get_p()    == p);
    CHECK(testParticle2.get_mass() == 0.0);
}
