// Copyright (C) 2018 Chun Shen
#include "doctest.h"
#include "CollisionEvent.h"
#include "RandomUlty.h"
#include <set>
#include <vector>
#include <algorithm>

using MCGlb::SpatialVec;
using MCGlb::MomentumVec;
using MCGlb::Nucleon;
using MCGlb::CollisionEvent;
using std::shared_ptr;

TEST_CASE("Test constructor") {
    SpatialVec x_coll = {4.0, 1.0, 2.0, 3.0};
    SpatialVec x1 = {0.0, 0.5, 2.0, 3.0};
    SpatialVec x2 = {0.0, 1.5, 2.0, 3.0};
    MomentumVec p1 = {0.0};
    MomentumVec p2 = {0.0};
    std::shared_ptr<MCGlb::RandomUtil::Random> ran_gen_ptr(
                                    new MCGlb::RandomUtil::Random(0, 0., 1.0));
    shared_ptr<Nucleon> proj(new Nucleon(x1, p1, ran_gen_ptr));
    shared_ptr<Nucleon> targ(new Nucleon(x2, p2, ran_gen_ptr));
    CollisionEvent test_event(x_coll, proj, targ);
    CHECK(test_event.get_proj_collided_times() == 0);
    CHECK(test_event.get_targ_collided_times() == 0);
    CHECK(test_event.get_collision_time() == 4.0);
}

TEST_CASE("Test the copy") {
    SpatialVec x_coll = {4.0, 1.0, 2.0, 3.0};
    SpatialVec x1 = {0.0, 0.5, 2.0, 3.0};
    SpatialVec x2 = {0.0, 1.5, 2.0, 3.0};
    MomentumVec p1 = {0.0};
    MomentumVec p2 = {0.0};
    std::shared_ptr<MCGlb::RandomUtil::Random> ran_gen_ptr(
                                    new MCGlb::RandomUtil::Random(0, 0., 1.0));
    shared_ptr<Nucleon> proj(new Nucleon(x1, p1, ran_gen_ptr));
    shared_ptr<Nucleon> targ(new Nucleon(x2, p2, ran_gen_ptr));
    CollisionEvent test_event(x_coll, proj, targ);
    auto test_event2 = test_event;
    CHECK(test_event2.get_collision_time() == 4.0);
}

TEST_CASE("Test the collision validation function") {
    SpatialVec x_coll = {4.0, 1.0, 2.0, 3.0};
    SpatialVec x1 = {0.0, 0.5, 2.0, 3.0};
    SpatialVec x2 = {0.0, 1.5, 2.0, 3.0};
    MomentumVec p1 = {0.0};
    MomentumVec p2 = {0.0};
    std::shared_ptr<MCGlb::RandomUtil::Random> ran_gen_ptr(
                                    new MCGlb::RandomUtil::Random(0, 0., 1.0));
    shared_ptr<Nucleon> proj(new Nucleon(x1, p1, ran_gen_ptr));
    shared_ptr<Nucleon> targ(new Nucleon(x2, p2, ran_gen_ptr));
    {
        CollisionEvent test_event(x_coll, proj, targ);
        CHECK(test_event.is_valid() == true);
        targ->increment_collided_times();
        CHECK(test_event.is_valid() == false);
    }
    // make sure that proj and targ are still in the momery after the
    // collsion event is deleted
    targ->increment_collided_times();
    CHECK(targ->get_collided_times() == 2);
}

TEST_CASE("Test pushing collision events to a set with defined comparison function") {
    SpatialVec x_coll1 = {4.0, 1.0, 2.0, 3.0};
    SpatialVec x_coll2 = {2.0, 1.0, 2.0, 3.0};
    SpatialVec x1 = {0.0, 0.5, 2.0, 3.0};
    SpatialVec x2 = {0.0, 1.5, 2.0, 3.0};
    MomentumVec p1 = {0.0};
    MomentumVec p2 = {0.0};
    std::shared_ptr<MCGlb::RandomUtil::Random> ran_gen_ptr(
                                    new MCGlb::RandomUtil::Random(0, 0., 1.0));
    shared_ptr<Nucleon> proj(new Nucleon(x1, p1, ran_gen_ptr));
    shared_ptr<Nucleon> targ(new Nucleon(x2, p2, ran_gen_ptr));
    CollisionEvent test_event1(x_coll1, proj, targ);
    CollisionEvent test_event2(x_coll2, proj, targ);
    std::set<CollisionEvent> collision_schedule;
    collision_schedule.insert(test_event1);
    collision_schedule.insert(test_event2);
    CHECK(collision_schedule.size() == 2);
    CHECK((*collision_schedule.begin()).get_collision_time() == 2.0);
}

TEST_CASE("Test pushing collision events to a vector and sort") {
    SpatialVec x_coll1 = {4.0, 1.0, 2.0, 3.0};
    SpatialVec x_coll2 = {2.0, 1.0, 2.0, 3.0};
    SpatialVec x1 = {0.0, 0.5, 2.0, 3.0};
    SpatialVec x2 = {0.0, 1.5, 2.0, 3.0};
    MomentumVec p1 = {0.0};
    MomentumVec p2 = {0.0};
    std::shared_ptr<MCGlb::RandomUtil::Random> ran_gen_ptr(
                                    new MCGlb::RandomUtil::Random(0, 0., 1.0));
    shared_ptr<Nucleon> proj(new Nucleon(x1, p1, ran_gen_ptr));
    shared_ptr<Nucleon> targ(new Nucleon(x2, p2, ran_gen_ptr));
    CollisionEvent test_event1(x_coll1, proj, targ);
    CollisionEvent test_event2(x_coll2, proj, targ);
    std::vector<CollisionEvent> collision_schedule;
    collision_schedule.push_back(test_event1);
    collision_schedule.push_back(test_event2);
    CHECK(collision_schedule.size() == 2);
    std::sort(collision_schedule.begin(), collision_schedule.end());
    CHECK((*collision_schedule.begin()).get_collision_time() == 2.0);
}

TEST_CASE("Test get collision nucleon pointer") {
    SpatialVec x_coll1 = {4.0, 1.0, 2.0, 3.0};
    SpatialVec x1 = {0.0, 0.5, 2.0, 3.0};
    SpatialVec x2 = {0.0, 1.5, 2.0, 3.0};
    MomentumVec p1 = {0.0};
    MomentumVec p2 = {0.0};
    std::shared_ptr<MCGlb::RandomUtil::Random> ran_gen_ptr(
                                    new MCGlb::RandomUtil::Random(0, 0., 1.0));
    shared_ptr<Nucleon> proj(new Nucleon(x1, p1, ran_gen_ptr));
    shared_ptr<Nucleon> targ(new Nucleon(x2, p2, ran_gen_ptr));
    CollisionEvent test_event1(x_coll1, proj, targ);
    auto proj_ptr = test_event1.get_proj_nucleon_ptr().lock();
    CHECK(proj_ptr->get_x() == x1);
    proj_ptr->set_x(x2);
    CHECK(proj_ptr->get_x() == x2);
}
