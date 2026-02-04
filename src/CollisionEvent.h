// Copyright @ Chun Shen 2018

#ifndef SRC_COLLISIONEVENT_H_
#define SRC_COLLISIONEVENT_H_

#include <memory>

#include "Nucleon.h"
#include "data_structs.h"

namespace MCGlb {

class CollisionEvent {
  private:
    SpatialVec x_coll;
    std::weak_ptr<Nucleon> proj_nucleon;
    std::weak_ptr<Nucleon> targ_nucleon;
    int proj_collided_times;
    int targ_collided_times;
    int produced_n_strings_;

  public:
    CollisionEvent() = default;
    CollisionEvent(
        SpatialVec x_in, std::shared_ptr<Nucleon> proj_in,
        std::shared_ptr<Nucleon> targ_in);

    bool operator<(const CollisionEvent &event1) const {
        return (x_coll[0] < event1.get_collision_time());
    }

    SpatialVec get_collision_position() const { return (x_coll); }
    real get_collision_time() const { return (x_coll[0]); }
    std::weak_ptr<Nucleon> get_proj_nucleon_ptr() const {
        return (proj_nucleon);
    }
    std::weak_ptr<Nucleon> get_targ_nucleon_ptr() const {
        return (targ_nucleon);
    }
    int get_proj_collided_times() const { return (proj_collided_times); }
    int get_targ_collided_times() const { return (targ_collided_times); }
    bool is_valid() const;

    void set_produced_n_strings(int n_strings) {
        produced_n_strings_ = n_strings;
    }
    int get_num_strings() const { return (produced_n_strings_); }
};

struct compare_collision_time {
    bool operator()(
        const std::shared_ptr<CollisionEvent> lhs,
        const std::shared_ptr<CollisionEvent> rhs) const {
        return (lhs->get_collision_time() < rhs->get_collision_time());
    }
};

}  // namespace MCGlb

#endif  // SRC_COLLISIONEVENT_H_
