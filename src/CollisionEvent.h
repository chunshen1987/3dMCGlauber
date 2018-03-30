// Copyright @ Chun Shen 2018

#ifndef SRC_COLLISIONEVENT_H_
#define SRC_COLLISIONEVENT_H_

#include "data_structs.h"
#include "Nucleon.h"
#include <memory>

namespace MCGlb {

class CollisionEvent {
 private:
    SpatialVec x_coll;
    std::weak_ptr<Nucleon> proj_nucleon;
    std::weak_ptr<Nucleon> targ_nucleon;
    int proj_collided_times;
    int targ_collided_times;
    bool produced_a_string;

 public:
    CollisionEvent() = default;
    CollisionEvent(SpatialVec x_in, std::shared_ptr<Nucleon> proj_in,
                   std::shared_ptr<Nucleon> targ_in);
    
    bool operator< (const CollisionEvent &event1) const {
        return(x_coll[0] < event1.get_collision_time());
    }

    SpatialVec get_collision_position() const {return(x_coll);}
    real get_collision_time() const {return(x_coll[0]);}
    std::weak_ptr<Nucleon> get_proj_nucleon_ptr() const {return(proj_nucleon);}
    std::weak_ptr<Nucleon> get_targ_nucleon_ptr() const {return(targ_nucleon);}
    int get_proj_collided_times() const {return(proj_collided_times);}
    int get_targ_collided_times() const {return(targ_collided_times);}
    bool is_valid() const;

    void set_produced_a_string(bool flag) {produced_a_string = flag;}
    bool is_produced_a_string() const {return(produced_a_string);}
};

struct compare_collision_time {
    bool operator() (const std::shared_ptr<CollisionEvent> lhs,
                     const std::shared_ptr<CollisionEvent> rhs) const {
        return(lhs->get_collision_time() < rhs->get_collision_time());
    }
};

}

#endif  // SRC_COLLISIONEVENT_H_
