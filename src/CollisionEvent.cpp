// Copyright @ Chun Shen 2018

#include "CollisionEvent.h"

using std::shared_ptr;

namespace MCGlb {

CollisionEvent::CollisionEvent(SpatialVec x_coll_in,
                               shared_ptr<Nucleon> proj_in,
                               shared_ptr<Nucleon> targ_in) {
    proj_nucleon = proj_in;
    targ_nucleon = targ_in;
    x_coll = x_coll_in;
    proj_collided_times = proj_nucleon.lock()->get_collided_times();
    targ_collided_times = targ_nucleon.lock()->get_collided_times();
    produced_a_string = false;
}

bool CollisionEvent::is_valid() const {
    bool flag = true;
    if (proj_collided_times != proj_nucleon.lock()->get_collided_times())
        flag = false;
    if (targ_collided_times != targ_nucleon.lock()->get_collided_times())
        flag = false;
    return(flag);
}

}
