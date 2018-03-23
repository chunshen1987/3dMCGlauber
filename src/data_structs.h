// Copyright @ Chun Shen 2018

#ifndef SRC_DATA_STRUCTS_H_
#define SRC_DATA_STRUCTS_H_

#include <array>
#include <vector>

namespace MCGlb {

typedef double real;
typedef std::array<real, 4> SpatialVec;
typedef std::array<real, 4> MomentumVec;

struct Quark {
    SpatialVec x;
    MomentumVec p;
    real mass;
};

struct Nucleon {
    SpatialVec x;
    MomentumVec p;
    real mass;
    std::vector<Quark> quarkList;
};

}


#endif   // SRC_DATA_STRUCTS_H_
