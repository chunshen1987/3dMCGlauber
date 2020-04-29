// Copyright @ Chun Shen 2018

#include "Random.h"

namespace RandomUtil {

Random::Random(int seed, double min, double max) :
    rand_uniform_dist_(min, max) {
    seed_ = seed;
    if (seed == -1) {
        std::random_device ran_dev;
        seed_ = ran_dev();
    }
    ran_generator_ = std::unique_ptr<std::mt19937>(new std::mt19937(seed_));
}

}

