// Copyright @ Chun Shen 2018

#include "RandomUlty.h"

namespace MCGlb {

namespace RandomUtil {

Random::Random(int seed, double min, double max) :
    rand_uniform_dist_(min, max), rand_normal_dist_(0., 1.),
    rand_gamma_dist_(1.0, 1.0) {
    seed_ = seed;
    Gamma_beta_ = 1.0;
    if (seed == -1) {
        std::random_device ran_dev;
        seed_ = ran_dev();
    }
    ran_generator_ = std::unique_ptr<std::mt19937>(new std::mt19937(seed_));
}


Random::Random(int seed, int min, int max) :
    rand_int_uniform_dist_(min, max) {
    seed_ = seed;
    if (seed == -1) {
        std::random_device ran_dev;
        seed_ = ran_dev();
    }
    ran_generator_ = std::unique_ptr<std::mt19937>(new std::mt19937(seed_));
}


Random::Random(int seed, double min, double max, double Gamma_beta) :
    rand_uniform_dist_(min, max), rand_normal_dist_(0., 1.),
    rand_gamma_dist_(Gamma_beta, Gamma_beta) {
    seed_ = seed;
    Gamma_beta_ = Gamma_beta;
    if (seed == -1) {
        std::random_device ran_dev;
        seed_ = ran_dev();
    }
    ran_generator_ = std::unique_ptr<std::mt19937>(new std::mt19937(seed_));
}

}

}
