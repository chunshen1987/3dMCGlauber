// Copyright @ Chun Shen 2018

#ifndef SRC_RANDOM_H_
#define SRC_RANDOM_H_

#include <random>
#include <memory>

namespace RandomUtil {

class Random {
 private:
    int seed;
    std::random_device ran_dev;
    std::unique_ptr<std::mt19937> ran_generator;
    std::uniform_real_distribution<double> rand_uniform_dist;
    
 public:
    Random(int seed_in, double min = 0.0, double max = 1.0);
    double rand_uniform() {return(rand_uniform_dist(*ran_generator));}
    int get_seed() const {return(seed);}
};

}

#endif  // SRC_RANDOM_H_
