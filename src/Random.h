// Copyright @ Chun Shen 2018

#ifndef SRC_RANDOM_H_
#define SRC_RANDOM_H_

#include <random>
#include <memory>

namespace RandomUtil {

class Random {
 private:
    int seed_;
    std::unique_ptr<std::mt19937> ran_generator_;
    std::uniform_real_distribution<double> rand_uniform_dist_;
    std::uniform_int_distribution<int> rand_int_uniform_dist_;

 public:
    Random(int seed_in, double min = 0.0, double max = 1.0);
    Random(int seed_in, int min = 0, int max = 1);
    double rand_uniform() {return(rand_uniform_dist_(*ran_generator_));}
    int rand_int_uniform() {return(rand_int_uniform_dist_(*ran_generator_));}
    int get_seed() const {return(seed_);}
};

}

#endif  // SRC_RANDOM_H_
