// Copyright @ Chun Shen 2018

#ifndef SRC_RANDOM_H_
#define SRC_RANDOM_H_

#include <memory>
#include <random>

namespace RandomUtil {

class Random {
  private:
    unsigned int seed_;
    double Gamma_beta_;
    std::shared_ptr<std::mt19937> ran_generator_;
    std::uniform_real_distribution<double> rand_uniform_dist_;
    std::uniform_int_distribution<int> rand_int_uniform_dist_;
    std::normal_distribution<double> rand_normal_dist_;
    std::gamma_distribution<double> rand_gamma_dist_;

  public:
    Random(int seed_in, double min = 0.0, double max = 1.0);
    Random(int seed_in, int min = 0, int max = 1);
    Random(int seed_in, double min, double max, double Gamma_beta);
    double rand_uniform() { return (rand_uniform_dist_(*ran_generator_)); }
    double rand_normal(double mu, double sigma) {
        return (mu + sigma * rand_normal_dist_(*ran_generator_));
    }

    double rand_gamma_dis() { return (rand_gamma_dist_(*ran_generator_)); }

    int rand_int_uniform() { return (rand_int_uniform_dist_(*ran_generator_)); }
    unsigned int get_seed() const { return (seed_); }
    double get_Gamma_beta() const { return (Gamma_beta_); }
    std::shared_ptr<std::mt19937> getRanGenerator() { return (ran_generator_); }
};

}  // namespace RandomUtil

#endif  // SRC_RANDOM_H_
