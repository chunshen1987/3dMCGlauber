// Copyright @ Chun Shen 2018

#ifndef SRC_EVENT_GENERATOR_H_
#define SRC_EVENT_GENERATOR_H_

#include <memory>
#include <string>
#include <vector>
#include "Glauber.h"
#include "MakeDensity.h"
#include "Parameters.h"
#include "Random.h"
#include "pretty_ostream.h"

namespace MCGlb {

class EventGenerator {
 private:
    Parameters parameter_list_;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr_;
    std::unique_ptr<Glauber> mc_glauber_ptr_;
    std::unique_ptr<MakeDensity> density_maker_ptr_;
    bool statistics_only_;
    bool batchDensityOutput_;
    bool batchDensity2DOutput_;
    bool batchEccOutput_;
    bool initialEstOutput_;
    pretty_ostream messager;
    std::vector<float> cenEstMinBiasList_;
    float cenEstMin_;
    float cenEstMax_;

 public:
    EventGenerator() = default;
    EventGenerator(std::string input_filename, int argc, char* argv[],
                   int seed=0);
    ~EventGenerator() {};

    float computeCenEstimator(const int Npart, const int Ncoll,
                              const int Nstrings) const;
    void generateMinBiasEventList();
    void generate_events(int nev, int event_id_offset=0);
    bool event_of_interest_trigger(const int Npart, const int Ncoll,
                                   const int Nstrings) const;
};

};

#endif  // SRC_EVENT_GENERATOR_H_
