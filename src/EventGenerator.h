// Copyright @ Chun Shen 2018

#ifndef SRC_EVENT_GENERATOR_H_
#define SRC_EVENT_GENERATOR_H_

#include <memory>
#include <string>
#include "Glauber.h"
#include "Parameters.h"
#include "Random.h"
#include "pretty_ostream.h"

namespace MCGlb {

class EventGenerator {
 private:
    Parameters parameter_list_;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr_;
    std::unique_ptr<Glauber> mc_glauber_ptr_;
    bool statistics_only_;
    pretty_ostream messager;

 public:
    EventGenerator() = default;
    EventGenerator(std::string input_filename, int seed=0);
    ~EventGenerator() {};

    void generate_events(int nev, int event_id_offset=0);
    bool event_of_interest_trigger(int Npart, int Ncoll, int Nstrings);
};

};

#endif  // SRC_EVENT_GENERATOR_H_
