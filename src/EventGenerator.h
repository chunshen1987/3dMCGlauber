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
    int nev;
    Parameters parameter_list;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr;
    std::unique_ptr<Glauber> mc_glauber_ptr;
    bool statistics_only;
    pretty_ostream messager;

 public:
    EventGenerator() = default;
    EventGenerator(std::string input_filename, int nev);
    ~EventGenerator() {};

    void generate_events();
    bool event_of_interest_trigger(int Npart, int Ncoll, int Nstrings);
};

};

#endif  // SRC_EVENT_GENERATOR_H_
