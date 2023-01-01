// Copyright @ Chun Shen 2018

#ifndef SRC_EVENT_GENERATOR_H_
#define SRC_EVENT_GENERATOR_H_

#include <memory>
#include <string>
#include "Glauber.h"
#include "Parameters.h"
#include "RandomUlty.h"
#include "pretty_ostream.h"

namespace MCGlb {

class EventGenerator {
 private:
    Parameters parameter_list_;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr_;
    std::unique_ptr<Glauber> mc_glauber_ptr_;
    bool statistics_only_;
    pretty_ostream messager;
    int Ncoll_;
    int Npart_;
    std::vector<double> HardPartonPosAndMomProj;
    std::vector<double> HardPartonPosAndMomTarg;
    double proj_t, proj_x, proj_y, proj_z;
    double proj_E, proj_px, proj_py, proj_pz;
    double targ_t, targ_x, targ_y, targ_z;
    double targ_E, targ_px, targ_py, targ_pz;
    real ecm_;
 public:
    EventGenerator() = default;
    EventGenerator(std::string input_filename, int seed=0);
    ~EventGenerator() {};

    void generate_events(int nev, int event_id_offset=0);

    //! get the collisions information for the JETSCAPE framework
    void generate_pre_events();

    void New_Para_pointer(int seed);
    void set_parameter(string par, double val);
    void set_parameter(string par, float val);
    void set_parameter(string par, int val);

    //! after substracted the parton's momentum, 
    //! generate the 3D Glauber initial conditions for MUSIC
    void generate_strings();

    bool event_of_interest_trigger(int Npart, int Ncoll, int Nstrings);

    //! calculate the total nucleon density at Lab frame, unit is 1/fm^3
    double MCGlb_nucleon_density(double t, double x,
                                 double y, double z);
    //! calculate the target/projectile nucleon density at Lab frame,
    //
    //! unit is 1/fm^3
    double MCGlb_target_nucleon_density(double t, double x,
                                        double y, double z);

    double MCGlb_projectile_nucleon_density(double t, double x,
                                            double y, double z);

    void GetMomandPos_Proj(double &t, double &x, double &y, double &z,
                           double &E, double &px, double &py, double &pz);
    void GetMomandPos_Targ(double &t, double &x, double &y, double &z,
                           double &E, double &px, double &py, double &pz);
    void GetMom_Proj(double &E, double &px, double &py, double &pz);
    void GetMom_Targ(double &E, double &px, double &py, double &pz);
    void GetHardPos(double &t, double &x, double &y, double &z);

    std::vector<double> MCGlb_projectile_nucleon_z();
    std::vector<double> MCGlb_target_nucleon_z();
    std::vector<double> GetQuarkPosProj();
    std::vector<double> GetQuarkPosTarg();
    std::vector<double> GetRemMom_Proj();
    std::vector<double> GetRemMom_Targ();

    std::vector<CollisionEvent> get_CollisionEventvector();
};

};

#endif  // SRC_EVENT_GENERATOR_H_
