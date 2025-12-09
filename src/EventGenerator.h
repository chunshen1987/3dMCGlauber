// Copyright @ Chun Shen 2018

#ifndef SRC_EVENT_GENERATOR_H_
#define SRC_EVENT_GENERATOR_H_

#include <memory>
#include <string>
#include <vector>
#include "Glauber.h"
#include "MakeDensity.h"
#include "Parameters.h"
#include "RandomUlty.h"
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

    int Ncoll_;
    int Npart_;
    std::vector<double> proj_t_;
    std::vector<double> proj_x_;
    std::vector<double> proj_y_;
    std::vector<double> proj_z_;
    std::vector<double> proj_E_;
    std::vector<double> proj_px_;
    std::vector<double> proj_py_;
    std::vector<double> proj_pz_;
    std::vector<double> targ_t_;
    std::vector<double> targ_x_;
    std::vector<double> targ_y_;
    std::vector<double> targ_z_;
    std::vector<double> targ_E_;
    std::vector<double> targ_px_;
    std::vector<double> targ_py_;
    std::vector<double> targ_pz_;

    real ecm_;
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

    //! get the collisions information for the JETSCAPE framework
    void generate_pre_events();

    void set_parameter(string par, double val);
    void set_parameter(string par, float val);
    void set_parameter(string par, int val);
    void set_parameter(string par, string val);
    void print_parameter_list() const {
        parameter_list_.print_parameter_list();
    };

    void initializeGlauberModel() {
        mc_glauber_ptr_ = std::unique_ptr<Glauber>(
                    new Glauber(parameter_list_, ran_gen_ptr_));
    };

    //! after substracted the parton's momentum,
    //! generate the 3D Glauber initial conditions for MUSIC
    std::vector< std::vector<real> > generate_strings(bool wound_nucleons,
                                                      int event_id);

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

    void GetMomandPos_Proj(
        std::vector<double> &t, std::vector<double> &x,
        std::vector<double> &y, std::vector<double> &z,
        std::vector<double> &E, std::vector<double> &px,
        std::vector<double> &py, std::vector<double> &pz);
    void GetMomandPos_Targ(
        std::vector<double> &t, std::vector<double> &x,
        std::vector<double> &y, std::vector<double> &z,
        std::vector<double> &E, std::vector<double> &px,
        std::vector<double> &py, std::vector<double> &pz);

    std::vector<double> MCGlb_projectile_nucleon_z();
    std::vector<double> MCGlb_target_nucleon_z();
    std::vector<SpatialVec> MCGlb_projectile_nucleon_xyz();
    std::vector<SpatialVec> MCGlb_target_nucleon_xyz();
    std::vector<double> MCGlb_Proj_hot_spot_x();
    std::vector<double> MCGlb_Targ_hot_spot_x();
    std::vector<double> GetQuarkPosProj(double t, double x, double y, double z);
    std::vector<double> GetQuarkPosTarg(double t, double x, double y, double z);
    std::vector<double> GetRemMom_Proj();
    std::vector<double> GetRemMom_Targ();
    
    std::vector<int> Get_Hard_p_id();

    std::vector<CollisionEvent> get_CollisionEventvector();
    bool event_of_interest_trigger(const int Npart, const int Ncoll,
                                   const int Nstrings) const;
};

};

#endif  // SRC_EVENT_GENERATOR_H_
