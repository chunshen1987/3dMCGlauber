// Copyright @ Chun Shen 2018

#ifndef SRC_NUCLEUS_H_
#define SRC_NUCLEUS_H_

#include "data_structs.h"
#include "Nucleon.h"
#include "Random.h"
#include <vector>
#include <string>

namespace MCGlb {

class Nucleus {
 private:
    std::string name;
    int density_function_type;
    int A;
    int Z;
    WoodsSaxonParam WS_param_vec;       // rho, w, R, a
    real d_min;                         // minimum distance between nucleons

    std::vector<Nucleon> nucleon_list;
    
    std::unique_ptr<RandomUtil::Random> ran_gen_ptr;

 public:
    Nucleus() = default;
    Nucleus(std::string nucleus_name, int seed_in=-1, real d_min_in=0.9);
    ~Nucleus();

    void set_random_seed(int seed_in);
    int get_random_seed() {return(ran_gen_ptr->get_seed());}

    //! This function set Woods-Saxon parameters based on the nucleus name
    void set_nucleus_parameters(std::string nucleus_name);
    void set_woods_saxon_parameters(int A_in, int Z_in,
                                    real rho, real w, real R, real a,
                                    int density_function_type_in);
    void set_dmin (real d_min_in) {d_min = d_min_in;}
    int get_nucleus_A() const {return(A);}
    int get_nucleus_Z() const {return(Z);}
    WoodsSaxonParam get_woods_saxon_parameters() const {return(WS_param_vec);}

    //! This function generates the spatial and momentum configurations
    //! for the nucleus
    void generate_nucleus_3d_configuration();
    //! This function samples the nucleon spatial configuration for deuteron
    void generate_deuteron_configuration();
    //! The inverse CDF of the Hulthen function for deutron wavefunction
    real get_inverse_CDF_hulthen_function(real y);
    //! The Hulthen function for deutron wavefunction
    real hulthen_function_CDF(real r);

    //! This function samples a nucleon spatial configuration according to
    //! the Fermi Distribution
    void generate_nucleus_configuration_with_woods_saxon();
    real sample_r_from_woods_saxon();
    //! Fermi Distribution 
    real fermi_distribution(real r, real R_WS, real a_WS);

    int get_number_of_nucleons() {return(nucleon_list.size());}
    Nucleon get_nucleon(int idx) {return(nucleon_list.at(idx));}
    std::vector<Nucleon> get_nucleon_list() const {return(nucleon_list);}

    void shift_nucleus(SpatialVec x_shift);
    void recenter_nucleus();

};

}

#endif  // SRC_NUCLEUS_H_
