// Copyright @ Chun Shen 2018

#ifndef SRC_NUCLEUS_H_
#define SRC_NUCLEUS_H_

#include <memory>
#include <string>
#include <vector>

#include "Nucleon.h"
#include "Random.h"
#include "data_structs.h"

namespace MCGlb {

class Nucleus {
  private:
    std::string name;
    int density_function_type;
    int A_;
    int Z_;
    bool deformed_;
    bool confFromFile_;
    int lightNucleusOption_;
    WoodsSaxonParam WS_param_vec;  // rho, w, R, a, beta2, beta3, beta4, gamma
    real d_min_;                   // minimum distance between nucleons
    bool sample_valence_quarks;
    real Q2;  // Q2 when sampling valence quark
    real BG_;

    std::vector<std::shared_ptr<Nucleon>> nucleon_list_;
    std::vector<std::shared_ptr<Nucleon>> participant_list_;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr;

    bool nucleon_configuration_loaded_;
    std::vector<std::vector<std::vector<float>>> heavyIon_pos_;

    std::vector<std::array<float, 3>> proton_valence_quark_x_;
    std::vector<std::array<float, 3>> neutron_valence_quark_x_;
    std::vector<std::array<float, 2>> dipole_valence_quark_x_;

    int system_status_;
    int number_of_valence_quark_samples_;
    int N_sea_partons_;

  public:
    Nucleus() = default;
    Nucleus(
        std::string nucleus_name, std::shared_ptr<RandomUtil::Random> ran_gen,
        bool sample_valence_quarks = false, real BG = 4., real d_min = 0.9,
        bool deformed = true, bool confFromFile = false, int N_sea_partons = 1);
    ~Nucleus();

    std::string get_name() const { return (name); }
    int get_random_seed() const { return (ran_gen_ptr->get_seed()); }

    int readin_valence_quark_samples();

    void setLightNucleusOption(int option) { lightNucleusOption_ = option; }

    void set_valence_quark_Q2(real Q2_q) { Q2 = Q2_q; }
    //! This function set Woods-Saxon parameters based on the nucleus name
    void set_nucleus_parameters(std::string nucleus_name);
    void set_woods_saxon_parameters(
        int A_in, int Z_in, real rho, real w, real R, real a, real beta2,
        real beta3, real beta4, real gamma, int density_function_type_in);
    void setWoodsSaxonParameters(
        real rho, real w, real R, real a, real beta2, real beta3, real beta4,
        real gamma);
    void set_dmin(real d_min) { d_min_ = d_min; }
    real get_nucleon_minimum_distance() const { return (d_min_); }
    int get_nucleus_A() const { return (A_); }
    int get_nucleus_Z() const { return (Z_); }
    WoodsSaxonParam get_woods_saxon_parameters() const {
        return (WS_param_vec);
    }
    bool is_deformed() const { return (deformed_); }

    void add_a_participant(std::shared_ptr<Nucleon> ipart) {
        if (!ipart->is_wounded()) {
            // only at ipart one time
            participant_list_.push_back(ipart);
        }
    }

    //! This function generates the spatial and momentum configurations
    //! for the nucleus
    void generate_nucleus_3d_configuration();
    //! This function samples the nucleon spatial configuration for deuteron
    void generate_deuteron_configuration();
    //! The inverse CDF of the Hulthen function for deutron wavefunction
    real get_inverse_CDF_hulthen_function(real y) const;
    //! The Hulthen function for deutron wavefunction
    real hulthen_function_CDF(real r) const;

    //! Read in spatial configuration for triton
    void readin_nucleon_positions();

    //! This function samples the spatial configuration for triton
    int sample_nucleon_configuration();

    //! This function samples a nucleon spatial configuration according to
    //! the Fermi Distribution
    void generate_nucleus_configuration_with_woods_saxon();
    void generate_nucleus_configuration_with_deformed_woods_saxon();
    real sample_r_from_woods_saxon() const;
    real sample_r_from_deformed_woods_saxon() const;
    void sample_r_and_costheta_from_deformed_woods_saxon(
        real &phi, real &r, real &costheta) const;
    //! Fermi Distribution
    real fermi_distribution(real r, real R_WS, real a_WS) const;
    real getAvgWoodsSaxonDensity(real r) const;
    real spherical_harmonics(int l, real ct) const;
    real spherical_harmonics_Y22(int l, real ct, real phi) const;

    int get_number_of_nucleons() const { return (nucleon_list_.size()); }
    std::shared_ptr<Nucleon> get_nucleon(unsigned int idx) {
        return (nucleon_list_.at(idx));
    }
    std::vector<std::shared_ptr<Nucleon>> *get_nucleon_list() {
        return (&nucleon_list_);
    }
    int get_number_of_wounded_nucleons() const {
        return (static_cast<int>(participant_list_.size()));
    }
    std::shared_ptr<Nucleon> get_participant(unsigned int idx) {
        return (participant_list_.at(idx));
    }

    void shift_nucleus(SpatialVec x_shift);
    void recenter_nucleus();
    void rotate_nucleus(real phi, real theta);
    void rotate_nucleus_3D(real phi, real theta, real gamma);

    void accelerate_nucleus(real ecm, int direction);
    void accelerate_dipole(real ecm, int direction);
    void lorentz_contraction(real gamma);
    void set_nucleons_momentum_with_collision_energy(real beam_rapidity);
    void set_dipole_momentum_with_collision_energy(real beam_rapidity);
    real get_z_min() const;
    real get_z_max() const;

    void output_nucleon_positions(std::string filename) const;

    void sample_valence_quarks_inside_nucleons(real ecm, int direction);
    void add_soft_parton_ball(real ecm, int direction);

    void sample_fermi_momentum();

    void sample_quark_momentum_fraction(
        std::vector<real> &xQuark, const int number_of_quarks,
        const int nucleonType, const real mass, const real ecm) const;
    void sample_quark_momentum_fraction_in_dipole(
        std::vector<real> &xQuark, const int number_of_quarks,
        const real ecm) const;

    SpatialVec sample_valence_quark_position() const;
    real ExponentialDistribution(const real a, const real r) const;
};

}  // namespace MCGlb

#endif  // SRC_NUCLEUS_H_
