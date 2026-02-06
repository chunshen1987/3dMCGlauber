// Copyright @ Chun Shen 2018

#ifndef SRC_QCDSTRING_H_
#define SRC_QCDSTRING_H_

#include <memory>

#include "Nucleon.h"
#include "Quark.h"

using std::shared_ptr;
using std::weak_ptr;

namespace MCGlb {

class QCDString {
  private:
    SpatialVec x_production;
    real tau_form;
    real y_i_left, y_i_right;
    real y_f_left, y_f_right;
    real y_f_baryon_left, y_f_baryon_right;
    real y_f_electric_charge_left, y_f_electric_charge_right;

    real eta_s_left, eta_s_right;
    real eta_s_baryon_left, eta_s_baryon_right;
    real eta_s_electric_charge_left, eta_s_electric_charge_right;

    real baryon_charge_left_, baryon_charge_right_;
    real electric_charge_left_, electric_charge_right_;

    bool has_remnant_left_, has_remnant_right_;

    shared_ptr<Nucleon> proj;
    shared_ptr<Nucleon> targ;
    shared_ptr<Quark> proj_q;
    shared_ptr<Quark> targ_q;
    real m_over_sigma;
    real mass_;

  public:
    QCDString() = delete;

    QCDString(
        SpatialVec x_in, real tau_form_in, shared_ptr<Nucleon> proj_in,
        shared_ptr<Nucleon> targ_in, real m_over_sigma_in,
        real baryon_right_in = 0, real baryon_left_in = 0,
        real electric_charge_right_in = 0, real electric_charge_left_in = 0);

    QCDString(
        SpatialVec x_in, real tau_form_in, shared_ptr<Nucleon> proj_in,
        shared_ptr<Nucleon> targ_in, shared_ptr<Quark> proj_q_in,
        shared_ptr<Quark> targ_q_in, real m_over_sigma_in,
        real baryon_right_in = 0, real baryon_left_in = 0,
        real electric_charge_right_in = 0, real electric_charge_left_in = 0);

    QCDString(
        SpatialVec x_in, real tau_form_in, shared_ptr<Nucleon> proj_in,
        shared_ptr<Nucleon> targ_in, MomentumVec proj_p_in,
        MomentumVec targ_p_in, real m_over_sigma_in, real baryon_right_in = 0,
        real baryon_left_in = 0, real electric_charge_right_in = 0,
        real electric_charge_left_in = 0);

    real get_string_mass() const { return (mass_); }

    void set_tau_form(real tau_form_in) { tau_form = tau_form_in; }
    real get_tau_form() const { return (tau_form); }

    real get_m_over_sigma() const { return (m_over_sigma); }

    void set_x_production(SpatialVec x_in) { x_production = x_in; }
    SpatialVec get_x_production() const { return (x_production); }

    void set_initial_rapidities(real y_in_l, real y_in_r) {
        y_i_left = y_in_l;
        y_i_right = y_in_r;
    }
    real get_y_i_left() const { return (y_i_left); }
    real get_y_i_right() const { return (y_i_right); }

    void set_final_rapidities(real y_f_l, real y_f_r) { y_f_right = y_f_r; }
    real get_y_f_left() const { return (y_f_left); }
    real get_y_f_right() const { return (y_f_right); }

    void set_final_baryon_rapidities(real y_f_l, real y_f_r) {
        y_f_baryon_left = y_f_l;
        y_f_baryon_right = y_f_r;
        set_final_baryon_space_time_rapidities();
    }
    void set_final_electric_charge_rapidities(real y_f_l, real y_f_r) {
        y_f_electric_charge_left = y_f_l;
        y_f_electric_charge_right = y_f_r;
        set_final_electric_charge_space_time_rapidities();
    }

    real get_y_f_baryon_left() const { return (y_f_baryon_left); }
    real get_y_f_baryon_right() const { return (y_f_baryon_right); }

    real get_y_f_electric_charge_left() const {
        return (y_f_electric_charge_left);
    }
    real get_y_f_electric_charge_right() const {
        return (y_f_electric_charge_right);
    }

    void set_final_space_time_rapidities(real eta_s_l, real eta_s_r) {
        eta_s_left = eta_s_l;
        eta_s_right = eta_s_r;
    }
    real get_eta_s_left() const { return (eta_s_left); }
    real get_eta_s_right() const { return (eta_s_right); }

    void set_final_baryon_space_time_rapidities(real eta_s_l, real eta_s_r) {
        // set with inputs
        eta_s_baryon_left = eta_s_l;
        eta_s_baryon_right = eta_s_r;
    }
    void set_final_electric_charge_space_time_rapidities(
        real eta_s_l, real eta_s_r) {
        // set with inputs
        eta_s_electric_charge_left = eta_s_l;
        eta_s_electric_charge_right = eta_s_r;
    }

    //! This funciton sets baryon eta_s assuming linear rapidity profile
    void set_final_baryon_space_time_rapidities();

    //! This funciton sets baryon eta_s assuming linear rapidity profile
    void set_final_electric_charge_space_time_rapidities();

    real get_eta_s_electric_charge_left() const {
        return (eta_s_electric_charge_left);
    }
    real get_eta_s_electric_charge_right() const {
        return (eta_s_electric_charge_right);
    }

    real get_eta_s_baryon_left() const { return (eta_s_baryon_left); }
    real get_eta_s_baryon_right() const { return (eta_s_baryon_right); }

    void set_baryon_charge_left(real b_left) { baryon_charge_left_ = b_left; }

    void set_baryon_charge_right(real b_right) {
        baryon_charge_right_ = b_right;
    }
    real get_baryon_charge_left() const { return (baryon_charge_left_); }
    real get_baryon_charge_right() const { return (baryon_charge_right_); }

    void set_electric_charge_left(real eQ_left) {
        electric_charge_left_ = eQ_left;
    }
    void set_electric_charge_right(real eQ_right) {
        electric_charge_right_ = eQ_right;
    }
    real get_electric_charge_left() const { return (electric_charge_left_); }
    real get_electric_charge_right() const { return (electric_charge_right_); }

    void set_has_remnant_left(bool has_r_left) {
        has_remnant_left_ = has_r_left;
    }

    void set_has_remnant_right(bool has_r_right) {
        has_remnant_right_ = has_r_right;
    }

    bool get_has_remnant_left() const { return (has_remnant_left_); }
    bool get_has_remnant_right() const { return (has_remnant_right_); }

    shared_ptr<Nucleon> get_proj() { return (proj); }
    shared_ptr<Nucleon> get_targ() { return (targ); }
    shared_ptr<Quark> get_proj_q() { return (proj_q); }
    shared_ptr<Quark> get_targ_q() { return (targ_q); }

    void evolve_QCD_string();
    void evolve_QCD_string_with_free_streaming();
    void evolve_QCD_string_with_constant_deceleration();
    //! this function return the final eta_s_f for free-streaming the strings
    real get_freestreaming_eta_f(
        real delta_tau, real y_i, real t_0, real z_0) const;

    //! this function return the final eta_s_f for constant deceleration
    //! evolution for the strings
    real get_constant_decelerate_eta_f(
        real m_over_sigma_in, real delta_tau, real y_i, real t_0,
        real z_0) const;
};

}  // namespace MCGlb

#endif  // SRC_QCDSTRING_H_
