// Copyright @ Chun Shen 2018

#ifndef SRC_QCDSTRING_H_
#define SRC_QCDSTRING_H_

#include "Nucleon.h"
#include "Quark.h"
#include <memory>

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
    
    real eta_s_left, eta_s_right;
    real eta_s_baryon_left, eta_s_baryon_right;
    
    bool has_baryon_left, has_baryon_right;
    weak_ptr<Nucleon> proj;
    weak_ptr<Nucleon> targ;
    weak_ptr<Quark> proj_q;
    weak_ptr<Quark> targ_q;
    real m_over_sigma;

 public:
    QCDString() = default;
    QCDString(SpatialVec x_in, real tau_form,
              shared_ptr<Nucleon> proj_in, shared_ptr<Nucleon> targ_in,
              real string_tension_in,
              bool has_baryon_right, bool has_baryon_left);
    QCDString(SpatialVec x_in, real tau_form,
              shared_ptr<Nucleon> proj_in, shared_ptr<Nucleon> targ_in,
              shared_ptr<Quark> proj_q_in, shared_ptr<Quark> targ_q_in,
              real string_tension_in,
              bool has_baryon_right, bool has_baryon_left);

    void set_tau_form(real tau_form_in) {tau_form = tau_form_in;}
    real get_tau_form() const {return(tau_form);}

    real get_m_over_sigma() const {return(m_over_sigma);}

    void set_x_production(SpatialVec x_in) {x_production = x_in;}
    SpatialVec get_x_production() const {return(x_production);}

    void set_initial_rapidities(real y_in_l, real y_in_r) {
        y_i_left = y_in_l; y_i_right = y_in_r;
    }
    real get_y_i_left() const {return(y_i_left);}
    real get_y_i_right() const {return(y_i_right);}

    void set_final_rapidities(real y_f_l, real y_f_r) {
        y_f_left = y_f_l; y_f_right = y_f_r;
    }
    real get_y_f_left() const {return(y_f_left);}
    real get_y_f_right() const {return(y_f_right);}

    void set_final_baryon_rapidities(real y_f_l, real y_f_r) {
        y_f_baryon_left = y_f_l; y_f_baryon_right = y_f_r;
    }
    real get_y_f_baryon_left() const {return(y_f_baryon_left);}
    real get_y_f_baryon_right() const {return(y_f_baryon_right);}

    void set_final_space_time_rapidities(real eta_s_l, real eta_s_r) {
        eta_s_left = eta_s_l; eta_s_right = eta_s_r;
    }
    real get_eta_s_left() const {return(eta_s_left);}
    real get_eta_s_right() const {return(eta_s_right);}

    real get_eta_s_baryon_left() const {return(eta_s_baryon_left);}
    real get_eta_s_baryon_right() const {return(eta_s_baryon_right);}

    void set_has_baryon_left(bool has_b_left) {
        has_baryon_left = has_b_left;
    }

    void set_has_baryon_right(bool has_b_right) {
        has_baryon_right = has_b_right;
    }

    bool get_has_baryon_left() const {return(has_baryon_left);}
    bool get_has_baryon_right() const {return(has_baryon_right);}

    weak_ptr<Nucleon> get_proj() {return(proj);}
    weak_ptr<Nucleon> get_targ() {return(targ);}
    weak_ptr<Quark>   get_proj_q() {return(proj_q);}
    weak_ptr<Quark>   get_targ_q() {return(targ_q);}

    void evolve_QCD_string();
    void evolve_QCD_string_with_free_streaming();
    void evolve_QCD_string_with_constant_deceleration();
    //! this function return the final eta_s_f for free-streaming the strings
    real get_freestreaming_eta_f(real delta_tau, real y_i,
                                 real t_0, real z_0) const;

    //! this function return the final eta_s_f for constant deceleration
    //! evolution for the strings
    real get_constant_decelerate_eta_f(real m_over_sigma_in, real delta_tau,
                                       real y_i, real t_0, real z_0) const;

};

}

#endif  // SRC_QCDSTRING_H_
