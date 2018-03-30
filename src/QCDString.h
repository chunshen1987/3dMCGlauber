// Copyright @ Chun Shen 2018

#ifndef SRC_QCDSTRING_H_
#define SRC_QCDSTRING_H_

#include "Nucleon.h"
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
    real eta_s_left, eta_s_right;
    weak_ptr<Nucleon> proj;
    weak_ptr<Nucleon> targ;

 public:
    QCDString() = default;
    QCDString(SpatialVec x_in, real tau_form,
              shared_ptr<Nucleon> proj, shared_ptr<Nucleon> targ);

    void set_tau_form(real tau_form_in) {tau_form = tau_form_in;}
    real get_tau_form() const {return(tau_form);}

    void set_x_production(SpatialVec x_in) {x_production = x_in;}
    SpatialVec get_x_production() const {return(x_production);}

    void set_initial_rapidities(real y_in_l, real y_in_r) {
        y_i_left = y_in_l; y_i_right = y_in_r;
    }
    void set_final_rapidities(real y_f_l, real y_f_r) {
        y_f_left = y_f_l; y_f_right = y_f_r;
    }
    void set_final_space_time_rapidities(real eta_s_l, real eta_s_r) {
        eta_s_left = eta_s_l; eta_s_right = eta_s_r;
    }

};

}

#endif  // SRC_QCDSTRING_H_
