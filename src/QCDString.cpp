// Copyright @ Chun Shen 2018

#include "QCDString.h"

using std::shared_ptr;

namespace MCGlb {

QCDString::QCDString(SpatialVec x_in, real tau_form_in,
                     shared_ptr<Nucleon> proj_in, shared_ptr<Nucleon> targ_in,
                     real m_over_sigma_in,
                     bool has_baryon_right_in, bool has_baryon_left_in) {
    x_production       = x_in;
    tau_form           = tau_form_in;
    proj               = proj_in;
    targ               = targ_in;
    auto pvec          = targ.lock()->get_p();
    y_i_left           = atanh(pvec[3]/pvec[0]);
    pvec               = proj.lock()->get_p();
    y_i_right          = atanh(pvec[3]/pvec[0]);
    m_over_sigma       = m_over_sigma_in;
    has_baryon_right_  = has_baryon_right_in;
    has_baryon_left_   = has_baryon_left_in;
    eta_s_baryon_left  = 0.;
    eta_s_baryon_right = 0.;
    has_remnant_left_  = false;
    has_remnant_right_ = false;
}

QCDString::QCDString(SpatialVec x_in, real tau_form_in,
                     shared_ptr<Nucleon> proj_in, shared_ptr<Nucleon> targ_in,
                     shared_ptr<Quark> proj_q_in, shared_ptr<Quark> targ_q_in,
                     real m_over_sigma_in,
                     bool has_baryon_right_in, bool has_baryon_left_in) {
    x_production       = x_in;
    tau_form           = tau_form_in;
    proj               = proj_in;
    targ               = targ_in;
    proj_q             = proj_q_in;
    targ_q             = targ_q_in;
    y_i_left           = targ_q.lock()->get_rapidity();
    y_i_right          = proj_q.lock()->get_rapidity();
    m_over_sigma       = m_over_sigma_in;
    has_baryon_right_  = has_baryon_right_in;
    has_baryon_left_   = has_baryon_left_in;
    eta_s_baryon_left  = 0.;
    eta_s_baryon_right = 0.;
    has_remnant_left_  = false;
    has_remnant_right_ = false;
}

void QCDString::evolve_QCD_string() {
    if (m_over_sigma > 1e6) {
        evolve_QCD_string_with_free_streaming();
    } else {
        evolve_QCD_string_with_constant_deceleration();
    }
}

void QCDString::evolve_QCD_string_with_free_streaming() {
    // freestream the string by its formation time tau_form
    y_f_left    = y_i_left;
    y_f_right   = y_i_right;
    eta_s_left  = get_freestreaming_eta_f(tau_form, y_i_left, x_production[0],
                                          x_production[3]);
    eta_s_right = get_freestreaming_eta_f(tau_form, y_i_right, x_production[0],
                                          x_production[3]);
}

void QCDString::evolve_QCD_string_with_constant_deceleration() {
    real y_i_lrf = std::abs(y_i_right - y_i_left)/2.;
    real deceleration_factor = (tau_form*tau_form
                                /(2.*m_over_sigma*m_over_sigma));
    real y_loss = acosh(deceleration_factor + 1.);
    if (y_loss > y_i_lrf) {
        tau_form = m_over_sigma*sqrt(2.*(cosh(y_i_lrf) - 1.));
        y_loss   = y_i_lrf;
    }
    y_f_left    = y_i_left + y_loss;
    y_f_right   = y_i_right - y_loss;
    eta_s_left  = get_constant_decelerate_eta_f(-m_over_sigma, tau_form,
                                                y_i_left, x_production[0],
                                                x_production[3]);
    eta_s_right = get_constant_decelerate_eta_f(m_over_sigma, tau_form,
                                                y_i_right, x_production[0],
                                                x_production[3]);
}

real QCDString::get_freestreaming_eta_f(real delta_tau, real y_i,
                                        real t_0, real z_0) const {
    const real cosh_y_i = cosh(y_i);
    const real cosh_y_i_sq = cosh_y_i*cosh_y_i;
    const real temp_factor = t_0 - z_0*tanh(y_i);
    const real tau_0 = sqrt(t_0*t_0 - z_0*z_0);
    const real tau_f = tau_0 + delta_tau;
    const real dt = (cosh_y_i_sq*(
        - temp_factor + sqrt(temp_factor*temp_factor
                        + (tau_f*tau_f - tau_0*tau_0)/cosh_y_i_sq)));
    const real t_f = t_0 + dt;
    const real z_f = z_0 + dt*tanh(y_i);
    const real eta_s_f = 0.5*log((t_f + z_f)/(t_f - z_f));
    return(eta_s_f);
}

//! this function return the final eta_s_f for constant deceleration evolution
//! for the strings
real QCDString::get_constant_decelerate_eta_f(
    real m_over_sigma_in, real delta_tau, real y_i, real t_0, real z_0) const {
    const real cosh_y = cosh(y_i);
    const real sinh_y = sinh(y_i);
    const real temp = delta_tau/(2.*m_over_sigma_in);
    const real t_f = t_0 + delta_tau*(-temp*sinh_y + sqrt(temp*temp + 1)*cosh_y);
    const real z_f = z_0 + delta_tau*(-temp*cosh_y + sqrt(temp*temp + 1)*sinh_y);
    const real eta_s_f = 0.5*log((t_f + z_f)/(t_f - z_f));
    return(eta_s_f);
}


//! This funciton sets baryon eta_s assuming linear rapidity profile
void QCDString::set_final_baryon_space_time_rapidities() {
    const double slope = (eta_s_right - eta_s_left)/(y_f_right - y_f_left);
    eta_s_baryon_left  = eta_s_left + slope*(y_f_baryon_left  - y_f_left);
    eta_s_baryon_right = eta_s_left + slope*(y_f_baryon_right - y_f_left);
}


}
