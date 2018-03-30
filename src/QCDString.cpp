// Copyright @ Chun Shen 2018

#include "QCDString.h"

using std::shared_ptr;

namespace MCGlb {

QCDString::QCDString(SpatialVec x_in, real tau_form_in,
                     shared_ptr<Nucleon> proj_in, shared_ptr<Nucleon> targ_in,
                     real string_tension_in) {
    x_production = x_in;
    tau_form     = tau_form_in;
    proj         = proj_in;
    targ         = targ_in;
    auto pvec    = targ.lock()->get_p();
    y_i_left     = atanh(pvec[3]/pvec[0]);
    pvec         = proj.lock()->get_p();
    y_i_right    = atanh(pvec[3]/pvec[0]);
    string_tension = string_tension_in;
}

void QCDString::evolve_QCD_string() {

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


}
