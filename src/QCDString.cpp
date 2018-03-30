// Copyright @ Chun Shen 2018

#include "QCDString.h"

using std::shared_ptr;

namespace MCGlb {

QCDString::QCDString(SpatialVec x_in, real tau_form_in,
                     shared_ptr<Nucleon> proj_in, shared_ptr<Nucleon> targ_in) {
    x_production = x_in;
    tau_form     = tau_form_in;
    proj         = proj_in;
    targ         = targ_in;
    auto pvec    = targ.lock()->get_p();
    y_i_left     = atanh(pvec[3]/pvec[0]);
    pvec         = proj.lock()->get_p();
    y_i_right    = atanh(pvec[3]/pvec[0]);
}


}
