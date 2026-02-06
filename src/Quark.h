// Copyright @ Chun Shen 2018

#ifndef SRC_QUARK_H_
#define SRC_QUARK_H_

#include <cassert>

#include "Particle.h"

namespace MCGlb {

class Quark : public Particle {
  private:
    real pdf_x;
    real rapidity_q;
    bool remnant_set_ = false;
    int number_of_connections = 0;

    real baryon_ = 0;
    real charge_ = 0;
    real strange_ = 0;
    bool Q_used_ = false;

  public:
    Quark() = delete;
    Quark(SpatialVec x_in, MomentumVec p_in) {
        set_particle_variables(x_in, p_in);
        set_charges(0, 0, 0);
    }

    Quark(SpatialVec x_in, MomentumVec p_in, real mass_in) {
        set_particle_variables(x_in, p_in, mass_in);
        set_charges(0, 0, 0);
    }

    Quark(
        SpatialVec x_in, real pdf_x_in, real baryon_in = 0, real charge_in = 0,
        real strange_in = 0) {
        set_pdf_x(pdf_x_in);
        MomentumVec p_in = {0.0};
        set_particle_variables(x_in, p_in);
        set_charges(baryon_in, charge_in, strange_in);
    }

    Quark(
        SpatialVec x_in, MomentumVec p_in, real baryon_in, real charge_in,
        real strange_in) {
        set_particle_variables(x_in, p_in);
        set_charges(baryon_in, charge_in, strange_in);
    }

    Quark(
        SpatialVec x_in, MomentumVec p_in, real mass_in, real baryon_in,
        real charge_in, real strange_in) {
        set_particle_variables(x_in, p_in, mass_in);
        set_charges(baryon_in, charge_in, strange_in);
    }

    void set_pdf_x(real x_in) {
        assert(x_in >= 0.);
        assert(x_in <= 1.);
        pdf_x = x_in;
    }
    real get_pdf_x() const { return (pdf_x); }

    void set_rapidity(real rapidity_in) { rapidity_q = rapidity_in; }
    real get_rapidity() const { return (rapidity_q); }

    bool is_remnant_set() const { return (remnant_set_); }
    void set_remnant(bool remnant) { remnant_set_ = remnant; }

    void add_a_connection() { number_of_connections++; }
    int get_number_of_connections() const { return (number_of_connections); }

    void set_baryon(real baryon_in) { baryon_ = baryon_in; }
    real get_baryon() const { return (baryon_); }

    void set_charge(real charge_in) { charge_ = charge_in; }
    real get_charge() const { return (charge_); }

    void set_strange(real strange_in) { strange_ = strange_in; }
    real get_strange() const { return (strange_); }

    void set_charges(real baryon_in, real charge_in, real strange_in) {
        baryon_ = baryon_in;
        charge_ = charge_in;
        strange_ = strange_in;
    }

    void set_Q_used(bool used) { Q_used_ = used; }
    bool Q_was_used() const { return (Q_used_); }
};

}  // namespace MCGlb

#endif  // SRC_QUARK_H_
