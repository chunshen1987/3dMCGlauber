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

    real baryon = 0;
    real charge = 0;
    real strange = 0;
    bool baryon_used = false;
    bool remnant_carry_baryon_number_ = false;
    real remnant_baryon_number = 1.0/3.0;

  public:
    Quark() = default;
    Quark(SpatialVec x_in, MomentumVec p_in) {
        set_particle_variables(x_in, p_in);
        set_charges(0, 0, 0);
    }

    Quark(SpatialVec x_in, MomentumVec p_in, real mass_in) {
        set_particle_variables(x_in, p_in, mass_in);
        set_charges(0, 0, 0);
    }

     Quark(SpatialVec x_in, real pdf_x_in, real baryon_in = 0, real charge_in = 0, real strange_in = 0) {
        set_pdf_x(pdf_x_in);
        MomentumVec p_in = {0.0};
        set_particle_variables(x_in, p_in);
        set_charges(baryon_in, charge_in, strange_in);
    }

    Quark(SpatialVec x_in, MomentumVec p_in, real baryon_in, real charge_in, real strange_in) {
        set_particle_variables(x_in, p_in);
        set_charges(baryon_in, charge_in, strange_in);
    }

    Quark(SpatialVec x_in, MomentumVec p_in, real mass_in,  real baryon_in, real charge_in, real strange_in) {
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

    void set_baryon(real baryon_in) {baryon = baryon_in; }
    real get_baryon() const { return (baryon); }

    void set_charge(real charge_in) {charge = charge_in; }
    real get_charge() const { return (charge); }

    void set_strange(real strange_in) {strange = strange_in; }
    real get_strange() const { return (strange); }

    void set_charges(real baryon_in, real charge_in, real strange_in) {
        baryon = baryon_in;
        charge = charge_in;
        strange = strange_in;
    }

    void set_baryon_used(bool used) {baryon_used = used; }
    bool baryon_was_used() const { return (baryon_used);}
    bool is_remnant_carry_baryon_number() const {
        return (remnant_carry_baryon_number_);
    }
    real get_remnant_baryon_number() const {
        return (remnant_baryon_number);
    }
    void set_remnant_carry_baryon_number(bool remnant) {
        remnant_carry_baryon_number_ = remnant;
    }
    void set_remnant_baryon_number(real remnant) {
        remnant_baryon_number = remnant;
    }
};

}  // namespace MCGlb

#endif  // SRC_QUARK_H_
