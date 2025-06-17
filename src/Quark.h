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

    int baryon_part = 0; // represents 1/3: a baryon # of +-1 would be +-3.
    int charge_part = 0; // represents 1/3: a charge # of +-2/3 would be +-2.
    int strangeness = 0;

  public:
    Quark() = default;
    Quark(SpatialVec x_in, MomentumVec p_in) {
        set_particle_variables(x_in, p_in);
    }

    Quark(SpatialVec x_in, MomentumVec p_in, real mass_in) {
        set_particle_variables(x_in, p_in, mass_in);
    }

    Quark(SpatialVec x_in, real pdf_x_in) {
        set_pdf_x(pdf_x_in);
        MomentumVec p_in = {0.0};
        set_particle_variables(x_in, p_in);
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

    void set_baryon(int baryon_in) {baryon_part = baryon_in; }
    int get_baryon() const { return (baryon_part); }

    void set_charge(int charge_in) {charge_part = charge_in; }
    int get_charge() const { return (charge_part); }

    void set_strange(int strange_in) {strangeness = strange_in; }
    int get_strange() const { return (strangeness); } 
};

}  // namespace MCGlb

#endif  // SRC_QUARK_H_
