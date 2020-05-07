// This is a small program to generate the triplets of valence quarks
// with single valence quark distribution provided by nuclear PDF
// Copyright @ Chun Shen 2020

#include <memory>
#include <array>
#include <iostream>
#include <sstream>
#include "eps09.h"
#include "Random.h"
#include "LHAPDF/LHAPDF.h"

using RandomUtil::Random;
using std::shared_ptr;
using std::array;

const int number_of_quarks = 3;
const int number_of_samples = 100000;
const double EPS = 1e-15;

typedef struct{
    array<double, number_of_quarks> Xarr;
    double score;
} Triplet;


double sample_a_u_quark_momentum_fraction(
        const bool flag_NPDF, const shared_ptr<LHAPDF::PDF> pdf,
        const shared_ptr<Random> ran_gen_ptr, const double A) {
    double x;
    double xfu, xfubar, tmp, correction;
    double ruv = 1.;
    double rdv = 1.;
    const double Q2 = 1.;
    do {
        x = ran_gen_ptr->rand_uniform();
        if (flag_NPDF) {
            double ru, rd, rs, rc, rb, rg;
            eps09(2, 1, A, x, sqrt(Q2), ruv, rdv, ru, rd, rs,
                  rc, rb, rg);
        }

        xfubar     = pdf->xfxQ2(-2, x, Q2);
        xfu        = pdf->xfxQ2( 2, x, Q2);
        tmp        = ran_gen_ptr->rand_uniform();
        correction = 1.0;
    } while (tmp > ((xfu - xfubar)*ruv*correction));
    return(x);
}


double sample_a_d_quark_momentum_fraction(
        const bool flag_NPDF, const shared_ptr<LHAPDF::PDF> pdf,
        const shared_ptr<Random> ran_gen_ptr, const double A) {
    double x;
    double xfd, xfdbar, tmp, correction;
    double ruv = 1.;
    double rdv = 1.;
    const double Q2 = 1.;
    do {
        x = ran_gen_ptr->rand_uniform();
        if (flag_NPDF) {
            double ru, rd, rs, rc, rb, rg;
            eps09(2, 1, A, x, sqrt(Q2), ruv, rdv, ru, rd, rs,
                  rc, rb, rg);
        }
        // ruv seems to be always equal to rdv,
        // so I am fine not distinguishing proton and neutron here

        xfdbar     = pdf->xfxQ2(-1, x, Q2);
        xfd        = pdf->xfxQ2( 1, x, Q2);
        tmp        = ran_gen_ptr->rand_uniform();
        correction = 1.0;
    } while (tmp > ((xfd - xfdbar)*rdv*correction));
    return(x);
}


void compute_score(Triplet &triplet_i) {
    double sum = 0.;
    for (int iq = 0; iq < number_of_quarks; iq++) {
        sum += triplet_i.Xarr[iq];
    }
    if (sum > 1.)
        triplet_i.score = 0.;
    else
        triplet_i.score = sum;
}


void swap_two_quarks(Triplet &triplet_1, Triplet &triplet_2, const int q_id) {
    double swap = triplet_1.Xarr[q_id];
    triplet_1.Xarr[q_id] = triplet_2.Xarr[q_id];
    triplet_2.Xarr[q_id] = swap;

    // recompute the scores
    compute_score(triplet_1);
    compute_score(triplet_2);
}


void one_metropolis_step(
    shared_ptr<Random> ran_int_gen_1, shared_ptr<Random> ran_int_gen_2,
    array<Triplet, number_of_samples> &quark_samples, int flag_force_violation,
    double &delta_score, int &delta_violations) {

    int sample1 = ran_int_gen_1->rand_int_uniform();
    if (flag_force_violation == 1) {
        for (int ii = 0; ii < number_of_samples; ii++) {
            if (quark_samples[ii].score < EPS) {
                sample1 = ii;
                break;
            }
        }
    }

    int sample2 = 0;
    double score_curr, score_prev;
    double sample1_prev_score, sample2_prev_score;
    do {
        do {
            sample2 = ran_int_gen_1->rand_int_uniform();
        } while (sample2 == sample1);
        int quark_id = ran_int_gen_2->rand_int_uniform();
        sample1_prev_score = quark_samples[sample1].score;
        sample2_prev_score = quark_samples[sample2].score;

        score_prev = sample1_prev_score + sample2_prev_score;
        swap_two_quarks(quark_samples[sample1], quark_samples[sample2],
                        quark_id);
        score_curr = (  quark_samples[sample1].score
                      + quark_samples[sample2].score);

        if (score_curr < score_prev) {
            // swap back
            swap_two_quarks(quark_samples[sample1], quark_samples[sample2],
                            quark_id);
            score_curr = (  quark_samples[sample1].score
                          + quark_samples[sample2].score);
        }
    } while (flag_force_violation == 1 && score_curr < score_prev);

    delta_score = (score_curr - score_prev)/number_of_samples;

    delta_violations = 0;
    if (sample1_prev_score < EPS && quark_samples[sample1].score > EPS)
        delta_violations -= 1;
    if (sample1_prev_score > EPS && quark_samples[sample1].score < EPS)
        delta_violations += 1;
    if (sample2_prev_score < EPS && quark_samples[sample2].score > EPS)
        delta_violations -= 1;
    if (sample2_prev_score > EPS && quark_samples[sample2].score < EPS)
        delta_violations += 1;
}


double compute_total_score(
        const array<Triplet, number_of_samples> &quark_samples) {
    double sum = 0.;
    for (const auto &triplet_i : quark_samples) {
        sum += triplet_i.score;
    }
    sum /= quark_samples.size();
    return(sum);
}


int number_of_violations(
        const array<Triplet, number_of_samples> &quark_samples) {
    int n_voilations = 0;
    for (const auto &triplet_i : quark_samples) {
        if (triplet_i.score < EPS) {
            n_voilations++;
        }
    }
    return(n_voilations);
}


int main(int argc, char* argv[]) {
    // create LHA pdf object
    const std::string setname = "CT10nnlo";
    const int imem = 0;
    shared_ptr<LHAPDF::PDF> pdf;
    pdf = shared_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(setname, imem));

    int A = 1;          // proton
    if (argc > 1) {
        A = std::stoi(*(argv + 1));
    }
    // define nuclear pdf
    bool flag_NPDF = false;
    if (A == 197 || A == 208) flag_NPDF = true;

    // create quark sample lists for protons and neutrons
    array<Triplet, number_of_samples> proton_quark_samples;
    array<Triplet, number_of_samples> neutron_quark_samples;

    // create random number generator
    int ran_seed = -1;
    shared_ptr<Random> ran_gen_ptr;
    shared_ptr<Random> ran_int_gen_1, ran_int_gen_2;
    ran_gen_ptr = shared_ptr<Random>(new Random(ran_seed, 0.0, 1.0));
    ran_int_gen_1 = shared_ptr<Random>(
                    new Random(ran_seed, 0, number_of_samples - 1));
    ran_int_gen_2 = shared_ptr<Random>(
                    new Random(ran_seed, 0, number_of_quarks - 1));

    for (int i = 0; i < number_of_samples; i++) {
        array<double, number_of_quarks> d_x;
        array<double, number_of_quarks> u_x;
        for (int iq = 0; iq < number_of_quarks; iq++) {
            d_x[iq] = sample_a_d_quark_momentum_fraction(flag_NPDF, pdf,
                                                         ran_gen_ptr, A);
            u_x[iq] = sample_a_u_quark_momentum_fraction(flag_NPDF, pdf,
                                                         ran_gen_ptr, A);
        }
        proton_quark_samples[i].Xarr[0] = d_x[0];
        proton_quark_samples[i].Xarr[1] = u_x[0];
        proton_quark_samples[i].Xarr[2] = u_x[1];
        neutron_quark_samples[i].Xarr[0] = d_x[1];
        neutron_quark_samples[i].Xarr[1] = d_x[2];
        neutron_quark_samples[i].Xarr[2] = u_x[2];
        compute_score(proton_quark_samples[i]);
        compute_score(neutron_quark_samples[i]);
    }

    // begin Metropolis
    double proton_total_score = compute_total_score(proton_quark_samples);
    double neutron_total_score = compute_total_score(neutron_quark_samples);
    int proton_nviolations = number_of_violations(proton_quark_samples);
    int neutron_nviolations = number_of_violations(neutron_quark_samples);
    long long int iter = 0;
    while (proton_nviolations > 0) {
        if (iter % 1000000 == 0) {
            std::cout << "proton iter = " << iter << ": nviolations = "
                      << proton_nviolations
                      << ", score = " << proton_total_score << std::endl;
        }
        double delta_p;
        int delta_violation_p;

        if (proton_nviolations > number_of_samples*0.0005) {
            one_metropolis_step(ran_int_gen_1, ran_int_gen_2,
                                proton_quark_samples, 0,
                                delta_p, delta_violation_p);
        } else {
            one_metropolis_step(ran_int_gen_1, ran_int_gen_2,
                                proton_quark_samples, 1,
                                delta_p, delta_violation_p);
        }
        proton_total_score += delta_p;
        proton_nviolations += delta_violation_p;
        iter++;
    }
    std::cout << "proton iter = " << iter << ": nviolations = "
              << proton_nviolations
              << ", score = " << proton_total_score << std::endl;

    // output to file in binary
    std::stringstream of_p_name;
    of_p_name << "tables/proton_valence_quark_samples";
    if (A == 197) {
        of_p_name << "_NPDFAu.dat";
    } else if (A == 208) {
        of_p_name << "_NPDFPb.dat";
    } else {
        of_p_name << ".dat";
    }
    std::ofstream of_p(of_p_name.str().c_str(),
                       std::ios::out | std::ios::binary | std::ofstream::app);
    for (const auto triplet_i: proton_quark_samples) {
        for (int i = 0; i < number_of_quarks; i++) {
            float x_i = static_cast<float>(triplet_i.Xarr[i]);
            of_p.write((char*) &(x_i), sizeof(float));
        }
    }
    of_p.close();

    iter = 0;
    while (neutron_nviolations > 0) {
        if (iter % 1000000 == 0) {
            std::cout << "neutron iter = " << iter << ": nviolations = "
                      << neutron_nviolations
                      << ", score = " << neutron_total_score << std::endl;
        }
        double delta_n;
        int delta_violation_n;

        if (neutron_nviolations > number_of_samples*0.0005) {
            one_metropolis_step(ran_int_gen_1, ran_int_gen_2,
                                neutron_quark_samples, 0,
                                delta_n, delta_violation_n);
        } else {
            one_metropolis_step(ran_int_gen_1, ran_int_gen_2,
                                neutron_quark_samples, 1,
                                delta_n, delta_violation_n);
        }

        neutron_total_score += delta_n;
        neutron_nviolations += delta_violation_n;
        iter++;
    }
    std::cout << "neutron iter = " << iter << ": nviolations = "
              << neutron_nviolations
              << ", score = " << neutron_total_score << std::endl;
    std::stringstream of_n_name;
    of_n_name << "tables/neutron_valence_quark_samples";
    if (A == 197) {
        of_n_name << "_NPDFAu.dat";
    } else if (A == 208) {
        of_n_name << "_NPDFPb.dat";
    } else {
        of_n_name << ".dat";
    }
    std::ofstream of_n(of_n_name.str().c_str(),
                       std::ios::out | std::ios::binary | std::ofstream::app);
    for (const auto triplet_i: neutron_quark_samples) {
        for (int i = 0; i < number_of_quarks; i++) {
            float x_i = static_cast<float>(triplet_i.Xarr[i]);
            of_n.write((char*) &(x_i), sizeof(float));
        }
    }
    of_n.close();
    return 0;
}
