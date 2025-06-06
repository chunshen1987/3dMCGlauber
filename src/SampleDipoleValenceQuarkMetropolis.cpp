// This is a small program to generate the double quarks
// with single valence quark distribution p(x)=x^alpha(1-x)^beta
// Copyright @ Wenbin Zhao, Chun Shen 2021

#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#include "Random.h"
using RandomUtil::Random;
using std::array;
using std::shared_ptr;

const int two_quarks = 2;
const int number_of_samples = 200000;
const double acc_violation_fraction = 2e-3;
const double allow_violation_fraction = 0.0;
const long int ntol = 10000000;
const double EPS = 1e-15;

typedef struct {
    array<double, two_quarks> Xarr2;
    double score;
} Double_let;

double binary_search(
    std::vector<double> &yarr, std::vector<double> &xarr, double key) {
    int start = 0;
    int end = yarr.size();
    int mid = 0;
    while (end - start > 1) {
        mid = static_cast<int>((start + end) / 2);
        if (yarr[mid] > key) {
            end = mid;
        } else {
            start = mid;
        }
    }
    double ret =
        (xarr[start]
         + (xarr[start] - xarr[end]) / (yarr[start] - yarr[end])
               * (key - yarr[start]));
    return ret;
}

void compute_score(Double_let &doublet_i) {
    double sum = 0.;
    for (int iq = 0; iq < two_quarks; iq++) {
        sum += doublet_i.Xarr2[iq];
    }
    if (sum > 1.)
        doublet_i.score = 0.;
    else
        doublet_i.score = sum;
}

void swap_two_quarks(
    Double_let &doublet_1, Double_let &doublet_2, const int q_id) {
    double swap = doublet_1.Xarr2[q_id];
    doublet_1.Xarr2[q_id] = doublet_2.Xarr2[q_id];
    doublet_2.Xarr2[q_id] = swap;

    // recompute the scores
    compute_score(doublet_1);
    compute_score(doublet_2);
}

void one_metropolis_step(
    shared_ptr<Random> ran_int_gen_1, shared_ptr<Random> ran_int_gen_2,
    array<Double_let, number_of_samples> &quark_samples,
    int flag_force_violation, double &delta_score, int &delta_violations) {
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
    int itol = 0;
    do {
        do {
            sample2 = ran_int_gen_1->rand_int_uniform();
        } while (sample2 == sample1);
        int quark_id = ran_int_gen_2->rand_int_uniform();
        sample1_prev_score = quark_samples[sample1].score;
        sample2_prev_score = quark_samples[sample2].score;

        score_prev = sample1_prev_score + sample2_prev_score;
        swap_two_quarks(
            quark_samples[sample1], quark_samples[sample2], quark_id);
        score_curr =
            (quark_samples[sample1].score + quark_samples[sample2].score);

        if (score_curr < score_prev) {
            // swap back
            swap_two_quarks(
                quark_samples[sample1], quark_samples[sample2], quark_id);
            score_curr =
                (quark_samples[sample1].score + quark_samples[sample2].score);
        }
        itol++;
    } while (flag_force_violation == 1 && score_curr < score_prev && itol < 10);

    delta_score = (score_curr - score_prev) / number_of_samples;

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
    const array<Double_let, number_of_samples> &quark_samples) {
    double sum = 0.;
    for (const auto &doublet_i : quark_samples) {
        sum += doublet_i.score;
    }
    sum /= quark_samples.size();
    return (sum);
}

int number_of_violations(
    const array<Double_let, number_of_samples> &quark_samples) {
    int n_voilations = 0;
    for (const auto &doublet_i : quark_samples) {
        if (doublet_i.score < EPS) {
            n_voilations++;
        }
    }
    return (n_voilations);
}

int main(int argc, char *argv[]) {
    // the quark's PDF in the dipole, p(x)=x^alpha(1-x)^beta
    double dx = 0.01;
    int length = static_cast<int>(1. / dx) + 1;

    double Alpha = 2.0;
    double Beta = 2.0;
    std::vector<double> xarr(length, 0.);
    std::vector<double> CDF(length, 0.);

    for (unsigned int i = 1; i < CDF.size(); i++) {
        double loopx = (i - 0.5) * dx;
        double fbeta = pow(loopx, Alpha) * pow(1.0 - loopx, Beta);
        CDF[i] = (fbeta * dx + CDF[i - 1]);
        xarr[i] = i * dx;
    }

    double norm = std::beta(Alpha + 1.0, Beta + 1.0);
    for (unsigned int i = 1; i < CDF.size(); i++) {
        CDF[i] /= norm;
    }

    // create quark sample lists for protons and neutrons
    array<Double_let, number_of_samples> dipole_quark_samples;

    // create random number generator
    int ran_seed = -1;
    shared_ptr<Random> ran_gen_ptr;
    shared_ptr<Random> ran_int_gen_1, ran_int_gen_2;
    ran_gen_ptr = shared_ptr<Random>(new Random(ran_seed, 0.0, 1.0));
    ran_int_gen_1 =
        shared_ptr<Random>(new Random(ran_seed, 0, number_of_samples - 1));
    ran_int_gen_2 = shared_ptr<Random>(new Random(ran_seed, 0, two_quarks - 1));
    for (int i = 0; i < number_of_samples; i++) {
        array<double, two_quarks> quark_x;
        for (int iq = 0; iq < two_quarks; iq++) {
            double tmp = ran_gen_ptr->rand_uniform();
            quark_x[iq] = binary_search(CDF, xarr, tmp);
        }
        dipole_quark_samples[i].Xarr2[0] = quark_x[0];
        dipole_quark_samples[i].Xarr2[1] = quark_x[1];
        compute_score(dipole_quark_samples[i]);
    }

    // begin Metropolis
    double dipole_total_score = compute_total_score(dipole_quark_samples);
    int dipole_nviolations = number_of_violations(dipole_quark_samples);
    long long int iter = 0;
    long int itol = 0;
    while (dipole_nviolations > 0 && itol < ntol) {
        if (iter % 1000000 == 0) {
            std::cout << "dipole iter = " << iter
                      << ": nviolations = " << dipole_nviolations
                      << ", <sum_x> = " << dipole_total_score << std::endl;
        }
        double delta_p;
        int delta_violation_p;

        if (dipole_nviolations > number_of_samples * acc_violation_fraction) {
            one_metropolis_step(
                ran_int_gen_1, ran_int_gen_2, dipole_quark_samples, 0, delta_p,
                delta_violation_p);
        } else if (
            dipole_nviolations > number_of_samples * allow_violation_fraction) {
            one_metropolis_step(
                ran_int_gen_1, ran_int_gen_2, dipole_quark_samples, 1, delta_p,
                delta_violation_p);
            itol++;
        } else {
            break;
        }

        dipole_total_score += delta_p;
        dipole_nviolations += delta_violation_p;

        iter++;
    }

    std::cout << "dipole iter = " << iter
              << ": nviolations = " << dipole_nviolations
              << ", <sum_x> = " << dipole_total_score << std::endl;
    // output to file in binary
    std::stringstream of_p_name;
    of_p_name << "tables/dipole_valence_quark_samples.dat";
    std::ofstream of_p(
        of_p_name.str().c_str(),
        std::ios::out | std::ios::binary | std::ofstream::app);
    for (const auto doublet_i : dipole_quark_samples) {
        if (doublet_i.score > 0.) {
            for (int i = 0; i < two_quarks; i++) {
                float x_i = static_cast<float>(doublet_i.Xarr2[i]);
                of_p.write((char *)&(x_i), sizeof(float));
            }
        }
    }
    of_p.close();
    return 0;
}
