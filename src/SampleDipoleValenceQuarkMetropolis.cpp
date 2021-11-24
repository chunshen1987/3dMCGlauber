// This is a small program to generate the double quarks
// with single valence quark distribution p(x)=x^alpha(1-x)^beta
// Copyright @ Wenbin Zhao, Chun Shen 2021

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
const int two_quarks = 2;
const int number_of_samples = 100000;
const double acc_violation_fraction = 5e-4;
const double allow_violation_fraction = 0;
const long int ntol = 5000000;
const double EPS = 1e-15;

typedef struct{
    array<double, number_of_quarks> Xarr;
    double score;
} Triplet;

typedef struct{
    array<double, two_quarks> Xarr2;
    double score;
} Double_let;

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

// sample a (anti-)quark's momentum fraction of the meson or dipole
double sample_a_quark_momentum_fraction_in_diople(
         double CDF[], int size,const shared_ptr<Random> ran_gen_ptr,double dx) {
    double x;       
    int ndivided=20; 
    int index1=size/ndivided;
    int startid;
    int simble=0;
    int simble2=0;
    double tmp = ran_gen_ptr->rand_uniform();
    if(tmp<=CDF[index1]){
      startid=index1*0;
      simble2=1;
    }
    if(tmp>CDF[index1*ndivided-index1]){
      startid=index1*ndivided-index1;
      simble2=1;
    }
    if(simble2==0){
      for (int k=1;k<ndivided-1;k++){
          if(tmp>CDF[index1*k] && tmp<=CDF[index1*k+index1]){startid=index1*k;simble2=1;}
          if(simble2==1)break;
      }
    }
    
    for (int i = startid; i < index1+startid; i++) {
        if(tmp<CDF[i+1]){
           x=i*1.0*dx+dx/2.0;
           simble=1;
        }
        if(simble==1)break;
    }
    if(simble==0)x=0.99;
    //std::cout<<"priliminary x "<<x<<" "<<simble<<std::endl;
    return (x);
}


void compute_score(Double_let &triplet_i) {
    double sum = 0.;
    for (int iq = 0; iq < two_quarks; iq++) {
        sum += triplet_i.Xarr2[iq];
    }
    if (sum > 1.)
        triplet_i.score = 0.;
    else
        triplet_i.score = sum;
}


void swap_two_quarks(Double_let &triplet_1, Double_let &triplet_2, const int q_id) {
    double swap = triplet_1.Xarr2[q_id];
    triplet_1.Xarr2[q_id] = triplet_2.Xarr2[q_id];
    triplet_2.Xarr2[q_id] = swap;

    // recompute the scores
    compute_score(triplet_1);
    compute_score(triplet_2);
}


void one_metropolis_step(
    shared_ptr<Random> ran_int_gen_1, shared_ptr<Random> ran_int_gen_2,
    array<Double_let, number_of_samples> &quark_samples, int flag_force_violation,
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
        const array<Double_let, number_of_samples> &quark_samples) {
    double sum = 0.;
    for (const auto &triplet_i : quark_samples) {
        sum += triplet_i.score;
    }
    sum /= quark_samples.size();
    return(sum);
}


int number_of_violations(
        const array<Double_let, number_of_samples> &quark_samples) {
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
    // the quark's PDF in the dipole, p(x)=x^alpha(1-x)^beta
    double dx=0.002;
    int lenght=1/dx;
    double CDF[lenght+100]={0.0};
    double alpha=2.0;
    double beta=2.0;
    double loopx=0.0;
    int index=0;
    while(loopx<1.0){
        double dipole_px=pow(loopx,alpha)*pow(1-loopx,beta)*dx;
        if(index==0){
            CDF[index]=dipole_px;
        }else{
            CDF[index]=dipole_px+CDF[index-1];
        }
        
        loopx=loopx+dx;
        index++;
    }
    for(int i=0;i<index;i++){
        CDF[i]=CDF[i]/CDF[index-1];
        //std::cout<<" "<<i<<" "<<CDF[i]<<std::endl;
    }
    
    // define nuclear pdf
    bool flag_NPDF = false;
    if (A == 197 || A == 208) flag_NPDF = true;

    // create quark sample lists for protons and neutrons
    array<Double_let, number_of_samples> dipole_quark_samples;
    //array<Triplet, number_of_samples> neutron_quark_samples;

    // create random number generator
    int ran_seed = -1;
    shared_ptr<Random> ran_gen_ptr;
    shared_ptr<Random> ran_int_gen_1, ran_int_gen_2;
    ran_gen_ptr = shared_ptr<Random>(new Random(ran_seed, 0.0, 1.0));
    ran_int_gen_1 = shared_ptr<Random>(
                    new Random(ran_seed, 0, number_of_samples - 1));
    ran_int_gen_2 = shared_ptr<Random>(
                    new Random(ran_seed, 0, two_quarks - 1));

    for (int i = 0; i < number_of_samples; i++) {
        array<double, two_quarks> d_x;
        for (int iq = 0; iq < two_quarks; iq++) {
            d_x[iq] = sample_a_quark_momentum_fraction_in_diople(CDF, index,
                                                         ran_gen_ptr, dx);
        }
        dipole_quark_samples[i].Xarr2[0] = d_x[0];
        dipole_quark_samples[i].Xarr2[1] = d_x[1];
        compute_score(dipole_quark_samples[i]);
    }

    // begin Metropolis
    double dipole_total_score = compute_total_score(dipole_quark_samples);
    int dipole_nviolations = number_of_violations(dipole_quark_samples);
    long long int iter = 0;
    long int itol = 0;
    while (dipole_nviolations > 0 && itol < ntol) {
        
        if (iter % 1000000 == 0) {
            std::cout << "dipole iter = " << iter << ": nviolations = "
                      << dipole_nviolations
                      << ", <sum_x> = " << dipole_total_score << std::endl;
        }
        
        double delta_p;
        int delta_violation_p;

        if (dipole_nviolations > number_of_samples*acc_violation_fraction) {
            one_metropolis_step(ran_int_gen_1, ran_int_gen_2,
                                dipole_quark_samples, 0,
                                delta_p, delta_violation_p);
        } else if (dipole_nviolations
                   > number_of_samples*allow_violation_fraction) {
            one_metropolis_step(ran_int_gen_1, ran_int_gen_2,
                                dipole_quark_samples, 1,
                                delta_p, delta_violation_p);
            itol++;
        } else {
            break;
        }

        dipole_total_score += delta_p;
        dipole_nviolations += delta_violation_p;

        iter++;
    }
    
    std::cout << "dipole iter = " << iter << ": nviolations = "
              << dipole_nviolations
              << ", <sum_x> = " << dipole_total_score << std::endl;
    
    // output to file in binary
    std::stringstream of_p_name;
    of_p_name << "tables/dipole_valence_quark_samples";
    if (A == 197) {
        of_p_name << "_NPDFAu.dat";
    } else if (A == 208) {
        of_p_name << "_NPDFPb.dat";
    } else {
        of_p_name << ".dat";
    }
    std::ofstream of_p(of_p_name.str().c_str(),
                       std::ios::out | std::ios::binary | std::ofstream::app);
    for (const auto triplet_i: dipole_quark_samples) {
        if (triplet_i.score < 1.) {
            for (int i = 0; i < two_quarks; i++) {
                float x_i = static_cast<float>(triplet_i.Xarr2[i]);
                of_p.write((char*) &(x_i), sizeof(float));
            }
        }
    }
    of_p.close();
    return 0;
}
