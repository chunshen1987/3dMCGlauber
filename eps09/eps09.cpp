//=================================================================================
//
//	File: eps09.cxx
//  Author: Thomas Ullrich (thomas.ullrich@bnl.gov)
//  Last modified: June 16, 2009
//
//  This is the C++ version of he original EPS09.f. The original description
//  from EPS09.f is below. In the translation to C++ I tried to stay as close
//  as possible to the original. All variable names are kept the same and the
//  arrays are handled as in the Fortran version (hence all arrays are longer
//  by one and the first index (0) is ignored. Error messages are the same as
//  in the original.
//
//  Syntax:
//
//  void eps09(int order, int pset, int AAA,
//             double xxx, double QQQ, double& ruv,
//             double& rdv, double& ru, double& rd,
//             double& rs, double& rc, double& rb,
//             double& rg);
//    
//  For details see below.
//=================================================================================
//
// An interface for the scale dependent nuclear modifications
// 		R_f^A(x,Q) = f_A(x,Q)/f_p(x,Q) 
// where f_A is the distribution of the parton flavour f for a PROTON in a
// nucleus A, and f_p is the corresponding parton distribution in the 
// free proton.
//  
// When using this interface, please refer to:
//  
// K.J. Eskola, H. Paukkunen and C.A. Salgado,
// "EPS09 - a New Generation of NLO and LO Nuclear Parton Distribution Functions,"
// Published as JHEP04(2009) 065.
// Eprint: arXiv:0902.4154 [hep-ph].
//
// Questions & comments to:
//   hannu.paukkunen@phys.jyu.fi
//   kari.eskola@phys.jyu.fi
//   carlos.salgado@usc.es
// 
//=================================================================================
// Instructions:
//
// For given input values of
//
//     order: 1=LO, 2=NLO   ; integer
//     pset : 1...31        ; integer
//            1     = central fit
//            2,3   = error sets S{+1}, S{-1}
//            4,5   = error sets S{+2}, S{-2}
//            ...   ...
//            30,31 = error sets {S+15}, {S-15}
//     A    : atomic number ; integer
//     x    : Bjorken-x     ; double precision
//     Q    : scale in GeV  ; double precision
//
// the command 
//
//   Call EPS09(order, pset, A, x, Q, ruv, rdv, ru, rd, rs, rc, rb, rg)
//
// returns the bound proton nuclear corrections R_f^A(x,Q)
// (in double precision) for
//	
//	ruv = up valence
//	rdv = down valence
//	ru  = up sea
//	rd  = down sea
//	rs  = strange
//	rc  = charm
//	rb  = bottom
//	rg  = gluons
//
// The nuclear corrections for bound neutrons can be obtained
// by the isospin symmetry, e.g. the total up quark distribution
// per nucleon in a nucleus A with Z protons is
//
//  u_A(x,Q) =    Z/A * [ruv*uV_p(x,Q) + ru*uSea_p(x,Q)] +
//            (A-Z)/A * [rdv*dV_p(x,Q) + rd*dSea_p(x,Q)]
//
// Note that the parametrization should only be applied at the
// kinematical domain
//
//             1e-6 <= x <= 1
//              1.3 <= Q <= 1000 GeV.
//
// No warning message is displayed if these limits are
// exceeded, and outside these boundaries the modifications
// are frozen to the boundary values, i.e
//
//   for Q > 1000, the modifications at Q=1000 are returned,
//   for Q < 1.3,  the modifications at Q=1.3 are returned,
//   for x > 1,    the modifications at x=1 are returned
//   for x < 1e-6, the modifications at x=1e-6 are returned,
//
// The data used by the program for required order
// and atomic number A, are stored in separate files
//
//   LO : EPS09LOR_A
//   NLO: EPS09NLOR_A
//
// which must be located in the working directory.
//
// The error bands for absolute cross-sections and for
// their nuclear ratios should be computed as explained
// in Secs. 2.5 and 4 of arXiv:0902.4154 [hep-ph]. For
// the absolute cross sections, both the errors in the
// free-proton PDFs f_p(x,Q) and the errors in
// the modifications R_f^A(x,Q) should be accounted for.
// For the nuclear ratios, it is sufficient to account only
// for the errors in the modifications R_f^A(x,Q).
//
//=================================================================================
#include <iostream>
#include <fstream>
#include <cmath>
#include "eps09.h"

#include <cstdlib>

using namespace std;
    
void eps09(int order, int pset, int AAA,
           double xxx, double QQQ, double& ruv,
           double& rdv, double& ru, double& rd,
           double& rs, double& rc, double& rb,
           double& rg) {

    const double Q2min = 1.69;
    const double Q2max = 1000000;
    const double Qsteps = 50;
    
    double LSTEP, x, Q, Q2;
    double x_i=0.000001, arg[4+1], fu[4+1], res, fg[3+1];
    double result[9+1], dummy;
    double realQ, n_x;

    char filenimi[50];

    int xlinsteps=25, xlogsteps=25;
    int k, p, t, Qpoint, xpoint;
    int setnumber, j, A;

    static double allvalues[31+1][8+1][50+1][50+1];
    //static int psetlast = -10;
    static int Alast = -10, orderlast = -10;

    //
    // Stop if the set specifications are wrong ones
    //

    if (order != 1 && order != 2) {
        cerr << "Wrong order!" << endl;
        cerr << "LO : order = 1" << endl;
        cerr << "NLO: order = 2" << endl;
        exit(1);
    }

    if (pset  < 1 || pset > 31) {
        cerr << "Wrong set!" << endl;
        cerr << "Central set: pset = 1" << endl;
        cerr << "Error sets : pset = 2...31" << endl;
        exit(1);
    }

    //
    // Make sure not to change any
    // specifications given by the user
    //

    A  = AAA;
    x  = xxx;
    Q  = QQQ;
    Q2 = Q*Q; 

    //
    // Freeze x if it's < 10E-6 or > 1
    //

    if (x < x_i) x = x_i;
    if (x > 1) x = 1;

    //
    // Freeze Q^2 if it's < 1.69 or > 10E+6
    //

    if (Q2 < Q2min) Q2 = Q2min;     
    if (Q2 > Q2max) Q2 = Q2max;

    //
    // If the set specifications have been changed, read the tables again
    //

    if (A != Alast || order != orderlast) {
        //
        // Read the table
        //
        if (order == 1) 
            snprintf(filenimi, 50, "./eps09/EPS09LOR_%d", A);
        else
            snprintf(filenimi, 50, "./eps09/EPS09NLOR_%d", A);

        ifstream ifs(filenimi);

        if (!ifs) {
            cerr << "Missing file: " << filenimi << endl;
            exit(1);
        }

        for (setnumber = 1; setnumber <=31; setnumber++) {
            for (k = 0; k <=50; k++) {               
                ifs >> dummy;
                for (t = 0; t <=50; t++)
                    for (p=1; p<=8; p++) ifs >> allvalues[setnumber][p][k][t];
            }   
        }

        ifs.close();

        //psetlast  = pset;
        Alast     = A;
        orderlast = order;

    }

    //
    // Find out the position in the loglog Q^2-grid
    //

    realQ  = Qsteps * (log(log(Q2)/log(Q2min)))/
                      (log(log(Q2max)/log(Q2min)));
    Qpoint = static_cast<int>(realQ);

    if (Qpoint <= 0) Qpoint = 1;

    if (Qpoint >= static_cast<int>(Qsteps+0.5)-1)
        Qpoint = static_cast<int>(Qsteps+0.5)-1;

    LSTEP = (1.0/(xlogsteps)) * log(0.1/x_i);

    //
    // Interpolate the grids 
    //

    for (t=1; t<=8; t++) {

        // Find the position in the x-grid

        if (x <= 0.1) 
            n_x  = ((1.0/LSTEP) * log(x/x_i));
        else
            n_x    = ((x-0.1)*xlinsteps/(1.0-0.1) + xlogsteps);
        xpoint = static_cast<int>(n_x);
        
        if (xpoint <= 0) xpoint = 1;
        
        if (t == 1 || t == 2) 
            if (xpoint >= (xlinsteps+xlogsteps)-4)
                xpoint = (xlinsteps+xlogsteps)-4;
        
        if (t == 3 || t == 4  || t == 5  || t == 6 || t == 7) 
            if (xpoint >= (xlinsteps+xlogsteps)-7)
                xpoint = (xlinsteps+xlogsteps)-7;
        
        if (t == 8)
            if (xpoint >= (xlinsteps+xlogsteps)-4)
                xpoint = (xlinsteps+xlogsteps)-4;
        
        for (k= 1; k<=4; k++) {
            if (xpoint-2+k < xlogsteps)
                arg[k] = (x_i) * exp(LSTEP * (xpoint-2+k));
            else
                arg[k] = 0.1 + (xpoint-2+k-xlogsteps) * (1-0.1)/xlinsteps;
        }
        
        for (j= 1; j<=3; j++) {
            fu[1] = allvalues[pset][t][Qpoint-2+j][xpoint-1];
            fu[2] = allvalues[pset][t][Qpoint-2+j][xpoint];
            fu[3] = allvalues[pset][t][Qpoint-2+j][xpoint+1];
            fu[4] = allvalues[pset][t][Qpoint-2+j][xpoint+2];
            luovi(fu,arg,4,x,res);
            fg[j] = res;
        }
        
        arg[1] = Qpoint-1;
        arg[2] = Qpoint;
        arg[3] = Qpoint+1;

        luovi(fg,arg,3,realQ,res);
        
        result[t] = res;
    }
    
    ruv = result[1] > 0 ? result[1] : 0;
    rdv = result[2] > 0 ? result[2] : 0;
    ru  = result[3] > 0 ? result[3] : 0;
    rd  = result[4] > 0 ? result[4] : 0;
    rs  = result[5] > 0 ? result[5] : 0;
    rc  = result[6] > 0 ? result[6] : 0;
    rb  = result[7] > 0 ? result[7] : 0;
    rg  = result[8] > 0 ? result[8] : 0;
}

//
// Modified version of Cern Library
// interpolation routine E100
//
void luovi(const double* f, const double* arg, int mmm, double z, double &sum)
{
    double *cof = new double[mmm+1];
    int jndex, index;
    
    int mm = mmm < 20 ? mmm : 20;
    int m = mm - 1;
    for (int i=1; i<=mm; i++) cof[i] = f[i];
    for (int i=1; i<=m; i++) {
        for (int j=i; j<=m; j++){
            jndex = mm - j;
            index = jndex + i;
            cof[index] = (cof[index]-cof[index-1])/(arg[index]-arg[jndex]);
        }
    }
    sum = cof[mm];
    for (int i=1; i<=m; i++){
        index = mm - i;
        sum = (z-arg[index])*sum + cof[index];
    }
    delete [] cof;
}
