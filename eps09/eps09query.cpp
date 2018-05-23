//=================================================================================
//
//	File: eps09query.cxx
//  Author: Thomas Ullrich (thomas.ullrich@bnl.gov)
//  Last modified: June 16, 2009
//
//  Command line tool to query values out of EPS09.
//  Syntax:
//
//  eps09query order pset A x Q
//
//  For details on arguments see eps09.cxx
//
//=================================================================================
#include <iostream>
#include <stdlib.h>
#include "eps09.h"
#define PR(x) cout << #x << " = " << (x) << endl;
using namespace std;

int main(int argc, char **argv) 
{
    if (argc != 6) {
        cerr << "Usage: " << argv[0] << " order pset A x Q" << endl;
        return 2;
    }
    
    int order = atoi(argv[1]);
    int pset = atoi(argv[2]);
    int A = atoi(argv[3]);
    double x = atof(argv[4]);
    double Q = atof(argv[5]);
    double ruv, rdv, ru, rd, rs, rc, rb, rg;
    
    eps09(order, pset, A, x, Q, 
          ruv, rdv, ru, rd, rs, rc, rb, rg);
    
    PR(ruv);
    PR(rdv);
    PR(ru);
    PR(rd);
    PR(rs);
    PR(rc);
    PR(rb);
    PR(rg);
    
    return 0;
}
