// This file is part of Resummino.
//
// Copyright 2008-2011 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Utils for other modules.

#ifndef UTILS_H_
#define UTILS_H_

#define RESUMMINO_VERSION "1.0.7"

#include <complex>
#include <cstdio>

using namespace std;

// Returns x^2.
static inline double pow2(double x)
{
    return x * x;
}

// Returns x^2.
static inline complex<double> pow2(complex<double> x)
{
    return x * x;
}

// Returns x^4.
static inline double pow4(double x)
{
    double y = pow2(x);
    return y * y;
}

static inline int  pdg_particle_charge(int pdg_code) {
  switch (pdg_code) {
  case 1000012: // ~nu_eL
  case 1000014: // ~nu_muL
  case 1000016: // ~nu_tau
  case 12: // nu_eL
  case 14: // nu_muL
  case 16: // nu_tauL
  case 1000022: // ~chi_10
  case 1000023: // ~chi_20
  case 1000025: // ~chi_30
  case 1000035: // ~chi_40
  case 1000021: // ~g
  case -1000012: // ~nu_eL
  case -1000014: // ~nu_muL
  case -1000016: // ~nu_tau
  case -12: // nu_eL
  case -14: // nu_muL
  case -16: // nu_tauL
  case -1000022: // ~chi_10
  case -1000023: // ~chi_20
  case -1000025: // ~chi_30
  case -1000035: // ~chi_40
  case -1000021: // ~g
    return 0;
  case 1000011: // ~e_L
  case 2000011: // ~e_R
  case 1000013: // ~mu_L
  case 2000013: // ~mu_R
  case 1000015: // ~tau_1
  case 2000015: // ~tau_2
  case 11: // e
  case 13: // mu
  case 15: // tau
  case -1000024: // ~chi_1-
  case -1000037: // ~chi_2-
  case -1000002: // ~u_L*
  case -2000002: // ~u_R*
  case -1000004: // ~c_L*
  case -2000004: // ~c_R*
  case -1000006: // ~t_1*
  case -2000006: // ~t_2*
  case 1000001: // ~d_L
  case 2000001: // ~d_R
  case 1000003: // ~s_L
  case 2000003: // ~s_R
  case 1000005: // ~b_1
  case 2000005: // ~b_2
    return -1;
  case -1000011: // ~e_L
  case -2000011: // ~e_R
  case -1000013: // ~mu_L
  case -2000013: // ~mu_R
  case -1000015: // ~tau_1
  case -2000015: // ~tau_2
  case -11: // e
  case -13: // mu
  case -15: // tau
  case 1000024: // ~chi_1-
  case 1000037: // ~chi_2-
  case 1000002: // ~u_L*
  case 2000002: // ~u_R*
  case 1000004: // ~c_L*
  case 2000004: // ~c_R*
  case 1000006: // ~t_1*
  case 2000006: // ~t_2*
  case -1000001: // ~d_L
  case -2000001: // ~d_R
  case -1000003: // ~s_L
  case -2000003: // ~s_R
  case -1000005: // ~b_1
  case -2000005: // ~b_2
    return 1;
  default:
    fprintf(stderr, "error: Wrong PDG code.\n");
    exit(1);
  }
}

#endif
