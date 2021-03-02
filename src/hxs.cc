// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2015 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "params.h"
#include "dipoles.h"
#include "gsl_all.h"
#include "integration_method.h"
#include "pdf.h"
#include "pxs.h"
#include "utils.h"

// Minimum pt cut
// (useful to check e.g. against madgraph; or to check dipoles)
// (use a pt cut and then check separately the DIPOLE and REAL close to pt->0;
// you should get the same results with a different sign)
// #define PTMINCUT
// #define PTMIN 0.001

// Maximum pt cut
// #define PTMAX 0.01
// #define PTMAXCUT

#define DIPOLE // activate or deactivate the 2->3 dipole
#define REAL   // activate or deactivate the real emission

// turn off if you want to neglect the resonant TT and UU diagrams
// for associated gaugino-gluino production.
// (useful to just check the IR divergent contributions).
#define ONSUB

using namespace std;

// two-particle phase space used for LO and virtual corrections
double dPS2(double &xa, double &xb, double &t, double *x, Parameters *params) {
  const double sh = params->sh;

  double m1;
  double m2;

  SET_MASS;

  const double m1s = pow2(m1);
  const double m2s = pow2(m2);

  // Conversion factor from GeV^-2 to pb.
  double djac = 389379304.0;

  // Integration variables xa and xb.
  double xamin = pow2(m1 + m2) / sh;
  double xamax = 1.0;
  double xbmin = xamin;
  double xbmax = 1.0;

  // for massless final state particles like leptons
  if (is_lepton_lepton(params->out1, params->out2)) {
    xamin = pow2(params->Minv_min) / sh;
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
    xbmin = xamin / xa;
    xb = (xbmax - xbmin) * x[1] + xbmin;
    
    if ((xa * xb * sh < pow2(params->Minv_min)) || (xa * xb * sh > pow2(params->Minv_max))) {
        djac = 0.0;
        return 0.0;
    }
    
    djac *= xbmax - xbmin;

    // We can only integrate from 0 to 1, so we substitute
    // int_a^b f(x) dx = int_0^1 f[a + (b-a)y] (b-a) y.
  } else {
    if (xamin == 0.0) {
      xa = (xamax - xamin) * x[0] + xamin;
      djac *= xamax - xamin;
      xbmin /= xa;
      xb = (xbmax - xbmin) * x[1] + xbmin;
      djac *= xbmax - xbmin;

      // Another substitution.
    } else {
      xa = xamin * pow(xamax / xamin, x[0]);
      djac *= xa * log(xamax / xamin);
      xbmin /= xa;
      xb = xbmin * pow(xbmax / xbmin, x[1]);
      djac *= xb * log(xbmax / xbmin);
    }
  }

  // the partonic com energy
  const double s = xa * xb * sh;

  // Integration variable t. (kln is sqrt of the kaellen function!)
  const double tmin = -.5 * (s - m1s - m2s + kln(s, m1s, m2s));
  const double tmax = tmin + kln(s, m1s, m2s);
  t = (tmax - tmin) * x[2] + tmin;
  djac *= tmax - tmin;

  // Phase Space factor.
  djac *= 0.125 * M_1_PI / s;

  return djac;
}

// Phase space for real collinear emission.
double dPS2c(double &djacdelta, double &djacplus, double &xa, double &xb,
             double &xc, double &t, double &tc, double *x, Parameters *params) {
  const double sh = params->sh;

  double m1;
  double m2;

  // macro to set masses m1 and m2 for the different processes
  // (see utils.h)
  SET_MASS;

  const double m1s = pow2(m1);
  const double m2s = pow2(m2);

  // Conversion factor from GeV^-2 to pb.
  double djac = 389379304.0;

  // Integration variables xa and xb.
  double xamin = pow2(m1 + m2) / sh;
  double xamax = 1.0;
  double xbmin = xamin;
  double xbmax = 1.0;
  double xcmin = xbmin; // -> xc_min = (m1 + m2)^2 / s := z ; needed for correct
                        // use of plus-distribution
  double xcmax = 1.0;

  // mapping to interval [0:1]
  if (is_lepton_lepton(params->out1, params->out2)) {
    xamin = pow2(params->Minv_min) / sh;
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
    xbmin = xamin / xa;
    xb = (xbmax - xbmin) * x[1] + xbmin;
    
    if ((xa * xb * sh < pow2(params->Minv_min)) || (xa * xb * sh > pow2(params->Minv_max))) {
        djac = 0.0;
        djacdelta = 0.0;
        djacplus = 0.0;
        return 0;
    }
    
    djac *= xbmax - xbmin;
    xcmin = xamin / (xa * xb);
    xc = (xcmax - xcmin) * x[2] + xcmin;
    djac *= xcmax - xcmin;
  } else {

    if (xamin == 0.0) {
      xa = (xamax - xamin) * x[0] + xamin;
      djac *= xamax - xamin;
      xbmin /= xa;
      xb = (xbmax - xbmin) * x[1] + xbmin;
      djac *= xbmax - xbmin;
      djacdelta = djac; // note that djacdelta is independent of x[2], the
                        // actual xc integration
      xcmin /= xa * xb;
      xc = (xcmax - xcmin) * x[2] + xcmin;
      djac *= xcmax - xcmin;

    } else {
      xa = xamin * pow(xamax / xamin, x[0]);
      djac *= xa * log(xamax / xamin);
      xbmin /= xa;
      xb = xbmin * pow(xbmax / xbmin, x[1]);
      djac *= xb * log(xbmax / xbmin);
      djacdelta = djac;
      xcmin /= xa * xb;
      xc = xcmin * pow(xcmax / xcmin, x[2]);
      djac *= xc * log(xcmax / xcmin);
    }
  }

  djacplus = djac;
  const double s = xa * xb * sh;

  // Integration variable t. (for the delta and plus distribution part)
  const double tmin = -0.5 * (s - m1s - m2s + kln(s, m1s, m2s));
  const double tmax = tmin + kln(s, m1s, m2s);
  t = (tmax - tmin) * x[3] + tmin;
  djacdelta *= tmax - tmin;
  djacplus *= tmax - tmin;

  const double sc = xc * s;

  // Integration variable tc. (for the integration of the collinear phase space)
  const double tcmin = -0.5 * (sc - m1s - m2s + kln(sc, m1s, m2s));
  const double tcmax = tcmin + kln(sc, m1s, m2s);
  tc = (tcmax - tcmin) * x[3] + tcmin;
  djac *= tcmax - tcmin;

  // Phase Space factor.
  // Here we have three different phase space factors.
  // djacdelta is for the integrand which is proportional to delta(1-xc),
  // djacplus for the subtraction term of the plus distribution
  // and djac for the collinear phase space.
  // Note: the delta(1-xc) part is independent of xc, thus, integrating it over
  // x[2]~xc
  // from 0 to 1 gives just a factor 1
  djacdelta *= 0.125 * M_1_PI / s;
  djacplus *= 0.125 * M_1_PI / s;
  djac *= 0.125 * M_1_PI / sc;

  return djac;
}

// Three-particle phase space.
double dPS3(double &xa, double &xb, double &M2, double &pt2, double &Tp,
            double &Pp, double *x, Parameters *params) {
  const double sh = params->sh;

  double m1;
  double m2;

  SET_MASS;

  const double m1s = pow2(m1);
  const double m2s = pow2(m2);

  // Conversion factor from GeV^-2 to pb.
  double djac = 389379304.0;

  // Integration variable xa, xb and M2
  double xamin = pow2(m1 + m2) / sh;
  double xamax = 1.0;
  double xbmin = xamin;
  double xbmax = 1.0;
  double M2min = xamin * sh;
  double M2max = sh;

  // Massless outgoing particles
  if (is_lepton_lepton(params->out1, params->out2)) {
    xamin = pow2(params->Minv_min) / sh;
    xbmin = xamin;
    M2min = pow2(params->Minv_min);
    M2max = pow2(params->Minv_max);
    
    if ((xa * xb * sh < pow2(params->Minv_min)) || (xa * xb * sh > pow2(params->Minv_max))) {
        djac = 0.0;
        return 0;
    }
    
  }
  if (xamin == 0.0) {
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
    xbmin /= xa;
    xb = (xbmax - xbmin) * x[1] + xbmin;
    djac *= xbmax - xbmin;
    M2max *= xa * xb;
    M2 = (M2max - M2min) * x[2] + M2min;
    djac *= M2max - M2min;
  } else {
    xa = xamin * pow(xamax / xamin, x[0]);
    djac *= xa * log(xamax / xamin);
    xbmin /= xa;
    xb = xbmin * pow(xbmax / xbmin, x[1]);
    djac *= xb * log(xbmax / xbmin);
    M2max *= xa * xb;
    M2 = M2min * pow(M2max / M2min, x[2]);
    djac *= M2 * log(M2max / M2min);
  }

  const double s = xa * xb * sh;

  double pt2min = 0.0;
  double pt2max = 0.25 * pow2(s - M2) / s;

// some used preprocessor variables to add a pt cutoff
// we set pt2 with our functions (&pt2)
#ifdef PTMAXCUT
  if (pt2max > PTMAX * PTMAX) {
    pt2max = PTMAX * PTMAX;
  }
#endif
#ifdef PTMINCUT
  pt2min = PTMIN * PTMIN;
#endif

  // set the pt2 value; the check is performed for the cutoff.
  pt2 = (pt2max - pt2min) * x[3] + pt2min;
  if (pt2min > pt2max || pt2 > pt2max || pt2 < pt2min) {
    djac = 0.0;
    return 0;
  }

  djac *= pt2max - pt2min;

  // Integration variable Tp.
  const double Tpmin = 0.0;
  const double Tpmax = M_PI;
  Tp = (Tpmax - Tpmin) * x[4] + Tpmin;
  djac *= Tpmax - Tpmin;

  // Integration variable Pp.
  const double Ppmin = 0.0;
  const double Ppmax = M_PI;
  djac *= 2.0;
  Pp = (Ppmax - Ppmin) * x[5] + Ppmin;
  djac *= Ppmax - Ppmin;

  // Phase Space factor
  djac *= kln(M2, m1s, m2s) * sin(Tp) / sqrt(pow2(s - M2) - 4.0 * s * pt2) /
          M2 / pow4(4.0 * M_PI);

  return djac;
}

// Three-particle phase space used for on-shell subtraction for the
// associated production of gauginos and gluinos.
// Here s2 = s23 can be on-shell. s1 is the usual M^2.
// (Reference: Particle kinematics - Byckling, Kajantie)
// This phase space is only used for real quark emission where
// squarks can be produced on-shell and need to be subtracted.
// This is the so called helicity-frame.
double dPS3_ONSHELL23(double &xa, double &xb, double &s2, double &t1,
                      double &s1, double &phi, double *x, Parameters *params) {

  // hadronic com energy
  const double sh = params->sh;

  // final state masses
  double m1;
  double m2;

  // set final state masses with macro
  SET_MASS;

  const double m1s = pow2(m1);
  const double m2s = pow2(m2);

  // initialize jacobian
  double djac = 1.0;
  // Conversion factor from GeV^-2 to pb.
  djac = djac * 389379304.0;

  // Integration variable xa, xb
  double xamin = pow2(m1 + m2) / sh;
  double xamax = 1.0;
  double xbmin = xamin;
  double xbmax = 1.0;

  if (is_lepton_lepton(params->out1, params->out2)) {
    xamin = pow2(params->Minv_min) / sh;
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
    xbmin = xamin / xa;
    xb = (xbmax - xbmin) * x[1] + xbmin;
    
    if ((xa * xb * sh < pow2(params->Minv_min)) || (xa * xb * sh > pow2(params->Minv_max))) {
        djac = 0.0;
        return 0;
    }
    djac *= xbmax - xbmin;
  } else if (xamin == 0.0) {
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
    xbmin /= xa;
    xb = (xbmax - xbmin) * x[1] + xbmin;
    djac = djac * (xbmax - xbmin);
  } else {
    xa = xamin * pow(xamax / xamin, x[0]);
    djac *= xa * log(xamax / xamin);
    xbmin /= xa;
    xb = xbmin * pow(xbmax / xbmin, x[1]);
    djac = djac * xb * log(xbmax / xbmin);
  }

  // partonic com energy
  const double s = xa * xb * sh;

  // Integration variable s2 = s23 = 2.0 p2p3 + m2s;
  double s2min = m2s;
  double s2max = pow2(sqrt(s) - m1);

  // different mappings of s2 to [0:1] (linear, log, breit)
  // if (params->deg_squarks == false) {
  //   double z4p = atan((s2max - params->mSQs[0])/(params->mSQs[0] * 1.0E-2));
  //   double z4m = atan((s2min - params->mSQs[0])/(params->mSQs[0] * 1.0E-2));
  //   double y = (z4p - z4m) * x[2] + z4m;
  //   s2 = params->mSQs[0] + params->mSQs[0] * 1.0E-2 * tan(y);
  //   djac = djac * params->mSQs[0] * 1.0E-2 * (z4p - z4m)
  //     * 1.0/pow2(cos(z4m - x[2] * z4m + x[2] * z4p));
  // } else {
  if (s2max == 0.0 || s2min == 0.0) { // actually can never happen :D
    // linear maping
    s2 = (s2max - s2min) * x[2] + s2min;
    djac = djac * (s2max - s2min);
  } else {
    // log mapping
    s2 = s2min * pow(s2max / s2min, x[2]);
    djac *= s2 * log(s2max / s2min);
  }
  //    }

  // Integration variable t1 = (pa - p2)^2
  double t1min;
  double t1max;
  t1min = -.5 * (s - m1s - s2 + kln(s, s2, m1s));
  t1max = t1min + kln(s, s2, m1s);
  t1 = (t1max - t1min) * x[3] + t1min;
  djac = djac * (t1max - t1min);

  // Integration over s1 = s12 = M2 = (p1 + p2)^2.
  double s1min = m1s + m2s + 1.0 / (2.0 * s2) * (s - s2 - m1s) * (s2 + m2s) -
                 1.0 / (2.0 * s2) * (kln(s, s2, m1s) * kln(s2, m2s, 0));
  double s1max = m1s + m2s + 1.0 / (2.0 * s2) * (s - s2 - m1s) * (s2 + m2s) +
                 1.0 / (2.0 * s2) * (kln(s, s2, m1s) * kln(s2, m2s, 0));
  s1 = (s1max - s1min) * x[4] + s1min;
  djac = djac * (s1max - s1min);

  // Integration over azimuthal angle.
  double phimin = 0.0;
  double phimax = 2.0 * M_PI;

  phi = (phimax - phimin) * x[5] + phimin;
  djac = djac * (phimax - phimin);

  // Normalization factor.
  djac = djac / pow((2.0 * M_PI), 5);

  double sqrt_lambda = s;

  // Final Jacobian.
  djac = djac * M_PI / (2.0 * sqrt_lambda) / (4.0 * kln(s, s2, m1s));

  // Checks if phase space point is physical.
  if (Rkln(s, s2, m1s) < 0.0 || Rkln(s2, m2s, 0) < 0.0) {
    djac = 0.0;
  }

  return djac;
}

// Three-particle phase space used for on-shell subtraction for the
// associated production of gauginos and gluinos.
// Here s2 = s13 can be on-shell. s1 = s12 is the usual M2.
// Similar to dPS3_ONSHELL23, but p1 <-> p2.
// this phase space is introduced to easily integrate the other resonant region;
// otherwise it is hidden in the angular integration!
double dPS3_ONSHELL13(double &xa, double &xb, double &s2, double &t1,
                      double &s1, double &phi, double *x, Parameters *params) {

  // hadronic com energy
  const double sh = params->sh;

  // final state masses
  double m1;
  double m2;

  // macro to set final state masses
  SET_MASS;

  const double m1s = pow2(m1);
  const double m2s = pow2(m2);

  double djac = 1.0;
  // Conversion factor from GeV^-2 to pb.
  djac = djac * 389379304.0;

  // Integration variable xa, xb
  double xamin = pow2(m1 + m2) / sh;
  double xamax = 1.0;
  double xbmin = xamin;
  double xbmax = 1.0;

  if (is_lepton_lepton(params->out1, params->out2)) {
    xamin = pow2(params->Minv_min) / sh;
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
    xbmin = xamin / xa;
    xb = (xbmax - xbmin) * x[1] + xbmin;
    
    if ((xa * xb * sh < pow2(params->Minv_min)) || (xa * xb * sh > pow2(params->Minv_max))) {
//         djac = 0.0;
//         return 0;
        return 0.0;
    }
    
    djac *= xbmax - xbmin;
  } else if (xamin == 0.0) {
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
    xbmin /= xa;
    xb = (xbmax - xbmin) * x[1] + xbmin;
    djac = djac * (xbmax - xbmin);
  } else {
    xa = xamin * pow(xamax / xamin, x[0]);
    djac *= xa * log(xamax / xamin);
    xbmin /= xa;
    xb = xbmin * pow(xbmax / xbmin, x[1]);
    djac = djac * xb * log(xbmax / xbmin);
  }

  const double s = xa * xb * sh;

  double s2min = m1s;
  double s2max = pow2(sqrt(s) - m2);

  // different mapping; breit wigner only of squarks are mass degenerate
  // if (params->deg_squarks == false) {
  //   double z4p = atan((s2max - params->mSQs[0])/(params->mSQs[0] * 1.0E-2));
  //   double z4m = atan((s2min - params->mSQs[0])/(params->mSQs[0] * 1.0E-2));
  //   double y = (z4p - z4m) * x[2] + z4m;
  //   s2 = params->mSQs[0] + params->mSQs[0] * 1.0E-2 * tan(y);
  //   djac = djac * params->mSQs[0] * 1.0E-2 * (z4p - z4m)
  //     * 1.0/pow2(cos(z4m - x[2] * z4m + x[2] * z4p));
  // } else {
  if (s2max == 0.0 || s2min == 0.0) { // actually can never happen :D
    s2 = (s2max - s2min) * x[2] + s2min;
    djac = djac * (s2max - s2min);
  } else {
    // log mapping
    s2 = s2min * pow(s2max / s2min, x[2]);
    djac *= s2 * log(s2max / s2min);
  }
  //}

  double t1min;
  double t1max;
  t1min = -.5 * (s - m2s - s2 + kln(s, s2, m2s));
  t1max = t1min + kln(s, s2, m2s);
  t1 = (t1max - t1min) * x[3] + t1min;
  djac = djac * (t1max - t1min);

  double s1min = m1s + m2s + 1.0 / (2.0 * s2) * (s - s2 - m2s) * (s2 + m1s) -
                 1.0 / (2.0 * s2) * (kln(s, s2, m2s) * kln(s2, m1s, 0));
  double s1max = m1s + m2s + 1.0 / (2.0 * s2) * (s - s2 - m2s) * (s2 + m1s) +
                 1.0 / (2.0 * s2) * (kln(s, s2, m2s) * kln(s2, m1s, 0));
  s1 = (s1max - s1min) * x[4] + s1min;
  djac = djac * (s1max - s1min);

  // Integration over azimuthal angle.
  double phimin = 0.0;
  double phimax = 2.0 * M_PI;

  phi = (phimax - phimin) * x[5] + phimin;
  djac = djac * (phimax - phimin);

  // Normalization factor.
  djac = djac / pow((2.0 * M_PI), 5);

  double sqrt_lambda = s;
  // Final Jacobian.
  djac = djac * M_PI / (2.0 * sqrt_lambda) / (4.0 * kln(s, s2, m2s));

  if (Rkln(s, s2, m2s) < 0.0 || Rkln(s2, m1s, 0) < 0.0) {
    djac = 0.0;
  }

  return djac;
}

// Setting up the different integrands.

// Born
double IB(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj; // conversion of type void* into IOS*

  // Setting the needed variables xa,xb,t and s for a specific integration point
  // x (array).
  double xa, xb, t;
  const double djac = dPS2(xa, xb, t, x, params);  
  const double s = xa * xb * params->sh;

  // Symmetry factor for two neutralinos in the final state
  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 &&
      params->out1 / 4 == 0) {
    dij = 0.5;
  }

  // PDFs of gluons and quarks for incoming particles a and b.
  double ga, qa[2][6];
  double gb, qb[2][6];

  // fill the PDF arrays
  pdfX(ga, qa, xa, params->mufs);
  pdfX(gb, qb, xb, params->mufs);

  // Initialize the cross section.
  double sig = 0.0;

  // Initial state quark and gluon:
  // for associated squark gaugino pair production
  if (is_squark_gaugino(params->out1, params->out2)) {
    if (params->out1 > 30) {
      for (int i0 = 0; i0 < 5; i0++) {
        int i1 = 0;
        params->in1 = i0;
        params->in2 = i1;

        sig += (qa[0][i0] * gb + qb[0][i0] * ga) * born_gasq(s, t, params);
      }
    } else if (params->out2 > 30) {
      for (int i1 = 0; i1 < 5; i1++) {
        int i0 = 0;
        params->in1 = i0;
        params->in2 = i1;

        sig += (qa[1][i1] * gb + qb[1][i1] * ga) * born_gasq(s, t, params);
      }
    }

    double g3s = std::norm(params->gqq[0][0].R);
    return 4.0 * M_PI * aS(params->murs, params->set) * sig * djac * 1.0 /
           (2.0 * s * 96.0);
  }

  // Initial state: Quark and antiquark.
  // for drell yan like processes and associated gaugino-gluino production.
  // Sums over all possible initial states.
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      // proton-proton -> ic = 0; proton-antiproton -> ic = 1;
      if (is_charge_conserved(params->in1, params->in2, params->out1,
                              params->out2)) {
        if (is_gaugino_gluino(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] +
                  qb[params->ic][i0] * qa[1][i1]) *
                 born_gagl(s, t, params);
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] +
                  qb[params->ic][i0] * qa[1][i1]) *
                 born_sleptons(s, t, params);
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] +
                  qb[params->ic][i0] * qa[1][i1]) *
                 born_gauginos(s, t, params);
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] +
                  qb[params->ic][i0] * qa[1][i1]) *
                 born_leptons(s, t, params);
        }
      }
    }
  }

  // Flux, symmetry factor, spin and color average.
  if (is_gaugino_gluino(params->out1, params->out2)) {
    // LO with strong coupling.
    double g3s = std::norm(params->gqq[0][0].R); // needed to replace the fixed strong coupling by a running aS
    return 4.0 * M_PI * aS(params->murs, params->set) * sig * djac * 1.0 /
           (2.0 * s * 36.0) * 1.0 / g3s;
  } else {
    // LO without strong coupling.
    // LO color factor ->  3; average -> 1/(3*3*2*2); and 1/(2s) for the flux.
    // -> 1/24 as overall factor.
    return dij * sig * djac / (24.0 * s);
  }
}

// Virtual corrections. (completely analog to IB, but different partonic cross
// section)
double IV(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj; // conversion of type void* into IOS*
  double xa, xb, t;
  const double djac = dPS2(xa, xb, t, x, params);
  const double s = xa * xb * params->sh;

  // symmetry factor for identical final state particles
  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 &&
      params->out1 / 4 == 0) {
    dij = 0.5;
  }

  // PDFs of the gluons and the quarks for incoming particles a and b.
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double sig = 0.0;

  // initial state quark and gluon
  // for squark gaugino
  if (is_squark_gaugino(params->out1, params->out2)) {
    if (params->out1 > 30) {
      for (int i0 = 0; i0 < 5; i0++) {
        int i1 = 0;
        params->in1 = i0;
        params->in2 = i1;
        sig += (qa[0][i0] * gb + qb[0][i0] * ga) *
               (Virt_gaugino_squark(s, t, params) + DipI_gasq(s, t, params));
      }
    } else if (params->out2 > 30) {
      for (int i1 = 0; i1 < 5; i1++) {
        int i0 = 0;
        params->in1 = i0;
        params->in2 = i1;
        sig += (qa[1][i1] * gb + qb[1][i1] * ga) *
               (Virt_gaugino_squark(s, t, params) + DipI_gasq(s, t, params));
      }
    }

    double g3s = std::norm(params->gqq[0][0].R);
    return 4.0 * M_PI * aS(params->murs, params->set) * 4.0 * M_PI *
           aS(params->murs, params->set) * sig * djac * 1.0 / (2.0 * s * 96.0);
  }

  // Sum over all possible initial states.
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      params->in1 = i0;
      params->in2 = i1;
      if (is_charge_conserved(params->in1, params->in2, params->out1,
                              params->out2)) {
        if (is_gaugino_gluino(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] +
                  qb[params->ic][i0] * qa[1][i1]) *
                 (Virt_gaugino_gluino(s, t, params) + DipI_gagl(s, t, params));
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] +
                  qb[params->ic][i0] * qa[1][i1]) *
                 (Virt_sleptons(s, t, params) + DipI_sleptons(s, t, params));
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] +
                  qb[params->ic][i0] * qa[1][i1]) *
                 (Virt_gauginos(s, t, params) + DipI_gauginos(s, t, params));
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] +
                  qb[params->ic][i0] * qa[1][i1]) *
                 (Virt_leptons(s, t, params) + DipI_leptons(s, t, params));
        }
      }
    }
  }

  // Flux, symmetry factor, spin and color average.
  if (is_gaugino_gluino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    // pi^2 factors etc already included in pxs file.
    return (4.0 * M_PI * aS(params->murs, params->set)) *
           (4.0 * M_PI * aS(params->murs, params->set)) * dij * sig * djac *
           1.0 / (2.0 * s) * (1.0 / 36.0);
  } else {
    // Each virtual correction consists of a factor 1/(16 pi^2) due to the
    // scalar integral replacement.
    // Together with the factor 2 from 2*real(Mv * Mb) we get 2*g3s/16 pi^2 =
    // as/ 2*pi.
    return aS(params->murs, params->set) / (2.0 * M_PI) * dij * sig * djac *
           4.0 / (24.0 * s * 3.0);
  }
}

// Collinear emission.
double IC(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj;
  double djacdelta, djacplus, xa, xb, xc, t, tc;
  double djac = dPS2c(djacdelta, djacplus, xa, xb, xc, t, tc, x, params);
  const double s = xa * xb * params->sh;
  const double sc = xc * s;
  double z;

  // 2* pap1
  double sa1 = sja(s, t, params);
  // 2 * pbp1
  double sb1 = sjb(s, t, params);

  // 2 * pap1, where p1 belongs to the boosted phase space
  double sa1c = sja(sc, tc, params) / xc;
  // 2 * pbp1, where p1 belongs to the boosted phase space
  double sb1c = sjb(sc, tc, params) / xc;

  // This is z = (m1 + m2)^2/s (the actual lower integration limit of xc (see
  // above)).
  // In order to use the definition of the plus distribution we have to subtract
  // a finite z-dependent term.
  if (params->out1 >= 20 && params->out1 < 30) {
    z = pow2(params->Minv_min) / s;
  } else if (params->out1 >= 10 & params->out1 < 20) {
    z = pow2(params->mSL[params->out1 - 10] + params->mSL[params->out2 - 10]) /
        s;
  } else if (params->out1 == 30 & params->out2 < 10) {
    z = pow2(params->mGL + params->mCH[params->out2]) / s;
  } else if (params->out2 == 30 & params->out1 < 10) {
    z = pow2(params->mGL + params->mCH[params->out1]) / s;
  } else if (params->out1 > 30 & params->out2 < 10) {
    z = pow2(params->mSQ[params->out1 - 31] + params->mCH[params->out2]) / s;
  } else if (params->out2 > 30 & params->out1 < 10) {
    z = pow2(params->mSQ[params->out2 - 31] + params->mCH[params->out1]) / s;
  } else if (params->out1 < 10 & params->out2 < 10) {
    z = pow2(params->mCH[params->out1] + params->mCH[params->out2]) / s;
  }

  // Symmetry factor for two identical final state particles.
  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 &&
      params->out1 / 4 == 0) {
    dij = 0.5;
  }

  // PDFs
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  // initialize born (for delta and plus dist. and bornc)
  double born = 0.0;
  double bornc = 0.0;
  double sig = 0.0;

  // Sums over all possible initial states.
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      params->in1 = i0;
      params->in2 = i1;

      // set born and bornc for the different processes
      if (is_charge_conserved(params->in1, params->in2, params->out1,
                              params->out2)) {
        if (is_gaugino_gluino(params->out1, params->out2)) {
          born = born_gagl(s, t, params);
          bornc = born_gagl(sc, tc, params);
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          born = born_sleptons(s, t, params);
          bornc = born_sleptons(sc, tc, params);
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          born = born_gauginos(s, t, params);
          bornc = born_gauginos(sc, tc, params);
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          born = born_leptons(s, t, params);
          bornc = born_leptons(sc, tc, params);
        }

        // associated gaugino-gluino production
        // (treated differently due to massive colored final state gluino)
        if (is_gaugino_gluino(params->out1, params->out2)) {

          // Diagonal splitting.

          // P_bold_app_b_j is the expectation value of the insertion operator.
          // (it is factorization scale dependent)
          // The 1/xc is for the different flux factor -> 1/(2 sc) (see below).
          // quark with pa out of proton a and b.
          sig +=
              (qa[0][i0] * qb[1 - params->ic][i1] +
               qb[params->ic][i0] * qa[1][i1]) *

              // first two arguments are a, ap; then the two "spectators"
              // P_qq
              (P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                              PARTON_GLUINO, FUNCTION_ALL_WDELTA, xc, z, sa1c,
                              s, params->mufs, 1.0, INITIAL_AND_FINAL) *
                   djac * bornc * 1 / xc

               // subtraction term due to plus distribution
               -
               P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                              PARTON_GLUINO, FUNCTION_PLUS, xc, z, sa1, s,
                              params->mufs, 1.0, INITIAL_AND_FINAL) *
                   djacplus * born

               // parts proportional to the delta distribution
               +
               P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                              PARTON_GLUINO, FUNCTION_DELTA, xc, z, sa1, s,
                              params->mufs, 1.0, INITIAL_AND_FINAL) *
                   djacdelta * born);

          // antiquark with pb out of proton a and b
          // and P_qbqb (same expression as P_qq, but sa1->sb1)
          sig +=
              (qa[0][i0] * qb[1 - params->ic][i1] +
               qb[params->ic][i0] * qa[1][i1]) *

              (P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                              PARTON_GLUINO, FUNCTION_ALL_WDELTA, xc, z, sb1c,
                              s, params->mufs, 1.0, INITIAL_AND_FINAL) *
                   djac * bornc * 1 / xc

               -
               P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                              PARTON_GLUINO, FUNCTION_PLUS, xc, z, sb1, s,
                              params->mufs, 1.0, INITIAL_AND_FINAL) *
                   djacplus * born

               +
               P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                              PARTON_GLUINO, FUNCTION_DELTA, xc, z, sb1, s,
                              params->mufs, 1.0, INITIAL_AND_FINAL) *
                   djacdelta * born);

          // Similar as above, but now the expectation value of the K insertion
          // operator.
          // (not scale dependent)

          // quark with pa out of proton a and b
          sig += (qa[0][i0] * qb[1 - params->ic][i1] +
                  qb[params->ic][i0] * qa[1][i1]) *

                 (K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_GLUINO, FUNCTION_ALL_WDELTA, xc, z,
                                 params->mGL, sa1c, s, 1.0, INITIAL_AND_FINAL) *
                      djac * bornc * 1 / xc

                  -
                  K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_GLUINO, FUNCTION_PLUS, xc, z,
                                 params->mGL, sa1, s, 1.0, INITIAL_AND_FINAL) *
                      djacplus * born

                  +
                  K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_GLUINO, FUNCTION_DELTA, xc, z,
                                 params->mGL, sa1, s, 1.0, INITIAL_AND_FINAL) *
                      djacdelta * born);
          // antiquark with pb out of proton a and b
          sig += (qa[0][i0] * qb[1 - params->ic][i1] +
                  qb[params->ic][i0] * qa[1][i1]) *

                 (K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_GLUINO, FUNCTION_ALL_WDELTA, xc, z,
                                 params->mGL, sb1c, s, 1.0, INITIAL_AND_FINAL) *
                      djac * bornc * 1 / xc

                  -
                  K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_GLUINO, FUNCTION_PLUS, xc, z,
                                 params->mGL, sb1, s, 1.0, INITIAL_AND_FINAL) *
                      djacplus * born

                  +
                  K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_GLUINO, FUNCTION_DELTA, xc, z,
                                 params->mGL, sb1, s, 1.0, INITIAL_AND_FINAL) *
                      djacdelta * born);

          // Off-diagonal splitting.
          // no delta- or plus-distribution part!
          // gluon out of proton a (b) splits into an antiquark and a quark
          // going into the hard process

          // Pgq and Pgqb, sa1
          sig +=
              (qb[params->ic][i0] * ga + qb[1 - params->ic][i1] * ga) *

              (P_bold_aap_b_j(PARTON_GLUON, PARTON_QUARK, PARTON_ANTIQUARK,
                              PARTON_GLUINO, FUNCTION_ALL_WDELTA, xc, z, sa1c,
                              s, params->mufs, 1.0, INITIAL_AND_FINAL) *
               djac * bornc * 1 / xc

               );

          // Pgq and Pgqb, sb1
          sig +=
              (qa[0][i0] * gb + qa[1][i1] * gb) *

              (P_bold_aap_b_j(PARTON_GLUON, PARTON_QUARK, PARTON_ANTIQUARK,
                              PARTON_GLUINO, FUNCTION_ALL_WDELTA, xc, z, sb1c,
                              s, params->mufs, 1.0, INITIAL_AND_FINAL) *
               djac * bornc * 1 / xc

               );

          // Kgq, Kgqb and sa1
          sig += (qb[params->ic][i0] * ga + qb[1 - params->ic][i1] * ga) *

                 (K_bold_aap_b_j(PARTON_GLUON, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_GLUINO, FUNCTION_ALL_WDELTA, xc, z,
                                 params->mGL, sa1c, s, 1.0, INITIAL_AND_FINAL) *
                  djac * bornc * 1 / xc

                  );

          // Kgq, Kgqb and sb1
          sig += (qa[0][i0] * gb + qa[1][i1] * gb) *

                 (K_bold_aap_b_j(PARTON_GLUON, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_GLUINO, FUNCTION_ALL_WDELTA, xc, z,
                                 params->mGL, sb1c, s, 1.0, INITIAL_AND_FINAL) *
                  djac * bornc * 1 / xc

                  );

        }
        // cases without final state partons in LO
        else {

          // Collinear remainder for Drell-Yan-like processes.
          // (PARTON_NONE means no final state parton)
          // factor 2 due to symmetry (same P and K operators for gluon emission
          // from
          // the quark or antiquark)
          // bold Pqq

          // P_qq
          sig += 2 * (qa[0][i0] * qb[1 - params->ic][i1] +
                      qb[params->ic][i0] * qa[1][i1]) *
                 (P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_NONE, FUNCTION_ALL_WDELTA, xc, z, 0.0,
                                 s, params->mufs, 1.0, INITIAL) *
                      djac * bornc * 1 / xc -
                  P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_NONE, FUNCTION_PLUS, xc, z, 0.0, s,
                                 params->mufs, 1.0, INITIAL) *
                      djacplus * born +
                  P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_NONE, FUNCTION_DELTA, xc, z, 0.0, s,
                                 params->mufs, 1.0, INITIAL) *
                      djacdelta * born);

          // K_qq
          sig += 2 * (qa[0][i0] * qb[1 - params->ic][i1] +
                      qb[params->ic][i0] * qa[1][i1]) *
                 (K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_NONE, FUNCTION_ALL_WDELTA, xc, z, 0.0,
                                 0.0, s, 1.0, INITIAL) *
                      djac * bornc * 1 / xc -
                  K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_NONE, FUNCTION_PLUS, xc, z, 0.0, 0.0, s,
                                 1.0, INITIAL) *
                      djacplus * born +
                  K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_NONE, FUNCTION_DELTA, xc, z, 0.0, 0.0,
                                 s, 1.0, INITIAL) *
                      djacdelta * born);

          // Pgq and Pgqb
          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga + qa[1][i1] * gb +
                  qb[1 - params->ic][i1] * ga) *
                 (P_bold_aap_b_j(PARTON_GLUON, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_NONE, FUNCTION_ALL_WDELTA, xc, z, 0.0,
                                 s, params->mufs, 1.0, INITIAL) *
                  djac * bornc * 1 / xc);

          // Kgq and Kgqb
          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga + qa[1][i1] * gb +
                  qb[1 - params->ic][i1] * ga) *
                 (K_bold_aap_b_j(PARTON_GLUON, PARTON_QUARK, PARTON_ANTIQUARK,
                                 PARTON_NONE, FUNCTION_ALL_WDELTA, xc, z, 0.0,
                                 0.0, s, 1.0, INITIAL) *
                  djac * bornc * 1 / xc);
        }
      }
    }
  }

  if (is_gaugino_gluino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R); // LO with strong coupling
    return 1.0 / (36.0 * 2.0 * s) * sig * aS(params->murs, params->set) *
           (4.0 * M_PI * aS(params->murs, params->set)) / g3s;
  } else { // LO without strong coupling
    return dij * 3.0 / (36.0 * 2.0 * s) * sig * aS(params->murs, params->set);
  }
}

// Real gluon emission.
double IG(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj;
  double xa, xb, mi2, pt2, th, ph;
  const double djac = dPS3(xa, xb, mi2, pt2, th, ph, x, params);
  if (djac == 0) {
    return 0.0;
  }
  const double s = xa * xb * params->sh;

  // Numerical cutoff for stability.
  // (Some seeds could directly probe the divergent region and then you get a
  // nan - nan which is nan and not zero.)
  // Checked that for small pt values the dipoles completely reproduce the
  // 2->3 real emission diagrams.
  if (pt2 < 1.0E-6) {
    return 0.0;
  }

  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 &&
      params->out1 / 4 == 0) {
    dij = 0.5;
  }

  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double sig = 0.0;
  // Sums over all possible initial states.
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      params->in1 = i0;
      params->in2 = i1;

      if (is_charge_conserved(params->in1, params->in2, params->out1,
                              params->out2)) {

        // Real gluon emission for gaugino-gluino production.
        if (is_gaugino_gluino(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] +
                  qb[params->ic][i0] * qa[1][i1]) *
                 (
// Real emission partonic XS.
#ifdef REAL
                     real_gluon_gaugino_gluino(s, mi2, pt2, th, ph, 0, params) +
                     real_gluon_gaugino_gluino(s, mi2, pt2, th, ph, 1, params)
#endif
// Corresponding Catani-Seymour Dipoles dsigma^A.
// Initial state emitter (quark) and initial state spectator (antiquark)
#ifdef DIPOLE
                     - Dip_GLGA(INITIAL_INITIAL, PARTON_QUARK, PARTON_ANTIQUARK,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 0, params) -
                     Dip_GLGA(INITIAL_INITIAL, PARTON_QUARK, PARTON_ANTIQUARK,
                              PARTON_GLUON, s, mi2, pt2, th, ph, 1, params)

                     // Initial state emitter (antiquark) and initial state
                     // spectator (quark)
                     - Dip_GLGA(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 0, params) -
                     Dip_GLGA(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK,
                              PARTON_GLUON, s, mi2, pt2, th, ph, 1, params)

                     // Initial state emitter (quark) and final state spectator
                     // (gluino)
                     - Dip_GLGA(INITIAL_FINAL, PARTON_QUARK, PARTON_GLUINO,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 0, params) -
                     Dip_GLGA(INITIAL_FINAL, PARTON_QUARK, PARTON_GLUINO,
                              PARTON_GLUON, s, mi2, pt2, th, ph, 1, params)

                     // Initial state emitter (antiquark) and final state
                     // spectator (gluino)
                     - Dip_GLGA(INITIAL_FINAL, PARTON_ANTIQUARK, PARTON_GLUINO,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 0, params) -
                     Dip_GLGA(INITIAL_FINAL, PARTON_ANTIQUARK, PARTON_GLUINO,
                              PARTON_GLUON, s, mi2, pt2, th, ph, 1, params)

                     // Final state emitter (gluino) and initial state spectator
                     // (antiquark)
                     - Dip_GLGA(FINAL_INITIAL, PARTON_GLUINO, PARTON_ANTIQUARK,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 0, params) -
                     Dip_GLGA(FINAL_INITIAL, PARTON_GLUINO, PARTON_ANTIQUARK,
                              PARTON_GLUON, s, mi2, pt2, th, ph, 1, params)

                     // Final state emitter (gluino) and initial state spectator
                     // (quark)
                     - Dip_GLGA(FINAL_INITIAL, PARTON_GLUINO, PARTON_QUARK,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 0, params) -
                     Dip_GLGA(FINAL_INITIAL, PARTON_GLUINO, PARTON_QUARK,
                              PARTON_GLUON, s, mi2, pt2, th, ph, 1, params)
#endif
                         );
          // Real gluon emission for slepton pair production.
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          sig +=
              (qa[0][i0] * qb[1 - params->ic][i1] +
               qb[params->ic][i0] * qa[1][i1]) *
              (
#ifdef REAL
                  real_gluon_sleptons(s, mi2, pt2, th, ph, 0, params) +
                  real_gluon_sleptons(s, mi2, pt2, th, ph, 1, params)
#endif

#ifdef DIPOLE
                  -
                  (Dip_SLEPTONS(INITIAL_INITIAL, PARTON_QUARK, PARTON_ANTIQUARK,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 0, params) +
                   Dip_SLEPTONS(INITIAL_INITIAL, PARTON_QUARK, PARTON_ANTIQUARK,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 1, params)

                   +
                   Dip_SLEPTONS(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 0, params) +
                   Dip_SLEPTONS(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 1, params))
#endif
                      );
          // Real gluon emission for gaugino pair production.
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          sig +=
              (qa[0][i0] * qb[1 - params->ic][i1] +
               qb[params->ic][i0] * qa[1][i1]) *
              (
#ifdef REAL
                  real_gluon_gauginos(s, mi2, pt2, th, ph, 0, params) +
                  real_gluon_gauginos(s, mi2, pt2, th, ph, 1, params)
#endif

#ifdef DIPOLE
                  -
                  (Dip_GAUGINOS(INITIAL_INITIAL, PARTON_QUARK, PARTON_ANTIQUARK,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 0, params) +
                   Dip_GAUGINOS(INITIAL_INITIAL, PARTON_QUARK, PARTON_ANTIQUARK,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 1, params)

                   +
                   Dip_GAUGINOS(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 0, params) +
                   Dip_GAUGINOS(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK,
                                PARTON_GLUON, s, mi2, pt2, th, ph, 1, params))
#endif
                      );
          // Real gluon emission for lepton pair production.
          // (still using the "old" dipoles)
        }
        if (is_lepton_lepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] +
                  qb[params->ic][i0] * qa[1][i1]) *
                 (
#ifdef REAL
                     real_gluon_leptons(s, mi2, pt2, th, ph, 0, params) +
                     real_gluon_leptons(s, mi2, pt2, th, ph, 1, params)
#endif
#ifdef DIPOLE
                     - DipGA_leptons(s, mi2, pt2, th, ph, 0, params) -
                     DipGA_leptons(s, mi2, pt2, th, ph, 1, params) -
                     DipGB_leptons(s, mi2, pt2, th, ph, 0, params) -
                     DipGB_leptons(s, mi2, pt2, th, ph, 1, params)
#endif
                         );
        }
        if (isnan(sig)) {
          printf("pt2 = %.15f \n", pt2);
          exit(0);
        }
      }
    }
  }

  // real gluon emission for the associated squark gaugino production
  //(under construction)
  // squarks TODO: LASSE
  if (is_squark_gaugino(params->out1, params->out2)) {
    if (params->out1 > 30) {
      for (int i0 = 0; i0 < 5; i0++) {
        int i1 = 0;
        params->in1 = i0;
        params->in2 = i1;
        if (is_charge_conserved(params->in1, params->in2, params->out1,
                                params->out2)) {
          sig += (qa[0][i0] * gb + qb[0][i0] * ga) *
                 (real_gluon_gaugino_squark(s, mi2, pt2, th, ph, 0, params) +
                  real_gluon_gaugino_squark(s, mi2, pt2, th, ph, 1, params));
        }
      }
    }
  }

  // Flux, symmetry factor, spin and color average.
  if (is_gaugino_gluino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    // Factor 0.5 due to the phase space construction dPS3 with rapidity flags.
    return 0.5 * 1 / (72 * s) * (4 * M_PI * aS(params->murs, params->set)) *
           (4 * M_PI * aS(params->murs, params->set)) * sig * djac / g3s;
  } else {
    return dij * 0.5 * 4 / (72 * s) *
           (4 * M_PI * aS(params->murs, params->set)) * sig * djac;
  }
}

// Real quark and antiquark emission.
double IQ(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj; // conversion of type void* into IOS*
  double xa, xb, mi2, pt2, th, ph;
  const double djac = dPS3(xa, xb, mi2, pt2, th, ph, x, params);

  // Different three-particle phase space for on-shell subtraction
  // for the associated gaugino-gluino production.
  double s2_13, t1_13, s1_13, phi_13, djac13;
  double s2_23, t1_23, s1_23, phi_23, djac23;
  if (is_gaugino_gluino(params->out1, params->out2)) {
    djac13 = dPS3_ONSHELL13(xa, xb, s2_13, t1_13, s1_13, phi_13, x, params);
    djac23 = dPS3_ONSHELL23(xa, xb, s2_23, t1_23, s1_23, phi_23, x, params);
  }

  if (djac == 0.0) {
    return 0.0;
  }

#if !defined(DIPOLE) || !defined(REAL)
  djac13 = 0.0;
  djac23 = 0.0;
#endif

  const double s = xa * xb * params->sh;

  // Numerical cutoff for stability.
  // Due to the finite width and interference terms
  // of on- and off-shell contributions
  // non stable integration can occure
  // checked that dipole reproduces the real emission contribution
  // completely close to the divergent region
  if (pt2 < 1.0E-6) {
    return 0.0;
  }

#ifdef PTMINCUT
  if (pt2 < PTMIN * PTMIN) {
    return 0;
  }
#endif

#ifdef PTMAXCUT
  if (pt2 > PTMAX * PTMAX) {
    return 0;
  }
#endif

  // Symmetry factor.
  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 &&
      params->out1 / 4 == 0) {
    dij = 0.5;
  }

  // PDFs.
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double sig = 0.0;
  double sigOS = 0.0;

  // Color average for 2->3 and LO.
  double color_average_2to3 = 1.0 / 3.0 * 1.0 / 8.0;
  double color_average_born = 1.0 / 3.0 * 1.0 / 3.0;

  // Sums over all possible initial states.
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      params->in1 = i0;
      params->in2 = i1;

      // Tests charge conservation.
      if (is_charge_conserved(params->in1, params->in2, params->out1,
                              params->out2)) {

        // Associated gaugino-gluino production
        if (is_gaugino_gluino(params->out1, params->out2)) {

          // gluon with pb out of proton b and a emitting a quark.
          sig +=
              (qa[0][i0] * gb + qb[params->ic][i0] * ga) *

              (

// Partonic cross section for real quark emission.
#ifdef REAL
                  color_average_2to3 * (real_quark_gaugino_gluino(
                                            s, mi2, pt2, th, ph, 0, params) +
                                        real_quark_gaugino_gluino(
                                            s, mi2, pt2, th, ph, 1, params))
#endif

// Corresponding Dipoles.
#ifdef DIPOLE
                  -
                  color_average_born *
                      (Dip_GLGA(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK,
                                PARTON_QUARK, s, mi2, pt2, th, ph, 0, params) +
                       Dip_GLGA(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK,
                                PARTON_QUARK, s, mi2, pt2, th, ph, 1, params)

                       +
                       Dip_GLGA(INITIAL_FINAL, PARTON_GLUON, PARTON_GLUINO,
                                PARTON_QUARK, s, mi2, pt2, th, ph, 0, params) +
                       Dip_GLGA(INITIAL_FINAL, PARTON_GLUON, PARTON_GLUINO,
                                PARTON_QUARK, s, mi2, pt2, th, ph, 1, params))
#endif
                      );

          // gluon with pb out of proton b and a emitting an antiquark.
          sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                 (

// Partonic cross section for real antiquark emission.
#ifdef REAL
                     color_average_2to3 * (real_quarkb_gaugino_gluino(
                                               s, mi2, pt2, th, ph, 0, params) +
                                           real_quarkb_gaugino_gluino(
                                               s, mi2, pt2, th, ph, 1, params))
#endif

// Corresponding Dipoles.
#ifdef DIPOLE
                     -
                     color_average_born *
                         (Dip_GLGA(INITIAL_INITIAL, PARTON_GLUON,
                                   PARTON_ANTIQUARK, PARTON_ANTIQUARK, s, mi2,
                                   pt2, th, ph, 0, params) +
                          Dip_GLGA(INITIAL_INITIAL, PARTON_GLUON,
                                   PARTON_ANTIQUARK, PARTON_ANTIQUARK, s, mi2,
                                   pt2, th, ph, 1, params)

                          + Dip_GLGA(INITIAL_FINAL, PARTON_GLUON, PARTON_GLUINO,
                                     PARTON_ANTIQUARK, s, mi2, pt2, th, ph, 0,
                                     params) +
                          Dip_GLGA(INITIAL_FINAL, PARTON_GLUON, PARTON_GLUINO,
                                   PARTON_ANTIQUARK, s, mi2, pt2, th, ph, 1,
                                   params))
#endif

                         );

        } else if (is_slepton_slepton(params->out1, params->out2)) {

          sig +=
              (qa[0][i0] * gb + qb[params->ic][i0] * ga) *

              (
#ifdef REAL
                  real_quark_sleptons(s, mi2, pt2, th, ph, 0, params) +
                  real_quark_sleptons(s, mi2, pt2, th, ph, 1, params)
#endif

#ifdef DIPOLE
                  -
                  (Dip_SLEPTONS(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK,
                                PARTON_QUARK, s, mi2, pt2, th, ph, 0, params) +
                   Dip_SLEPTONS(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK,
                                PARTON_QUARK, s, mi2, pt2, th, ph, 1, params))
#endif
                      );

          sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *

                 (
#ifdef REAL
                     real_quarkb_sleptons(s, mi2, pt2, th, ph, 0, params) +
                     real_quarkb_sleptons(s, mi2, pt2, th, ph, 1, params)
#endif

#ifdef DIPOLE
                     - (Dip_SLEPTONS(INITIAL_INITIAL, PARTON_GLUON,
                                     PARTON_ANTIQUARK, PARTON_ANTIQUARK, s, mi2,
                                     pt2, th, ph, 0, params) +
                        Dip_SLEPTONS(INITIAL_INITIAL, PARTON_GLUON,
                                     PARTON_ANTIQUARK, PARTON_ANTIQUARK, s, mi2,
                                     pt2, th, ph, 1, params))
#endif
                         );

        } else if (is_gaugino_gaugino(params->out1, params->out2)) {

          sig +=
              (qa[0][i0] * gb + qb[params->ic][i0] * ga) *

              (
#ifdef REAL
                  real_quark_gauginos(s, mi2, pt2, th, ph, 0, params) +
                  real_quark_gauginos(s, mi2, pt2, th, ph, 1, params)
#endif
#ifdef DIPOLE
                  -
                  (Dip_GAUGINOS(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK,
                                PARTON_QUARK, s, mi2, pt2, th, ph, 0, params) +
                   Dip_GAUGINOS(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK,
                                PARTON_QUARK, s, mi2, pt2, th, ph, 1, params))
#endif
                      );

          sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *

                 (
#ifdef REAL
                     real_quarkb_gauginos(s, mi2, pt2, th, ph, 0, params) +
                     real_quarkb_gauginos(s, mi2, pt2, th, ph, 1, params)
#endif
#ifdef DIPOLE
                     - (Dip_GAUGINOS(INITIAL_INITIAL, PARTON_GLUON,
                                     PARTON_ANTIQUARK, PARTON_ANTIQUARK, s, mi2,
                                     pt2, th, ph, 0, params)

                        + Dip_GAUGINOS(INITIAL_INITIAL, PARTON_GLUON,
                                       PARTON_ANTIQUARK, PARTON_ANTIQUARK, s,
                                       mi2, pt2, th, ph, 1, params))
#endif
                         );

        } else if (is_lepton_lepton(params->out1, params->out2)) {

          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *

                 (
#ifdef REAL
                     real_quark_leptons(s, mi2, pt2, th, ph, 0, params) +
                     real_quark_leptons(s, mi2, pt2, th, ph, 1, params)
#endif
#ifdef DIPOLE
                     - DipQB_leptons(s, mi2, pt2, th, ph, 0, params) -
                     DipQB_leptons(s, mi2, pt2, th, ph, 1, params)
#endif
                         );

          sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *

                 (
#ifdef REAL
                     real_quarkb_leptons(s, mi2, pt2, th, ph, 0, params) +
                     real_quarkb_leptons(s, mi2, pt2, th, ph, 1, params)
#endif
#ifdef DIPOLE
                     - DipQA_leptons(s, mi2, pt2, th, ph, 0, params) -
                     DipQA_leptons(s, mi2, pt2, th, ph, 1, params)
#endif
                         );
        }
        if (isnan(sig)) {
          printf("pt2 = %.15f \n", pt2);
          exit(1);
        }
      }
    }
  }

  // Flux, symmetry factor, spin and color average.
  if (is_gaugino_gluino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    return (0.5 * 1.0 / 2.0 * 1.0 / 2.0 * 1.0 / (2.0 * s) *
            (4.0 * M_PI * aS(params->murs, params->set)) *
            (4.0 * M_PI * aS(params->murs, params->set)) * sig * djac / g3s);
  } else { // Drell-Yan like
    return dij * 0.5 * 4 / (96 * 2 * s) *
           (4 * M_PI * aS(params->murs, params->set)) * sig * djac;
  }
}

double IQ23(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj; // conversion of type void* into IOS*
  double xa, xb;
  // Different three-particle phase space for on-shell subtraction
  // for the associated gaugino-gluino production.
  double s2_23, t1_23, s1_23, phi_23, djac23;
  djac23 = dPS3_ONSHELL23(xa, xb, s2_23, t1_23, s1_23, phi_23, x, params);

#if !defined(DIPOLE) || !defined(REAL)
  djac23 = 0.0;
#endif

  const double s = xa * xb * params->sh;

  // PDFs.
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double sigOS = 0.0;

  // Color average for 2->3 and LO.
  double color_average_2to3 = 1.0 / 3.0 * 1.0 / 8.0;

  // Sums over all possible initial states.
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      params->in1 = i0;
      params->in2 = i1;

      // Tests charge conservation.
      if (is_charge_conserved(params->in1, params->in2, params->out1,
                              params->out2)) {

#ifdef ONSUB
        // On-shell subtraction for the associated gaugino-gluino production
        if (is_gaugino_gluino(params->out1, params->out2)) {

          if (abs(djac23) > 1E-15) {
            // Real quark emission where s23 = mSQs leads to resonance.
            sigOS += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
                     (color_average_2to3 *
                      (real_quark_gaugino_gluino_onshell_23(
                          s, s2_23, t1_23, s1_23, phi_23, 1, params)));

            // Real antiquark emission where s23 = mSQs leads to resonance.
            sigOS += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                     (color_average_2to3 *
                      (real_quarkb_gaugino_gluino_onshell_23(
                          s, s2_23, t1_23, s1_23, phi_23, 1, params)));
          }
        }
#endif
      }
    }
  }

  // Flux, symmetry factor, spin and color average.
  double g3s = std::norm(params->gqq[0][0].R);
  return (1.0 / 2.0 * 1.0 / 2.0 * 1.0 / (2.0 * s) *
          (4.0 * M_PI * aS(params->murs, params->set)) *
          (4.0 * M_PI * aS(params->murs, params->set)) * sigOS / g3s * djac23);
}

double IQ13(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj; // conversion of type void* into IOS*
  double xa, xb;
  // Different three-particle phase space for on-shell subtraction
  // for the associated gaugino-gluino production.
  double s2_13, t1_13, s1_13, phi_13, djac13;
  djac13 = dPS3_ONSHELL13(xa, xb, s2_13, t1_13, s1_13, phi_13, x, params);

#if !defined(DIPOLE) || !defined(REAL)
  djac13 = 0.0;
#endif

  const double s = xa * xb * params->sh;

  // PDFs.
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double sigOS = 0.0;

  // Color average for 2->3 and LO.
  double color_average_2to3 = 1.0 / 3.0 * 1.0 / 8.0;

  // Sums over all possible initial states.
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      params->in1 = i0;
      params->in2 = i1;

      // Tests charge conservation.
      if (is_charge_conserved(params->in1, params->in2, params->out1,
                              params->out2)) {

#ifdef ONSUB
        // On-shell subtraction for the associated gaugino-gluino production
        if (abs(djac13) > 1E-15) {
          // Real quark emission where s13 = mSQs leads to resonance.
          sigOS += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
                   (color_average_2to3 *
                    (real_quark_gaugino_gluino_onshell_13(
                        s, s2_13, t1_13, s1_13, phi_13, 2, params)));
          // Real antiquark emission where s13 = mSQs leads to resonance.
          sigOS += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                   (color_average_2to3 *
                    (real_quarkb_gaugino_gluino_onshell_13(
                        s, s2_13, t1_13, s1_13, phi_13, 2, params)));
        }
#endif
      }
    }
  }

  // Flux, symmetry factor, spin and color average.

  double g3s = std::norm(params->gqq[0][0].R);
  return (1.0 / 2.0 * 1.0 / 2.0 * 1.0 / (2.0 * s) *
          (4.0 * M_PI * aS(params->murs, params->set)) *
          (4.0 * M_PI * aS(params->murs, params->set)) * sigOS / g3s * djac13);
}

// Threshold resummation
double IR(double *x, size_t dim, void *prm) {
  Parameters *params = (Parameters *)prm;

  const double sh = params->sh;

  double m1;
  double m2;

  SET_MASS;

  const double m1s = pow2(m1);
  const double m2s = pow2(m2);

  // Conversion factor from GeV^-2 to pb.
  complex<double> djac(389379304.0, 0.0);

  // M2 integration.
  double M2min = pow2(m1 + m2);
  double M2max = sh;
  
  if (is_lepton_lepton(params->out1, params->out2)) {
    M2min = pow2(params->Minv_min);
    M2max = pow2(params->Minv_max);
  }
  
  double s = M2min * pow(M2max / M2min, x[2]);
  djac *= log(M2max / M2min);

  // Inverse Mellin.
  const complex<double> ephi(-M_SQRT1_2, M_SQRT1_2);
  const complex<double> nm =
      (M2min / sh / 0.09 + 0.09) - params->a1min - ephi * log(x[0]);
  djac *= ephi / x[0];

  // t integration.
  const double tmin = -0.5 * (s - m1s - m2s + kln(s, m1s, m2s));
  const double tmax = tmin + kln(s, m1s, m2s);
  double t = (tmax - tmin) * x[1] + tmin;

  djac *= tmax - tmin;
  djac /= 8.0 * M_PI * s;

  if (is_gaugino_gluino(params->out1, params->out2)) {
    // Ordinary NLL.
    return M_1_PI *
           imag(djac * pow(s / sh, -nm + 1.0) * Thadronic_nll_xs(nm, s, t, params));
  }
  // Collinear improved NLL.
  return M_1_PI *
         imag(djac * pow(s / sh, -nm + 1.0) * Thadronic_xs2(nm, s, t, params));
}

// Unimproved threshold resummation for Drell-Yan like processes.
double IR_nll_unimproved(double *x, size_t dim, void *prm) {
  Parameters *params = (Parameters *)prm;

  const double sh = params->sh;

  double m1;
  double m2;

  SET_MASS;

  const double m1s = pow2(m1);
  const double m2s = pow2(m2);

  // Conversion factor from GeV^-2 to pb.
  complex<double> djac(389379304.0, 0.0);

  // M2 integration.
  double M2min = pow2(m1 + m2);
  double M2max = sh;

  if (is_lepton_lepton(params->out1, params->out2)) {
    M2min = pow2(params->Minv_min);
    M2max = pow2(params->Minv_max);
  }

  double s = M2min * pow(M2max / M2min, x[2]);
  djac *= log(M2max / M2min);

  // Inverse Mellin.
  const complex<double> ephi(-M_SQRT1_2, M_SQRT1_2);
  const complex<double> nm =
      (M2min / sh / 0.09 + 0.09) - params->a1min - ephi * log(x[0]);
  djac *= ephi / x[0];

  // t integration.
  const double tmin = -0.5 * (s - m1s - m2s + kln(s, m1s, m2s));
  const double tmax = tmin + kln(s, m1s, m2s);
  double t = (tmax - tmin) * x[1] + tmin;
  djac *= tmax - tmin;
  djac /= 8.0 * M_PI * s;

  return M_1_PI *
         imag(djac * pow(s / sh, -nm + 1.0) * Thadronic_nll_xs(nm, s, t, params));
}

// Unimproved threshold resummation for Drell-Yan like processes.
double IR_nnll_unimproved(double *x, size_t dim, void *prm) {
  Parameters *params = (Parameters *)prm;

  const double sh = params->sh;

  double m1;
  double m2;

  SET_MASS;

  const double m1s = pow2(m1);
  const double m2s = pow2(m2);

  // Conversion factor from GeV^-2 to pb.
  complex<double> djac(389379304.0, 0.0);

  // M2 integration.
  double M2min = pow2(m1 + m2);
  double M2max = sh;

  if (is_lepton_lepton(params->out1, params->out2)) {
    M2min = pow2(params->Minv_min);
    M2max = pow2(params->Minv_max);
  }

  double s = M2min * pow(M2max / M2min, x[2]);
  djac *= log(M2max / M2min);

  // Inverse Mellin.
  const complex<double> ephi(-M_SQRT1_2, M_SQRT1_2);
  const complex<double> nm =
      (M2min / sh / 0.09 + 0.09) - params->a1min - ephi * log(x[0]);
  djac *= ephi / x[0];

  // t integration.
  const double tmin = -0.5 * (s - m1s - m2s + kln(s, m1s, m2s));
  const double tmax = tmin + kln(s, m1s, m2s);
  double t = (tmax - tmin) * x[1] + tmin;
  djac *= tmax - tmin;
  djac /= 8.0 * M_PI * s;

  return M_1_PI *
         imag(djac * pow(s / sh, -nm + 1.0) * Thadronic_nnll_xs(nm, s, t, params));
}

// Integration.
// If a specific integration is not stable increase the number of calls
// (3rd argument of Integration()).
void hadronic_xs(double &res, double &err, int Flag, Parameters *params) {

  switch (Flag) {
  case 0: // LO cross section.
    cout << endl;
    cout << "********************" << endl;
    cout << "* LO cross section *" << endl;
    cout << "********************" << endl;
    Integration(&IB, 3, 20000, params->precision, 1e-12, res, err, params, 0.1,
                5, 5);
    break;
  case 1: // virtual corrections + integrated dipole
    cout << endl;
    cout << "*******************************************" << endl;
    cout << "* virtual corrections + integrated dipole *" << endl;
    cout << "*******************************************" << endl;
    Integration(&IV, 3, 1200, params->precision, 1e-12, res, err, params, 0.2);
    break;
  case 2: // collinear remainder.
    cout << endl;
    cout << "*****************************" << endl;
    cout << "* collinear remainder (P+K) *" << endl;
    cout << "*****************************" << endl;
    Integration(&IC, 4, 25000, params->precision, 1e-12, res, err, params);
    break;
  case 3: // real gluon emission - dipole
    cout << endl;
    cout << "********************************" << endl;
    cout << "* real gluon emission - dipole *" << endl;
    cout << "********************************" << endl;
    Integration(&IG, 6, 25000, params->precision, 1e-12, res, err, params);
    break;
  case 4: // real quark emission - dipole - on-shell counterterm.
    cout << endl;
    cout << "**************************************************" << endl;
    cout << "* real light quark emission - dipole - (onshell) *" << endl;
    cout << "**************************************************" << endl;
    if (is_gaugino_gluino(params->out1, params->out2)) {
      double res13, res23, err13, err23;
      cout << endl;
      Integration(&IQ, 6, 50000, params->precision, 1e-12, res, err, params,
                  0.2, 5, 5);
      Integration(&IQ23, 6, 100000, params->precision, 1e-12, res23, err23,
                  params, 0.4, 5, 5);
      int calls = 50000;
      ;
      for (int i = 0; i < 12; i++) {
        if (i == 8 || i == 11) {
          continue; // no increase of calls needed
        }
        if (params->mGL < params->mSQ[i]) {
          calls = 150000; // on-shell remainder sizable
          break;
        }
      }
      Integration(&IQ13, 6, calls, params->precision, 1e-12, res13, err13,
                  params, 0.2, 5, 5);
      // add real quark emission contributions
      res = res + res13 + res23;
      err = sqrt(pow(err, 2) + pow(err13, 2) + pow(err23, 2));
    } else {
      Integration(&IQ, 6, 50000, params->precision, 1e-12, res, err, params,
                  0.1, 5, 5);
    }
    break;
  case 5: // threshold resummation (improved for DY; ordinary for
          // gaugino-gluino)
    cout << endl;
    cout << "************************************************" << endl;
    cout << "* threshold resummation (improved) - expansion *" << endl;
    cout << "************************************************" << endl;
    Integration(&IR, 3, 35000, params->precision, 1e-12, res, err, params, 0.2,
                5, 5);
    break;
  case 6: // ordinary threshold resummation (including nll logs)
    cout << endl;
    cout << "******************************************************" << endl;
    cout << "* threshold resummation (unimproved nll) - expansion *" << endl;
    cout << "******************************************************" << endl;
    Integration(&IR_nll_unimproved, 3, 35000, params->precision, 1e-12, res, err,
                params, 0.2, 5, 5);
    break;
  case 7: // ordinary threshold resummation (including nnll logs)
    cout << endl;
    cout << "*******************************************************" << endl;
    cout << "* threshold resummation (unimproved nnll) - expansion *" << endl;
    cout << "*******************************************************" << endl;
    Integration(&IR_nnll_unimproved, 3, 35000, params->precision, 1e-12, res, err,
                params, 0.2, 5, 5);
    break;
  default:
    cout << "hadronic_xs: Flag=" << Flag << endl;
    exit(0);
  }
  return;
}
