// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Command-line interface.

// File incldudes catani-seymour dipoles
// based on hep-ph/0201036 and hep-ph/9605323.

// TODO: LASSE (add dipoles and restructure)

#include "dipoles.h"
#include "gsl_all.h"
#include "utils.h"
#include <iostream>

// color factors
#define cA 3.0
#define cF (4.0 / 3.0)
#define I2R (0.5) // SU(3) normalization factor
#define nf 5.0    // number of light flavors

// Flavor constants given in hep-ph/0201036 p.53 eq. C.12 and p.32 eq. 5.91 and
// p.37 eq. 6.17
#define GAMMA_Q (3.0 / 2.0 * cF)
#define GAMMA_GL (3.0 / 2.0 * cA)
#define GAMMA_SQ (2.0 * cF)
#define GAMMA_G (11.0 / 6.0 * cA - 2.0 / 3.0 * I2R * nf)
#define K_Q (7.0 / 2.0 - pow2(M_PI) / 6.0) * cF
#define K_GL (7.0 / 2.0 - pow2(M_PI) / 6.0) * cA
#define K_SQ (4.0 - pow2(M_PI) / 6.0) * cF
#define K_G ((67.0 / 18.0 - pow2(M_PI) / 6.0) * cA - 10.0 / 9.0 * I2R * nf)

// Regulated Altarelli Parisi functions (hep-ph/0201036 p. 31, eq. 5.89).
double P_aap_reg(PartonType a, PartonType ap, double x) {
  if ((a == PARTON_QUARK && ap == PARTON_QUARK) ||
      (a == PARTON_ANTIQUARK && ap == PARTON_ANTIQUARK)) {
    return -cF * (1.0 + x); // P_qq_reg
  }
  if ((a == PARTON_QUARK && ap == PARTON_GLUON) ||
      (a == PARTON_ANTIQUARK && ap == PARTON_GLUON)) {
    return cF * (1.0 + pow2(1.0 - x)) / x; // P_qg_reg = P_qapg_reg
  }
  if ((a == PARTON_GLUON && ap == PARTON_QUARK) ||
      (a == PARTON_GLUON && ap == PARTON_ANTIQUARK)) {
    return I2R * (pow2(x) + pow2(1.0 - x)); // P_gq_reg = P_qqap_reg
  }
  if (a == PARTON_GLUON && ap == PARTON_GLUON) {
    return 2.0 * cA * ((1.0 - x) / x - 1.0 + x * (1.0 - x)); // P_gg_reg
  }
}

// General expression for the Altarelli Parisi splitting function
// (hep-ph/0201036 p. 32, eq. 5.94).
double P_aap(PartonType a, PartonType ap, FunctionType function_type, double x,
             double z) {

  // for splitting to same flavor a == ap; (Kronecker Delta)
  int delta_aap;
  if (a == ap) {
    delta_aap = 1;
  } else
    delta_aap = 0;

  // Ta * Ta for a == quark or a == gluon (Expectation values of the color
  // operators)
  double color;
  if (a == PARTON_QUARK || a == PARTON_ANTIQUARK) {
    color = cF;
  } else
    color = cA;

  // gamma_a for quark or gluon
  double gamma_a = 0.0;
  if (a == PARTON_QUARK || a == PARTON_ANTIQUARK) {
    gamma_a = GAMMA_Q;
  } else if (a == PARTON_GLUON) {
    gamma_a = GAMMA_G;
  } else if (a == PARTON_GLUINO) {
    gamma_a = GAMMA_GL;
  }

  switch (function_type) {
  case FUNCTION_ALL_WDELTA:
    // all parts which are not proportional to delta(1-x)
    return P_aap_reg(a, ap, x) + delta_aap * 2.0 * color * 1.0 / (1.0 - x);
  case FUNCTION_DELTA:
    // part proportional to delta
    // (because we integrate from z to 1 we have to subtract the integral
    // from 0 to z over 2*color*1/(1-x) leading to -2.0 * color * log(1-z) )
    return delta_aap * (gamma_a + 2.0 * color * log(1 - z));
  case FUNCTION_PLUS:
    // the plus distribution part (if there are prefactors
    // which include x and do not belong to the plus-dist. we set them to 1)
    return delta_aap * (2.0 * color * 1.0 / (1.0 - x));
  }
}

// Epsilon dependent parts of d-dimensional AP splitting functions
// (hep-ph/0201036 p. 32 eq 5.93).
double Php_aap(PartonType a, PartonType ap, double x) {
  if ((a == PARTON_QUARK && ap == PARTON_QUARK) ||
      (a == PARTON_ANTIQUARK && ap == PARTON_ANTIQUARK)) {
    return cF * (1 - x);
  } else if ((a == PARTON_QUARK || a == PARTON_ANTIQUARK) &&
             ap == PARTON_GLUON) {
    return cF * x;
  } else if (a == PARTON_GLUON &&
             (ap == PARTON_QUARK || ap == PARTON_ANTIQUARK)) {
    return 2.0 * I2R * x * (1.0 - x);
  } else if (a == PARTON_GLUON && ap == PARTON_GLUON) {
    return 0;
  }
}

// General expression for K_ab_bar (hep-ph/0201036 p. 44, eq. 6.56).
double K_aap_bar(PartonType a, PartonType ap, FunctionType function_type,
                 double x, double z) {

  // for splitting to same flavor a == ap;
  int delta_aap;
  if (a == ap) {
    delta_aap = 1.0;
  } else
    delta_aap = 0.0;

  // Ta * Ta for a == quark or a == gluon
  double color;
  if (a == PARTON_QUARK || a == PARTON_ANTIQUARK) {
    color = cF;
  } else
    color = cA;

  // gamma_a for quark or gluon
  double gamma_a = 0.0;
  if (a == PARTON_QUARK || a == PARTON_ANTIQUARK) {
    gamma_a = GAMMA_Q;
  } else if (a == PARTON_GLUON) {
    gamma_a = GAMMA_G;
  } else if (a == PARTON_GLUINO) {
    gamma_a = GAMMA_GL;
  }

  // K_a for quark or gluon
  double K_a = 0.0;
  if (a == PARTON_QUARK || a == PARTON_ANTIQUARK) {
    K_a = K_Q;
  } else if (a == PARTON_GLUON) {
    K_a = K_G;
  } else if (a == PARTON_GLUINO) {
    K_a = K_GL;
  }

  switch (function_type) {
  case FUNCTION_ALL_WDELTA:
    return P_aap_reg(a, ap, x) * log((1.0 - x) / x) + Php_aap(a, ap, x) +
           delta_aap * color * 2.0 / (1.0 - x) * log((1 - x) / x);
  case FUNCTION_DELTA:
    //-integral cF * 2.0/(1.0 - x) * log((1-x)/x) from 0 to z
    // (solved with mathematica and a substitution 1-x -> y to "help"
    // mathematica)
    return delta_aap *
           (-color / 3.0 * (pow(M_PI, 2) - 3. * pow(log(1 - z), 2) -
                            6. * gsl_sf_dilog(1 - z)) -
            (gamma_a + K_a -
             5.0 / 6.0 * pow2(M_PI) * color)); // explicit delta part
  case FUNCTION_PLUS:
    return delta_aap * color * 2.0 / (1.0 - x) * log((1 - x) / x);
  }
}

// J_a_ij (hep-ph/0201036 p. 27, eq. 5.58)
double J_a_ij(PartonType i, PartonType j, double muj,
              FunctionType function_type, double x, double z) {
  // J_gQ (J_gQb)
  if (i == PARTON_GLUON &&
      (j == PARTON_QUARK || j == PARTON_ANTIQUARK || j == PARTON_GLUINO)) {
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return (1.0 - x) / (2.0 * pow2(1.0 - x + pow2(muj))) -
             2.0 / (1.0 - x) * (1.0 + log(1.0 - x + pow2(muj))) +
             2.0 / (1.0 - x) * log(2.0 + pow2(muj) - x);
    case FUNCTION_DELTA: // again substitution used in mathematica 1-x -> y to
                         // obtain a "nice" result
      return -(-((pow(muj, 2) * z) /
                 ((1 + pow(muj, 2)) * (1 + pow(muj, 2) - z))) +
               4. * (1 + 2 * log(muj)) * log(1 - z) -
               4. * log(1 + pow(muj, 2)) * log(1 - z) -
               log(1 - z / (1. + pow(muj, 2))) +
               4. * gsl_sf_dilog(-pow(muj, -2)) -
               4. * gsl_sf_dilog((-1 + z) / pow(muj, 2))) /
             2.;
    case FUNCTION_PLUS:
      return (1.0 - x) / (2.0 * pow2((1.0 - x + pow2(muj)))) -
             2.0 / (1.0 - x) * (1.0 + log(1.0 - x + pow2(muj))) +
             2.0 / (1.0 - x) * log(2.0 + pow2(muj) - 1.0);
      // the last log is not part of the plus distribution, therefore there is a
      // 1 instead of an x
    }
  }
}

// (hep-ph/0201036 p. 45, eq. 6.57 - 6.60)
// No general expression.
// later u have to add Kappa_ab_sq in C.15 p.54 (TODO: LASSE)
double Kappa_aap_j(PartonType a, PartonType ap, PartonType j, double mj,
                   double sja, FunctionType function_type, double x, double z) {

  if (a == PARTON_GLUON && ap == PARTON_QUARK &&
      (j == PARTON_QUARK || j == PARTON_GLUINO)) { // Kappa_gq_q or Kappa_gq_gl
    return 0;
  } else if (a == PARTON_QUARK && ap == PARTON_QUARK &&
             (j == PARTON_QUARK ||
              j == PARTON_GLUINO)) { // Kappa_qq_q or Kappa_qq_gl

    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return 2.0 * (log(1 - x) / (1 - x) - log(2 - x) / (1 - x)) +
             J_a_ij(PARTON_GLUON, j, mj / sqrt(sja), FUNCTION_ALL_WDELTA, x,
                    z) +
             2.0 * 1.0 / (1.0 - x) *
                 log((2.0 - x) * sja / ((2.0 - x) * sja + pow2(mj)));
    case FUNCTION_DELTA:
      return (log(1 - z) * (2 * log(sja / (pow(mj, 2) + sja)) + log(1 - z))) -
             GAMMA_Q / cF + pow2(mj) / sja * log(pow2(mj) / (sja + pow2(mj))) +
             0.5 * pow2(mj) / (sja + pow2(mj)) +
             J_a_ij(PARTON_GLUON, j, mj / sqrt(sja), FUNCTION_DELTA, x, z);
    case FUNCTION_PLUS:
      return 2.0 * (log(1 - x) / (1 - x)) +
             J_a_ij(PARTON_GLUON, j, mj / sqrt(sja), FUNCTION_PLUS, x, z) +
             2.0 * 1 / (1 - x) * log(sja / (sja + pow2(mj)));
    }
  } else if (a == PARTON_QUARK && ap == PARTON_GLUON &&
             (j == PARTON_QUARK ||
              j == PARTON_GLUINO)) { // Kappa_qg_q or Kappa_qg_gl
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return 2.0 * cF / cA * pow2(mj) / (x * sja) *
             log(pow2(mj) / ((1 - x) * sja + pow2(mj)));
    case FUNCTION_DELTA:
      return 0.0;
    case FUNCTION_PLUS:
      return 0.0;
    }
  } else if (a == PARTON_GLUON && ap == PARTON_GLUON &&
             (j == PARTON_QUARK ||
              j == PARTON_GLUINO)) { // Kappa_gg_q or Kappa_gg_gl
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return 2.0 * pow2(mj) / (x * sja) *
                 log(pow2(mj) / ((1 - x) * sja) + pow2(mj)) // Kappa_qg_q
             + 2.0 * (log(1 - x) / (1 - x) - log(2 - x) / (1 - x)) +
             J_a_ij(PARTON_GLUON, PARTON_QUARK, mj / sqrt(sja),
                    FUNCTION_ALL_WDELTA, x, z) +
             2.0 * 1 / (1 - x) *
                 log((2 - x) * sja / ((2 - x) * sja + pow2(mj))); // Kappa_gg_q
    case FUNCTION_DELTA:
      return (log(1 - z) * (2 * log(sja / (pow(mj, 2) + sja)) + log(1 - z))) -
             GAMMA_Q / cF + pow2(mj) / sja * log(pow2(mj) / (sja + pow2(mj))) +
             0.5 * pow2(mj) / (sja + pow2(mj)) +
             J_a_ij(PARTON_GLUON, PARTON_QUARK, mj / sqrt(sja), FUNCTION_DELTA,
                    x, z);
    case FUNCTION_PLUS:
      return 2.0 * (log(1 - x) / (1 - x)) +
             J_a_ij(PARTON_GLUON, PARTON_QUARK, mj / sqrt(sja), FUNCTION_PLUS,
                    x, z) +
             2.0 * 1 / (1 - x) * log(sja / (sja + pow2(mj)));
    }
  }
}

// For the color operators in the formula see hep-ph/9605323 p.97 eq A.4.
// P_aaprime in (hep-ph/0201036 p. 46, eq. 6.67 and 6.53).
double P_bold_aap_b_j(PartonType a, PartonType ap, PartonType b, PartonType j,
                      FunctionType function_type, double x, double z,
                      double sja, double sab, double mufs, double alphas,
                      EmissionType emission) {

  // some cuttoff to avoid numerical issues
  if ((1.0 - x) < 1.0E-10) {
    return 0.0;
  }

  double color1 = 0.0; // will be Tj * Tap/Tap^2
  double color2 = 0.0; // will be Tb * Tap/Tap^2

  if (ap == PARTON_QUARK && b == PARTON_ANTIQUARK && j == PARTON_QUARK) {
    color1 = -1.0 / 2.0;
    color2 = -1.0 / 2.0;
  } else if (ap == PARTON_QUARK && b == PARTON_ANTIQUARK &&
             j == PARTON_GLUINO) {
    color1 = -1.0 / cF * cA / 2.0;
    color2 = +1.0 / cF * (cA - 2 * cF) / 2.0;
  } else if (ap == PARTON_QUARK && b == PARTON_ANTIQUARK && j == PARTON_NONE) {
    color1 = 0.0;
    color2 = -1.0;
  }

  if (emission == INITIAL_AND_FINAL) {
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * color1 *
                 log(mufs / (x * sja)) +
             alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * color2 *
                 log(mufs / (x * sab));
    case FUNCTION_DELTA:
      return alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * color1 *
                 log(mufs / sja) +
             alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * color2 *
                 log(mufs / sab);
    case FUNCTION_PLUS:
      return alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * color1 *
                 log(mufs / sja) +
             alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * color2 *
                 log(mufs / sab);
    }

  } else if (emission == INITIAL) {
    // color2 = -1;
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return color2 * alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) *
             log(mufs / (x * sab));
    case FUNCTION_DELTA:
      return color2 * alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) *
             log(mufs / sab);
    case FUNCTION_PLUS:
      return color2 * alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) *
             log(mufs / sab);
    }
  }
}

// p.44 and 46
// eq 6.55 + 6.68
double K_bold_aap_b_j(PartonType a, PartonType ap, PartonType b, PartonType j,
                      FunctionType function_type, double x, double z, double mj,
                      double sja, double sab, double alphas,
                      EmissionType emission) {

  // cuttoff to avoid numerical instabilities
  if ((1.0 - x) < 1.0E-10) {
    return 0.0;
  }

  double gamma_a = 0.0;
  if (a == PARTON_QUARK || a == PARTON_ANTIQUARK) {
    gamma_a = GAMMA_Q;
  } else if (a == PARTON_GLUON) {
    gamma_a = GAMMA_G;
  }

  int delta_aap;
  if (a == ap) {
    delta_aap = 1.0;
  } else
    delta_aap = 0.0;

  // Tj Tap/Tap^2, Tj Tap, Tb Tap/Tap^2, Tb Tap
  double color_factor1, color_factor2, color_factor3, color_factor4;
  if (ap == PARTON_QUARK && b == PARTON_ANTIQUARK && j == PARTON_QUARK) {
    color_factor1 = -1.0 / 2.0;
    color_factor2 = -1.0 / 2.0 * cF;
    color_factor3 = -1.0 / 2.0;
    color_factor4 = -1.0 / 2.0 * cF;
  } else if (ap == PARTON_QUARK && b == PARTON_ANTIQUARK &&
             j == PARTON_GLUINO) {
    color_factor1 = -cA / 2.0 * 1.0 / cF;
    color_factor2 = -cA / 2.0;
    color_factor3 = (cA - 2.0 * cF) / 2.0 * 1.0 / cF;
    color_factor4 = (cA - 2.0 * cF) / 2.0;
  } else if (ap == PARTON_QUARK && b == PARTON_ANTIQUARK && j == PARTON_NONE) {
    color_factor1 = 0.0;
    color_factor2 = 0.0;
    color_factor3 = -1.0;
    color_factor4 = -cF;
  }
  if (emission == INITIAL_AND_FINAL) {
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return alphas / (2 * M_PI) *
                 (K_aap_bar(a, ap, function_type, x, z) -
                  color_factor2 *
                      Kappa_aap_j(a, ap, j, mj, sja, function_type, x, z) -
                  color_factor1 *
                      (P_aap_reg(a, ap, x) *
                       log((1 - x) * sja / ((1 - x) * sja + pow2(mj))))) -
             alphas / (2 * M_PI) *
                 (color_factor3 * P_aap_reg(a, ap, x) * log(1 - x) +
                  color_factor4 * delta_aap * 2.0 * log(1 - x) / (1 - x));
    case FUNCTION_DELTA:
      return -alphas / (2 * M_PI) * color_factor1 * delta_aap * gamma_a *
                 (log((sja - 2 * mj * sqrt(sja + pow2(mj)) + 2 * pow2(mj)) /
                      sja) +
                  2 * mj / (sqrt(sja + pow2(mj)) + mj)) -
             alphas / (2 * M_PI) * color_factor4 * delta_aap *
                 (+pow(log(1 - z), 2) - pow2(M_PI) / 3) +
             alphas / (2 * M_PI) *
                 (K_aap_bar(a, ap, function_type, x, z) -
                  color_factor2 *
                      Kappa_aap_j(a, ap, j, mj, sja, function_type, x, z));
    case FUNCTION_PLUS:
      return alphas / (2 * M_PI) *
             (K_aap_bar(a, ap, function_type, x, z) -
              color_factor2 *
                  Kappa_aap_j(a, ap, j, mj, sja, function_type, x, z) -
              color_factor4 * delta_aap * 2.0 * log(1 - x) / (1 - x));
    }
  }
  if (emission == INITIAL) {
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return alphas / (2 * M_PI) * (K_aap_bar(a, ap, function_type, x, z)) -
             alphas / (2 * M_PI) *
                 (color_factor3 * P_aap_reg(a, ap, x) * log(1 - x) +
                  color_factor4 * delta_aap * 2.0 * log(1 - x) / (1 - x));
    case FUNCTION_DELTA:
      return alphas / (2 * M_PI) * (-color_factor4) * delta_aap *
                 (pow(log(1 - z), 2) - pow2(M_PI) / 3) +
             alphas / (2 * M_PI) * (K_aap_bar(a, ap, function_type, x, z));
    case FUNCTION_PLUS:
      return alphas / (2 * M_PI) *
             (K_aap_bar(a, ap, function_type, x, z) +
              (-color_factor4) * delta_aap * 2.0 * log(1 - x) / (1 - x));
    }
  }
}

// Splitting functions V for initial state emitter ai and specator b (V_ai_b)
// hep-ph/9605323 p.40 eq. 5.145ff
double V_qg_b(double x, double alpha) {
  return 8.0 * M_PI * alpha * (2.0 / (1.0 - x) - (1.0 + x)); // cF
}

// gluon emits quark or antiquark and the counterpart goes into the hard process
double V_gq_b(double x, double alpha) {
  return 8.0 * M_PI * alpha * (1.0 - 2.0 * x * (1.0 - x)); // TR
}

// this is for a gluon going into the hard process (not needed yet; check and
// fix)
double V_qq_b(double x, double alpha, double papi, double papb, double pbpi) {
  return 8.0 * M_PI * alpha *
         (-4.0 * x + (1 - x) / x * (2.0 * papb / (papi * pbpi))); // cF
}

double V_gg_b(double x, double alpha, double papi, double papb, double pbpi) {
  return 16.0 * alpha * (-4.0 * (x / (1.0 - x) + x * (1.0 - x)) +
                         (1.0 - x) / x * papb / (papi * pbpi) * (-2.0) * papi /
                             papb * pbpi); // cA
}

// Splitting functions V for initial state emitter ai and final state spectator
// j (V_ai_j)
// hep-ph/0201036 p. 30 eq 5.81ff
double V_qg_j(double x, double alpha, double zj) {
  return 8.0 * M_PI * alpha * (2.0 / (2.0 - x - zj) - 1.0 - x); // cF
}

double V_gqb_j(double x, double alpha) {
  return 8 * M_PI * alpha * (1.0 - 2.0 * x * (1.0 - x)); // I2R
}

// or better 5.84 ? (think so, needed for squark prodction) (TODO: LASSE)
// needed for associated squark gaugino production (fix)
double V_qq_j(double x, double alpha, double zi, double zj, double pjpi) {
  return 8 * M_PI * alpha * (-4.0 * x +
                             (1.0 - x) / x * 2.0 * zi * zj / pjpi *
                                 (-2.0 * pjpi / (zi * zj))); // cF
}

// Splitting functions V for final state emitter and initial state spectator
// (V_a_ij)
// hep-ph/0201036 p. 26 eq 5.50ff.
double V_a_gQ(double x, double alpha, double zj, double mQ, double pjpi) {
  // (cF for quarks but cA for gluinos)
  return 8. * M_PI * alpha * (2. / (2. - x - zj) - 1. - zj - pow2(mQ) / pjpi);
}

// eq.5.52 expanded in epsilon
double V_a_QQB(double x, double alpha, double zplus, double zminus, double zi) {
  return 8. * M_PI * alpha * (1. - 2. * (zplus - zi) * (zi - zminus)); // TR
}

// Expressions without the color operators
// p. 39 eq. 5.136
// initial emitter and specator
double D_ai_b(PartonType a, PartonType i, double papi, double papb, double pbpi,
              double x, double alpha, double color) {
  if ((a == PARTON_QUARK || a == PARTON_ANTIQUARK) && i == PARTON_GLUON) {
    return -1.0 / (2. * papi) * 1.0 / x * color * V_qg_b(x, alpha);
  } else if ((a == PARTON_GLUON && i == PARTON_QUARK) ||
             (a == PARTON_GLUON && i == PARTON_ANTIQUARK)) {
    return -1.0 / (2.0 * papi) * 1.0 / x * color * V_gq_b(x, alpha);
  } else if (a == PARTON_QUARK && i == PARTON_QUARK) {
    return -1.0 / (2.0 * papi) * 1.0 / x * color *
           V_qq_b(x, alpha, papi, papb, pbpi);
  }
}

// initial emitter and final spectator
// p.29 eq.5.71
double D_ai_j(PartonType a, PartonType i, double zi, double zj, double papi,
              double pjpi, double x, double alpha, double color) {
  if ((a == PARTON_QUARK || a == PARTON_ANTIQUARK) && i == PARTON_GLUON) {
    return -1.0 / (2.0 * papi) * 1.0 / x * color * V_qg_j(x, alpha, zj);
  } else if ((a == PARTON_GLUON && i == PARTON_ANTIQUARK) ||
             (a == PARTON_GLUON && i == PARTON_QUARK)) {
    return -1.0 / (2.0 * papi) * 1.0 / x * color * V_gqb_j(x, alpha);
  } else if (a == PARTON_QUARK && i == PARTON_QUARK) {
    return -1.0 / (2.0 * papi) * 1.0 / x * color *
           V_qq_j(x, alpha, zi, zj, pjpi);
  }
}

// Final emitter and initial spectator
// p.25 eq 5.40
double D_ij_a(PartonType i, PartonType j, double mi, double mj, double mij,
              double pjpi, double x, double alpha, double zj, double mQ,
              double zi, double zplus, double zminus, double color) {
  if (i == PARTON_GLUON && (j == PARTON_QUARK || j == PARTON_GLUINO)) {
    return -1.0 / (pow2(mi) + pow2(mj) + 2.0 * pjpi - pow2(mij)) * 1.0 / x *
           color * V_a_gQ(x, alpha, zj, mQ, pjpi);
  } else if (i == PARTON_GLUON && j == PARTON_QUARK) {
    // not needed now and maybe wrong (check!) (needed for squark gaugino
    // production)
    return -1.0 / (pow2(mi) + pow2(mj) + 2 * pjpi - pow2(mij)) * 1.0 / x *
           color * V_a_QQB(x, alpha, zplus, zminus, zi);
  }
}
