// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Computes the virtual part of the partonic cross section.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "kinematics.h"
#include "params.h"
#include "pxs.h"
#include "utils.h"

// IEPS 0 (finite part), 1 (1/eps coefficient) and 2 (1/eps^2 coefficient).
#define IEPS 0

// Finite shift from DREG to DRED.
#define DRBAR

// Decouple squarks, gluinos and top from the aS running.
// (see arxiv 9610490v1 p. 17 above eq. 19)
#define DECOUPLE_HEAVY

#define QSELF // external quark self-enery

#define SS

// left- and right-handed "blind" virtual corrections
double qcd_gasq(const double S, const double T, Parameters *params) {

  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = propagator_charge(aa, bb, ii, jj, CHANNEL_S);
  int qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  FI *ff = new FI();

  BORN_KINEMATICS;

  if (ii > 30) {
    ii = ii - 31;
  } else {
    jj = jj - 31;
  }

  struct Coupling Cb[4] = {0, 0};

#ifdef SS
  ff->SetPropagator(0, 0, 0, 0);
  // p2 is gaugino -> Squark production
  if (params->out2 <= 10) {
    Cb[0] = params->CHSQq[jj][ii][aa];
    Cb[1] = params->CHSQq[jj][ii][aa];
    ff->SetBCoupling(Cb);
    virt +=
        2.0 * ff->Mss_qqg1_SQGA(0.0, 0.0, 2.0 * ff->papb, 0.0, 0.0, 0.0, IEPS);

    // p1 is gaugino -> Antisquark production
  } else if (params->out1 <= 10) {
    Cb[0] = params->CHqSQ[ii][bb][jj];
    Cb[1] = params->CHqSQ[ii][bb][jj];
    ff->SetBCoupling(Cb);
    virt +=
        2.0 * ff->Mss_qqg1_SQGA(0.0, 0.0, 2.0 * ff->papb, 0.0, 0.0, 0.0, IEPS);
    // virt += ff->Mss_SQGA2();
  }

#endif

  delete ff;
  return virt;
}

double Virt_gaugino_squark(const double S, const double T, Parameters *params) {

  double result = 0.0;

  result += qcd_gasq(S, T, params);

  return result;
}

double DipI_gasq(const double S, const double T, Parameters *params) {
  return 0.0;
}
