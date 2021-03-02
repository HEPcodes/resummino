// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Partonic cross section.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <stdio.h>

#include "dipoles.h"
#include "kinematics.h"
#include "params.h"
#include "utils.h"

// TODO: LASSE
double real_gluon_gaugino_squark(const double S, const double M2,
                                 const double PT2, const double TH,
                                 const double PH, const int YS,
                                 Parameters *params) {
  double born = 0.0;
  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  // Set different kinematics for 2 -> 3 process
  FI *ff = new FI();
  JET_KINEMATICS

  if (ii > 30) {
    ii = ii - 31;
  } else {
    jj = jj - 31;
  }

  struct Coupling Cw[4] = {0, 0, 0, 0};
  Cw[0] = params->CHSQq[jj][ii][aa];
  Cw[1] = params->gSQSQ[ii][ii];
  Cw[2] = params->CHSQq[jj][ii][aa];
  Cw[3] = params->gSQSQ[ii][ii];

  ff->SetPropagator(params->mSQ[ii], params->mSQ[ii], params->mSQ[ii] * 1.0e-2,
                    params->mSQ[ii] * 1.0e-2);

  ff->SetWCoupling(Cw);
  born += ff->MG_SQGA1();

  delete ff;
  return born;
}
