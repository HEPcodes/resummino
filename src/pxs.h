// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Computes the partonic cross section at LO, NLO and NLL level.

#ifndef PXS_H_
#define PXS_H_

#include <complex>
#include "prm.h"

using namespace std;

// LO cross section
double Born(const double, const double, Parameters*);
double Born1(const double, const double, Parameters*);
double Born2(const double, const double, Parameters*);

// Virtual
double Virt(const double, const double, Parameters*);
double Virt2(const double, const double, Parameters*);

// Dipole substraction
double DipI(const double, const double, Parameters*);

// Real
double RealG(const double, const double, const double,
             const double, const double, const int, Parameters*);
double RealQ(const double, const double, const double,
             const double, const double, const int, Parameters*);
double RealQB(const double, const double, const double,
              const double, const double, const int, Parameters*);
double DipGA(const double, const double, const double,
             const double, const double, const int, Parameters*);
double DipGB(const double, const double, const double,
             const double, const double, const int, Parameters*);
double DipQA(const double, const double, const double,
             const double, const double, const int, Parameters*);
double DipQB(const double, const double, const double,
             const double, const double, const int, Parameters*);

// Resummation
complex<double> Thadronic_xs(const complex<double>, const double, const double,
                             Parameters*);
complex<double> Thadronic_xs2(const complex<double>, const double, const double,
                              Parameters*);
complex<double> JtXS2(const complex<double>, const complex<double>,
                      const double, const double, Parameters *);
complex<double> PtXS(const complex<double>, const complex<double>,
                     const double, const double, Parameters *);
complex<double> JtXS(const complex<double>, const complex<double>,
                     const double, const double, Parameters *);

#endif
