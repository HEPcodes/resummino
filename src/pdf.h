// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Returns the PDFs for gluons and quarks and the corresponding \alpha_S.
// Uses LHAPDF to retrieve the PDFs in the usual space, and for conjugate
// spaces, it performs a numerical fit.

#ifndef PDF_H_
#define PDF_H_

#include <complex>
#include <cmath>

using namespace std;

#include "LHAPDF/LHAPDF.h"

// Returns the value of \alpha_S/2\pi at energy Q.
static inline double aS(const double Qsq, const int set)
{
    return LHAPDF::alphasPDF(std::sqrt(Qsq)) / (2.0 * M_PI);
}


void pdfX(double&, double[2][6], double, double);

void pdfN(std::complex<double> &, std::complex<double>[2][6],
          const std::complex<double>, double[8][8]);
void pdfEvolve(complex<double> &gg, complex<double> qq[2][6],
               complex<double> nm, complex<double> lmbd);
void pdfFit(double&, double[8][8], double, double);

#endif
