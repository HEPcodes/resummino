// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Mathematical functions.

#ifndef MATH_H_
#define MATH_H_

#include <complex>

using namespace std;

complex<double> Gamma(const complex<double> x);
complex<double> Psi(complex<double> x);
complex<double> Beta(const complex<double>, const complex<double>);

#endif
