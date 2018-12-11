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

#define RESUMMINO_VERSION "1.0.0"

#include <complex>

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

#endif
