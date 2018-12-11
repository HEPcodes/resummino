// This file is part of Resummino.
//
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Hadron PDF and the corresponding \alpha_S.

#ifndef HXS_H_
#define HXS_H_

void hadronic_xs(double&, double&, double&, int, int, Parameters*);
void hadronic_xs_dPT2(double&, double&, double&, int, int, Parameters*);
void hadronic_xs_dlnM2(double&, double&, double&, int, int, Parameters*);
double* get_masses(Parameters*);
#endif
