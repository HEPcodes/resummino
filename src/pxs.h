// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

// Partonic cross sections at LO (Born), NLO and NLL.
//
// The functions with names ending in eps{1,2} contain the terms proportional to
// \epsilon^{1,2} in the D-dimensional matrix elements, where D=4-2\epsilon.

#ifndef PXS_H_
#define PXS_H_

#include "params.h"
#include <complex>

using namespace std;

// LO cross section.
double born_gauginos(const double, const double, Parameters *);
double born_gauginos_eps1(const double, const double, Parameters *);
double born_gauginos_eps2(const double, const double, Parameters *);

double born_sleptons(const double, const double, Parameters *);

double born_leptons(const double, const double, Parameters *);
double born_leptons_eps1(const double, const double, Parameters *);
double born_leptons_eps2(const double, const double, Parameters *);

double born_gagl(const double, const double, Parameters *); // Gaugino-gluino.

double born_gasq(const double, const double, Parameters *); // Gaugino-squark.

// Virtual corrections.
double Virt_gauginos(const double, const double, Parameters *);
double Virt2_gauginos(const double, const double, Parameters *);

double Virt_sleptons(const double, const double, Parameters *);
double Virt2_sleptons(const double, const double, Parameters *);

double Virt_leptons(const double, const double, Parameters *);
double Virt2_leptons(const double, const double, Parameters *);

double Virt_gaugino_gluino(const double, const double, Parameters *);

double Virt_gaugino_squark(const double, const double, Parameters *);

// Dipole substraction for virtual corrections.
double DipI_gauginos(const double, const double, Parameters *);
double DipI_sleptons(const double, const double, Parameters *);
double DipI_leptons(const double, const double, Parameters *);
double DipI_gagl(const double, const double, Parameters *);
double DipI_gasq(const double, const double, Parameters *);

// Real gaugino
double real_gluon_gauginos(const double, const double, const double,
                           const double, const double, const int, Parameters *);
double real_quark_gauginos(const double, const double, const double,
                           const double, const double, const int, Parameters *);
double real_quarkb_gauginos(const double, const double, const double,
                            const double, const double, const int,
                            Parameters *);

double real_quark_gaugino_gluino_onshell(const double, const double,
                                         const double, const double,
                                         const double, const int, Parameters *);
double real_quarkb_gaugino_gluino_onshell(const double, const double,
                                          const double, const double,
                                          const double, const int,
                                          Parameters *);

double real_quark_gaugino_gluino_onshell_13(const double, const double,
                                            const double, const double,
                                            const double, int, Parameters *);
double real_quarkb_gaugino_gluino_onshell_13(const double, const double,
                                             const double, const double,
                                             const double, int, Parameters *);

double real_quark_gaugino_gluino_onshell_23(const double, const double,
                                            const double, const double,
                                            const double, int, Parameters *);
double real_quarkb_gaugino_gluino_onshell_23(const double, const double,
                                             const double, const double,
                                             const double, int, Parameters *);

// Real gaugino gluino
double real_gluon_gaugino_gluino(const double, const double, const double,
                                 const double, const double, const int,
                                 Parameters *);
double real_quark_gaugino_gluino(const double, const double, const double,
                                 const double, const double, const int,
                                 Parameters *);
double real_quarkb_gaugino_gluino(const double, const double, const double,
                                  const double, const double, const int,
                                  Parameters *);

double real_quark_gaugino_gluino_os(const double, const double, const double,
                                    const double, const double, const int,
                                    Parameters *);
double real_quarkb_gaugino_gluino_os(const double, const double, const double,
                                     const double, const double, const int,
                                     Parameters *);

double DipGAB_GLGA(const double, const double, const double, const double,
                   const double, const int, Parameters *);
double DipGBA_GLGA(const double, const double, const double, const double,
                   const double, const int, Parameters *);
double DipGA1_GLGA(const double, const double, const double, const double,
                   const double, const int, Parameters *);
double DipGB1_GLGA(const double, const double, const double, const double,
                   const double, const int, Parameters *);
double DipG1A_GLGA(const double, const double, const double, const double,
                   const double, const int, Parameters *);
double DipG1B_GLGA(const double, const double, const double, const double,
                   const double, const int, Parameters *);

// Real squark
double real_gluon_gaugino_squark(const double, const double, const double,
                                 const double, const double, const int,
                                 Parameters *);

// Real leptons
double real_gluon_leptons(const double, const double, const double,
                          const double, const double, const int, Parameters *);
double real_quark_leptons(const double, const double, const double,
                          const double, const double, const int, Parameters *);
double real_quarkb_leptons(const double, const double, const double,
                           const double, const double, const int, Parameters *);
double DipGA_leptons(const double, const double, const double, const double,
                     const double, const int, Parameters *);
double DipGB_leptons(const double, const double, const double, const double,
                     const double, const int, Parameters *);
double DipQA_leptons(const double, const double, const double, const double,
                     const double, const int, Parameters *);
double DipQB_leptons(const double, const double, const double, const double,
                     const double, const int, Parameters *);

// Real sleptons
double real_gluon_sleptons(const double, const double, const double,
                           const double, const double, const int, Parameters *);
double real_quark_sleptons(const double, const double, const double,
                           const double, const double, const int, Parameters *);
double real_quarkb_sleptons(const double, const double, const double,
                            const double, const double, const int,
                            Parameters *);
double DipGA_sleptons(const double, const double, const double, const double,
                      const double, const int, Parameters *);
double DipGB_sleptons(const double, const double, const double, const double,
                      const double, const int, Parameters *);
double DipQA_sleptons(const double, const double, const double, const double,
                      const double, const int, Parameters *);
double DipQB_sleptons(const double, const double, const double, const double,
                      const double, const int, Parameters *);

// Resummation
complex<double> Thadronic_nll_xs(const complex<double>, const double, const double,
                             Parameters *);
complex<double> Thadronic_nnll_xs(const complex<double>, const double, const double,
                             Parameters *);
complex<double> Thadronic_xs2(const complex<double>, const double, const double,
                              Parameters *);
complex<double> JtXS2(const complex<double>, const complex<double>,
                      const double, const double, Parameters *);
complex<double> PtXS(const complex<double>, const complex<double>, const double,
                     const double, Parameters *);
complex<double> JtXS(const complex<double>, const complex<double>, const double,
                     const double, Parameters *);

// sja
double sja(const double S, const double T, Parameters *params);
double sjb(const double S, const double T, Parameters *params);
double sab(const double S, const double T, Parameters *params);

#endif
