// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Resummation module.

#include <cmath>
#include <complex>
#include <iostream>
#include "gsl_all.h"
#include "prm.h"
#include <stdio.h>
#include "utils.h"
#include "pxs.h"
#include "pdf.h"
#include "mth.h"

using namespace std;

#define NFLVR 5 // Number of active flavours [ 5 RECOMMENDED ]

// Threshold resummation.
complex<double> ThG(const complex<double> LL, const double S,
                    const double Muf2, const double Mur2, const int set)
{
    // Coeff. of the beta function
    const int nf = NFLVR;
    const double dnf = (double)nf;
    const double beta0 = 5.5 - dnf / 3.;
    const double beta1 = 25.5 - 9.5 * dnf / 3.;

    // Sudakov coeff. for quarks
    const double A1 = 8. / 3. / beta0;
    const double A2 = .5 * A1 * (67. / 3. - pow(M_PI, 2) - 10. / 9.*dnf) / beta0;

    // Logarithms
    const complex<double> lmbd = aS(Mur2, set) * beta0 * LL;
    const complex<double> lnlmbd = log(1. - 2.*lmbd);

    // LL + NLL Sudakov
    const complex<double> g1 = A1 * (2. - (2. - 1. / lmbd) * lnlmbd);
    const complex<double> g2 = -A2 * (2.*lmbd + lnlmbd)
                               + A1 * (2.*lmbd + lnlmbd) * log(S / Mur2) - A1 * 2.*lmbd * log(S / Muf2)
                               + A1 * beta1 / pow(beta0, 2) * (2.*lmbd + lnlmbd * (1. + .5 * lnlmbd));

    return exp(LL * g1 + g2);
}

complex<double> Thadronic_xs(complex<double> NM, double S,
                             double T, Parameters *Params)
{
    // Majorana symmetry factor
    double dij = 1.;
    if (Params->out1 < 10 && Params->out1 == Params->out2 && Params->out1 / 4 == 0) {
        dij = .5;
    }

    // Getting N-space PDFs
    complex<double> gN, qN[2][6];
    pdfN(gN, qN, NM, Params->afit);
    const complex<double> lnNb = log(NM) + M_EULER;

    const double aa = aS(Params->murs, Params->set) * 4. / 3.;
    //complex<double> sig(0., 0.);
    //complex<double> brn(0., 0.);
    complex<double> RES(0., 0.);
    complex<double> EXP(0., 0.);
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            Params->in1 = i0;
            Params->in2 = i1;
            // Test charge conservation
            if ((Params->out1 < 10 && Params->in1 / 3 - Params->in2 / 3 == Params->out1 / 4 - Params->out2 / 4)
                    || (Params->out1 >= 10 && Params->out1 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 10) / 3 - (Params->out2 - 10) / 3)


                    || (Params->out1 >= 16 && Params->out1 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 13) / 3 - (Params->out2 - 10) / 3)


                    || (Params->out2 >= 16 && Params->out2 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 13) / 3 - (Params->out2 - 13) / 3)) {


                double brn = Born(S, T, Params);
                if (brn == 0.) {
                    continue;
                }
                double vir = pow(M_PI, 2) / 3.*brn + Virt2(S, T, Params) + DipI(S, T, Params);
                RES += (qN[0][i0] * qN[1 - Params->ic][i1] + qN[Params->ic][i0] * qN[1][i1]) *
                       brn * exp((
                                     +vir / brn + log(S / Params->mufs) * 3.
                                     + (lnNb - .5 * log(S / Params->mufs)) / NM * 4.
                                 ) * aa);
                EXP += (qN[0][i0] * qN[1 - Params->ic][i1] + qN[Params->ic][i0] * qN[1][i1]) *
                       (brn + vir * aa + (log(S / Params->mufs) * 3.
                                          + (lnNb - .5 * log(S / Params->mufs)) / NM * 4.
                                          + lnNb * (lnNb - log(S / Params->mufs)) * 4.) * brn * aa);

            }
        }
    }

    return (RES * ThG(lnNb, S, Params->murs, Params->murs, Params->set) - EXP) * dij / 24. / S;
}

complex<double> ThG2(const complex<double> LL,
                     double MIS, double MURS, const int set)
{
    // Coeff. of the beta function
    const int nf = NFLVR;
    const double dnf = (double)nf;
    const double beta0 = 5.5 - dnf / 3.;
    const double beta1 = 25.5 - 9.5 * dnf / 3.;

    // Sudakov coeff. for quarks
    const double A1 = 8. / 3. / beta0;
    const double B1 = -4. / beta0;
    const double A2 = .5 * A1 * (67. / 3. - pow(M_PI, 2) - 10. / 9.*dnf) / beta0;

    // Logarithms
    const complex<double> lmbd = aS(MURS, set) * beta0 * LL;
    const complex<double> lnlmbd = log(1. - lmbd);

    // LL + NLL Sudakov
    const complex<double> g1 = A1 * (1. + lnlmbd / lmbd);
    const complex<double> g2 = (lmbd + lnlmbd)
                               * (A1 * (log(MIS / MURS) + beta1 / pow(beta0, 2)) - A2 + B1)
                               - B1 * lmbd + .5 * A1 * beta1 * pow(lnlmbd / beta0, 2);

    return exp(LL * g1 + g2);
}

complex<double> Thadronic_xs2(complex<double> N, double S, double T, Parameters *Params)
{
    // Logarithm
    const complex<double> nb = N * exp(M_EULER);
    const complex<double> ll = log(nb) * 2.;

    const complex<double> Pqq = 4. / 3.*(1.5 + 1. / N / (N + 1.) - 2.*(Psi(N + 1.) + M_EULER));
    const complex<double> Pqg = .5 * (2. + N + pow(N, 2)) / N / (N + 1.) / (N + 2.);

    // Majorana symmetry factor
    double dij = 1.;
    if (Params->out1 < 10 && Params->out1 == Params->out2 && Params->out1 / 4 == 0) {
        dij = .5;
    }

    // Getting N-space PDFs
    complex<double> gE, qE[2][6], gN, qN[2][6];
    pdfN(gN, qN, N, Params->afit);
    pdfN(gE, qE, N, Params->afit);

    // Evolving factorization scale (NLL)
    const int nf = NFLVR;
    const double dnf = (double)nf;
    const double beta0 = 5.5 - dnf / 3.;
    const double aa = aS(Params->murs, Params->set);
    const complex<double> lmbd = aa * beta0 * log(Params->mufs / S * pow(nb, 2));
    pdfEvolve(gE, qE, N, lmbd);

    complex<double> RES(0., 0.);
    complex<double> EXP(0., 0.);
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            Params->in1 = i0;
            Params->in2 = i1;
            // Test charge conservation
            if ((Params->out1 < 10 && Params->in1 / 3 - Params->in2 / 3 == Params->out1 / 4 - Params->out2 / 4)
                    || (Params->out1 >= 10 && Params->out1 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 10) / 3 - (Params->out2 - 10) / 3)


                    || (Params->out1 >= 16 && Params->out1 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 13) / 3 - (Params->out2 - 10) / 3)


                    || (Params->out2 >= 16 && Params->out2 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 13) / 3 - (Params->out2 - 13) / 3)) {
                double brn = Born(S, T, Params);
                double vir = pow(M_PI, 2) / 3.*brn + Virt(S, T, Params) + DipI(S, T, Params);
                RES += (qE[0][i0] * qE[1 - Params->ic][i1] + qE[Params->ic][i0] * qE[1][i1]) * (brn + vir * aa * 4. / 3.);
                EXP += (qN[0][i0] * qN[1 - Params->ic][i1] + qN[Params->ic][i0] * qN[1][i1]) * (brn + vir * 4. / 3.*aa +
                        (ll * (-ll * 4. / 3. + 4.) - 2.*Pqq * (ll - log(S / Params->mufs))) * aa * brn);
                EXP += (qN[0][i0] + qN[1 - Params->ic][i1] + qN[Params->ic][i0] + qN[1][i1]) * gN *
                       (log(S / Params->mufs) - ll) * Pqg * aa * brn;
            }
        }
    }
    // Flux, symmetry factor, spin and color average
    return (RES * ThG2(ll, S, Params->murs, Params->set) - EXP) * dij / 24. / S;
}

// Transverse momentum resummation.
void CFact(complex<double> &g, complex<double> q[2][6],
           complex<double> N, double aSmu, complex<double> aSbb)
{
    complex<double> gdum = g;
    complex<double> qdum[2][6];
    for (size_t i0 = 0; i0 < 2; i0++) {
        for (size_t i1 = 0; i1 < 6; i1++) {
            qdum[i0][i1] = q[i0][i1];
        }
    }

    complex<double> Cqq, Cqg;
    Cqq = 1. + aSmu * 4. / 3. / N / (N + 1.); // aSmu at NLL / aSbb at NNLL
    Cqg = aSbb / (N + 1.) / (N + 2.); // aSbb at NLL / aSbb at NNLL

    for (size_t i0 = 0; i0 < 2; i0++) {
        for (size_t i1 = 0; i1 < 6; i1++) {
            q[i0][i1] = Cqq * qdum[i0][i1] + Cqg * gdum;
        }
    }
}

complex<double> PtG(const complex<double> LL,
                    double MIS, double MURS, const int set)
{
    // Coeff. of the beta function
    const int nf = NFLVR;
    const double dnf = (double)nf;
    const double beta0 = 5.5 - dnf / 3.;
    const double beta1 = 25.5 - 9.5 * dnf / 3.;

    // Sudakov coeff. for quarks
    const double A1 = 8. / 3. / beta0;
    const double B1 = -4. / beta0;
    const double A2 = .5 * A1 * (67. / 3. - pow(M_PI, 2) - 10. / 9.*dnf) / beta0;

    // Logarithms
    const complex<double> lmbd = aS(MURS, set) * beta0 * LL;
    const complex<double> lnlmbd = log(1. - lmbd);

    // LL + NLL Sudakov
    const complex<double> g1 = A1 * (1. + lnlmbd / lmbd);
    const complex<double> g2 = (lmbd / (1. - lmbd) + lnlmbd)
                               * (A1 * (log(MIS / MURS) + beta1 / pow(beta0, 2) * (1. + lnlmbd)) - A2 + B1)
                               - B1 * lmbd / (1. - lmbd) - .5 * A1 * beta1 * pow(lnlmbd / beta0, 2);

    return exp(LL * g1 + g2);
}

complex<double> PtXS(complex<double> N, complex<double> B,
                     double S, double T, Parameters *Params)
{
    // Logarithm
    const complex<double> bb = B * sqrt(S) * exp(M_EULER) * .5;
    const complex<double> ll = log(1. + bb * bb);

    const complex<double> Pqq =
        4. / 3.*(1.5 + 1. / N / (N + 1.) - 2.*(Psi(N + 1.) + M_EULER));
    const complex<double> Pqg =
        .5 * (2. + N + pow(N, 2)) / N / (N + 1.) / (N + 2.);

    // Majorana symmetry factor
    double dij = 1.;
    if (Params->out1 < 10 && Params->out1 == Params->out2 && Params->out1 / 4 == 0) {
        dij = .5;
    }

    // Getting N-space PDFs
    complex<double> gE, qE[2][6], gN, qN[2][6];
    pdfN(gN, qN, N, Params->afit);
    pdfN(gE, qE, N, Params->afit);

    // Evolving factorization scale [ NLL ]
    const int nf = NFLVR;
    const double dnf = (double)nf;
    const double beta0 = 5.5 - dnf / 3.;
    const double aa = aS(Params->murs, Params->set);
    const complex<double> lmbd = aa * beta0 * log(Params->mufs / S * (1. + bb * bb));
    pdfEvolve(gE, qE, N, lmbd);
    CFact(gE, qE, N, aa, aa / (1. - lmbd));

    complex<double> RES(0., 0.);
    complex<double> EXP(0., 0.);
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            Params->in1 = i0;
            Params->in2 = i1;
            // Test charge conservation
            if ((Params->out1 < 10 && Params->in1 / 3 - Params->in2 / 3 == Params->out1 / 4 - Params->out2 / 4)
                    || (Params->out1 >= 10 && Params->out1 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 10) / 3 - (Params->out2 - 10) / 3)


                    || (Params->out1 >= 16 && Params->out1 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 13) / 3 - (Params->out2 - 10) / 3)


                    || (Params->out2 >= 16 && Params->out2 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 13) / 3 - (Params->out2 - 13) / 3)) {

                double brn = Born(S, T, Params);
                double vir = Virt(S, T, Params) + DipI(S, T, Params);
                RES += (qE[0][i0] * qE[1 - Params->ic][i1] + qE[Params->ic][i0] * qE[1][i1]) *
                       (brn + vir * aa * 4. / 3.);
                EXP += (qN[0][i0] * qN[1 - Params->ic][i1] + qN[Params->ic][i0] * qN[1][i1]) *
                       ll * (-ll * 4. / 3. + 4. - 2.*Pqq) * aa * brn;
                EXP += (qN[0][i0] + qN[1 - Params->ic][i1] + qN[Params->ic][i0] + qN[1][i1]) * gN *
                       (-ll) * Pqg * aa * brn;
            }
        }
    }
    // Flux, symmetry factor, spin and color average
    return (RES * PtG(ll, S, Params->murs, Params->set) - EXP) * dij / 24. / S;
}

// Joint resummation
complex<double> JtG(const complex<double> LL,
                    const complex<double> LN, double MIS, double MURS, const int set)
{
    // Coeff. of the beta function
    const int nf = NFLVR;
    const double dnf = (double)nf;
    const double beta0 = 5.5 - dnf / 3.;
    const double beta1 = (25.5 - 9.5 * dnf / 3.) / pow(beta0, 2);

    // Sudakov coeff. for quarks
    const double A1 = 8. / 3. / beta0;
    const double B1 = -4. / beta0;
    const double A2 = .5 * A1 * (67. / 3. - pow(M_PI, 2) - 10. / 9.*dnf) / beta0;

    // Logarithms
    const complex<double> lmbdN = aS(MURS, set) * beta0 * LN;
    const complex<double> lmbd = aS(MURS, set) * beta0 * LL;
    const complex<double> lnlmbd = log(1. - lmbd);

    // LL + NLL Sudakov
    const complex<double> g1 = A1 * (1. + lnlmbd / lmbd);
    const complex<double> g2 = (lmbd / (1. - lmbd) * (1. - lmbdN) + lnlmbd)
                               * (A1 * log(MIS / MURS) - A2) + A1 * beta1 * (.5 * pow(lnlmbd, 2) +
                                       (1. - lmbdN) / (1. - lmbd) * (lmbd + lnlmbd)) + B1 * lnlmbd;

    return exp(LL * g1 + g2);
}

complex<double> JtXS(complex<double> N, complex<double> B,
                     double S, double T, Parameters *Params)
{
    // Logarithm
    const double eta = 1.;
    const complex<double> nb = N * exp(M_EULER);
    const complex<double> bb = B * sqrt(S) * exp(M_EULER) * .5;
    const complex<double> ll = log(bb + nb / (1. + bb / nb * eta)) * 2.;
    const complex<double> ln = log(nb) * 2.;

    const complex<double> Pqq =
        4. / 3.*(1.5 + 1. / N / (N + 1.) - 2.*(Psi(N + 1.) + M_EULER));
    const complex<double> Pqg =
        .5 * (2. + N + pow(N, 2)) / N / (N + 1.) / (N + 2.);

    // Majorana symmetry factor
    double dij = 1.;
    if (Params->out1 < 10 && Params->out1 == Params->out2 && Params->out1 / 4 == 0) {
        dij = .5;
    }

    // Getting N-space PDFs
    complex<double> gE, qE[2][6], gN, qN[2][6];
    pdfN(gN, qN, N, Params->afit);
    pdfN(gE, qE, N, Params->afit);

    // Evolving factorization scale [ NLL ]
    const int nf = NFLVR;
    const double dnf = (double)nf;
    const double beta0 = 5.5 - dnf / 3.;
    const double aa = aS(Params->murs, Params->set);
    const complex<double> lmbd = aa * beta0 *
                                 log(Params->mufs / S * pow(bb + nb / (1. + bb / nb * eta), 2));
    pdfEvolve(gE, qE, N, lmbd);
    CFact(gE, qE, N, aa, aa / (1. - lmbd));

    complex<double> RES(0., 0.);
    complex<double> EXP(0., 0.);
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            Params->in1 = i0;
            Params->in2 = i1;
            // Test charge conservation
            if ((Params->out1 < 10 && Params->in1 / 3 - Params->in2 / 3 == Params->out1 / 4 - Params->out2 / 4)
                    || (Params->out1 >= 10 && Params->out1 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 10) / 3 - (Params->out2 - 10) / 3)


                    || (Params->out1 >= 16 && Params->out1 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 13) / 3 - (Params->out2 - 10) / 3)


                    || (Params->out2 >= 16 && Params->out2 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 13) / 3 - (Params->out2 - 13) / 3)) {
                double brn = Born(S, T, Params);
                double vir = Virt(S, T, Params) + DipI(S, T, Params);
                RES += (qE[0][i0] * qE[1 - Params->ic][i1] + qE[Params->ic][i0] * qE[1][i1]) *
                       (brn + vir * aa * 4. / 3.);
                EXP += (qN[0][i0] * qN[1 - Params->ic][i1] + qN[Params->ic][i0] * qN[1][i1]) *
                       ll * (-ll * 4. / 3. + 4. - 2.*Pqq) * aa * brn;
                EXP += (qN[0][i0] + qN[1 - Params->ic][i1] + qN[Params->ic][i0] + qN[1][i1]) * gN *
                       (-ll) * Pqg * aa * brn;
            }
        }
    }
    // Flux, symmetry factor, spin and color average
    return (RES * JtG(ll, ln, S, Params->murs, Params->set) - EXP) * dij / 24. / S;
}

complex<double> JtXS2(complex<double> N, complex<double> B,
                      double S, double T, Parameters *Params)
{
    // Logarithm
    const double eta = 1.;
    const complex<double> nb = N * exp(M_EULER);
    const complex<double> bb = B * sqrt(S) * exp(M_EULER) * .5;
    const complex<double> ll = log(bb + nb / (1. + bb / nb * eta)) * 2.;
    const complex<double> ln = log(nb) * 2.;

    const complex<double> Pqq =
        4. / 3.*(1.5 + 1. / N / (N + 1.) - 2.*(Psi(N + 1.) + M_EULER));
    const complex<double> Pqg =
        .5 * (2. + N + pow(N, 2)) / N / (N + 1.) / (N + 2.);

    // Majorana symmetry factor
    double dij = 1.;
    if (Params->out1 < 10 && Params->out1 == Params->out2 && Params->out1 / 4 == 0) {
        dij = .5;
    }

    // Getting N-space PDFs
    complex<double> gE, qE[2][6], gN, qN[2][6];
    pdfN(gN, qN, N, Params->afit);
    pdfN(gE, qE, N, Params->afit);

    // Evolving factorization scale [ NLL ]
    const int nf = NFLVR;
    const double dnf = (double)nf;
    const double beta0 = 5.5 - dnf / 3.;
    const double aa = aS(Params->murs, Params->set);
    const complex<double> lmbd = aa * beta0 *
                                 log(Params->mufs / S * pow(bb + nb / (1. + bb / nb * eta), 2));
    pdfEvolve(gE, qE, N, lmbd);
    CFact(gE, qE, N, aa, aa / (1. - lmbd));

    complex<double> RES(0., 0.);
    complex<double> EXP(0., 0.);
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            Params->in1 = i0;
            Params->in2 = i1;
            // Test charge conservation
            if ((Params->out1 < 10 && Params->in1 / 3 - Params->in2 / 3 == Params->out1 / 4 - Params->out2 / 4)
                    || (Params->out1 >= 10 && Params->out1 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 10) / 3 - (Params->out2 - 10) / 3)


                    || (Params->out1 >= 16 && Params->out1 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 13) / 3 - (Params->out2 - 10) / 3)


                    || (Params->out2 >= 16 && Params->out2 < 20
                        && Params->in1 / 3 - Params->in2 / 3 == (Params->out1 - 13) / 3 - (Params->out2 - 13) / 3)) {

                double brn = Born(S, T, Params);
                double vir = Virt(S, T, Params) + DipI(S, T, Params);
                RES += (qE[0][i0] * qE[1 - Params->ic][i1] + qE[Params->ic][i0] * qE[1][i1]) *
                       (brn + vir * aa * 4. / 3.);
                EXP += (qN[0][i0] * qN[1 - Params->ic][i1] + qN[Params->ic][i0] * qN[1][i1]) *
                       (brn + vir * 4. / 3.*aa +
                        (ll * (-ll * 4. / 3. + 4.) - 2.*Pqq * (ll - log(S / Params->mufs))) * aa * brn);
                EXP += (qN[0][i0] + qN[1 - Params->ic][i1] + qN[Params->ic][i0] + qN[1][i1]) * gN *
                       (log(S / Params->mufs) - ll) * Pqg * aa * brn;
            }
        }
    }
    // Flux, symmetry factor, spin and color average
    return (RES * JtG(ll, ln, S, Params->murs, Params->set) - EXP) * dij / 24. / S;
}
