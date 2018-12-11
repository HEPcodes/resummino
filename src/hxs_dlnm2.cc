// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Computes the hadronic cross section invariant-mass distribution at LO,
// NLO (collinear, virtual, gluon emission and quark emission) and NLL.

#include <cmath>
#include <complex>
#include <iostream>
#include <cstdlib>
#include "gsl_all.h"
#include "prm.h"
#include "utils.h"
#include "pxs.h"
#include "pdf.h"

#define CMLLN 0.9 // C coeff. of the inverse Mellin transform
#define VBSSL 2.9 // V coeff. of the inverse Bessel transform
#define PT2CUT 0.0 // pt2 cut for the integration (non-zero for joint resummation)

static inline double kln(double x, double y, double z)
{
    return std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)
                     - 2.*x * y - 2.*y * z - 2.*z * x);
}

double dPS2_dlnM2(double &xa, double &xb, double &t, double *x, Parameters *Params)
{
    const double sh = Params->sh;
    //old: const double m1 = Params->mCH[Params->out1];
    //old: const double m2 = Params->mCH[Params->out2];

    double m1;
    double m2;

    if (Params->out1 < 10) {
        m1 = Params->mCH[Params->out1];
        m2 = Params->mCH[Params->out2];
    }

    else  if (Params->out1 >= 10 && Params->out1 < 20) {
        m1 = Params->mSL[Params->out1 - 10];
        m2 = Params->mSL[Params->out2 - 10];
    }


    const double m1s = std::pow(m1, 2);
    const double m2s = std::pow(m2, 2);
    const double s = Params->mis;

    if (s < std::pow(m1 + m2, 2)) {
        std::cout << "dPS2_dlnM2: M2 too small\n";
        exit(0);
    } else if (s > sh) {
        std::cout << "dPS2_dlnM2: M2 too large\n";
        exit(1);
    }

    // Jacobian initialization
    double djac = 389379304.;

    // Integration variable xa
    double xamin = s / sh;
    double xamax = 1.;
    if (xamin == 0.) {
        xa = (xamax - xamin) * x[0] + xamin;
        djac *= xamax - xamin;
    } else {
        xa = xamin * std::pow(xamax / xamin, x[0]);
        djac *= xa * std::log(xamax / xamin);
    }

    // Getting xb
    xb = xamin / xa;

    // Integration variable t
    const double tmin = -.5 * (s - m1s - m2s + kln(s, m1s, m2s));
    const double tmax = tmin + kln(s, m1s, m2s);
    t = (tmax - tmin) * x[1] + tmin;
    djac *= tmax - tmin;

    // Phase Space factor
    djac *= 0.125 * M_1_PI / s * xb;

    return djac;
}

double dPS2C_dlnM2(double &djacdelta, double &djacplus,
                   double &xa, double &xb, double &xc,
                   double &tc, double *x, Parameters *Params)
{
    const double sh = Params->sh;
    //old:  const double m1 = Params->mCH[Params->out1];
    //old:  const double m2 = Params->mCH[Params->out2];

    double m1;
    double m2;

    if (Params->out1 < 10) {
        m1 = Params->mCH[Params->out1];
        m2 = Params->mCH[Params->out2];
    }

    else  if (Params->out1 >= 10 && Params->out1 < 20) {
        m1 = Params->mSL[Params->out1 - 10];
        m2 = Params->mSL[Params->out2 - 10];
    }

    const double m1s = std::pow(m1, 2);
    const double m2s = std::pow(m2, 2);
    const double mis = Params->mis;

    if (mis < std::pow(m1 + m2, 2)) {
        std::cout << "dPS2d_dlnM2: M2 too small\n";
        exit(0);
    } else if (mis > sh) {
        std::cout << "dPS2d_dlnM2: M2 too large\n";
        exit(1);
    }

    // Jacobian initialization
    double djac = 389379304.;

    // Integration variable xa and xb
    double xamin = mis / sh;
    double xamax = 1.;
    double xcmin = xamin;
    double xcmax = 1.;
    if (xamin == 0.) {
        xa = (xamax - xamin) * x[0] + xamin;
        djac *= xamax - xamin;
        xcmin /= xa;
        djacdelta = djac;
        xc = (xcmax - xcmin) * x[1] + xcmin;
        djac *= xcmax - xcmin;
    } else {
        xa = xamin * std::pow(xamax / xamin, x[0]);
        djac *= xa * std::log(xamax / xamin);
        xcmin /= xa;
        djacdelta = djac;
        xc = xcmin * std::pow(xcmax / xcmin, x[1]);
        djac *= xc * std::log(xcmax / xcmin);
    }

    xb = xamin / xa / xc;

    const double sc = mis;

    // Integration variable tc
    const double tcmin = -.5 * (sc - m1s - m2s + kln(sc, m1s, m2s));
    const double tcmax = tcmin + kln(sc, m1s, m2s);
    tc = (tcmax - tcmin) * x[2] + tcmin;
    djacdelta *= tcmax - tcmin;
    djac     *= tcmax - tcmin;

    // Phase Space factor
    djacdelta *= 0.125 * M_1_PI / sc * xb * xc;
    djac     *= 0.125 * M_1_PI / sc * xb;
    djacplus  = djac * xc;

    return djac;
}

double dPS3_dlnM2(double &xa, double &xb, double &M2, double &pt2,
                  double &Tp, double &Pp, double *x, Parameters *Params)
{
    const double sh = Params->sh;
    //old: const double m1 = Params->mCH[Params->out1];
    //old: const double m2 = Params->mCH[Params->out2];

    double m1;
    double m2;

    if (Params->out1 < 10) {
        m1 = Params->mCH[Params->out1];
        m2 = Params->mCH[Params->out2];
    }

    else  if (Params->out1 >= 10 && Params->out1 < 20) {
        m1 = Params->mSL[Params->out1 - 10];
        m2 = Params->mSL[Params->out2 - 10];
    }

    const double m1s = std::pow(m1, 2);
    const double m2s = std::pow(m2, 2);
    M2 = Params->mis;

    if (M2 < std::pow(m1 + m2, 2)) {
        std::cout << "dPS3_dlnM2: M2 too small\n";
        exit(0);
    } else if (M2 > sh) {
        std::cout << "dPS3_dlnM2: M2 too large\n";
        exit(0);
    }

    // Jacobian initialization
    double djac = 389379304.;

    // Integration variable xa, xb
    double xamin = M2 / sh;
    double xamax = 1.;
    double xbmin = xamin;
    double xbmax = 1.;
    if (xamin == 0.) {
        xa = (xamax - xamin) * x[0] + xamin;
        djac *= xamax - xamin;
        xbmin /= xa;
        xb = (xbmax - xbmin) * x[1] + xbmin;
        djac *= xbmax - xbmin;
    } else {
        xa = xamin * std::pow(xamax / xamin, x[0]);
        djac *= xa * std::log(xamax / xamin);
        xbmin /= xa;
        xb = xbmin * std::pow(xbmax / xbmin, x[1]);
        djac *= xb * std::log(xbmax / xbmin);
    }

    const double s = xa * xb * sh;

    // Integration variable pt2
    const double pt2min = PT2CUT;
    const double pt2max = .25 * std::pow(s - M2, 2) / s;
    if (pt2max < pt2min) {
        return 0.;
    }
    pt2 = (pt2max - pt2min) * x[2] + pt2min;
    djac *= pt2max - pt2min;

    // Integration variable Tp
    const double Tpmin = 0.;
    const double Tpmax = M_PI;
    Tp = (Tpmax - Tpmin) * x[3] + Tpmin;
    djac *= Tpmax - Tpmin;

    // Integration variable Pp
    const double Ppmin = 0.;
    const double Ppmax = M_PI;
    djac *= 2.;
    Pp = (Ppmax - Ppmin) * x[4] + Ppmin;
    djac *= Ppmax - Ppmin;

    // Phase Space factor
    djac *= kln(M2, m1s, m2s) * std::sin(Tp)
            / std::sqrt(std::pow(s - M2, 2) - 4.*s * pt2) / std::pow(4.*M_PI, 4);

    return djac;
}

double IB_dlnM2(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj; // conversion of type void* into IOS*
    double xa, xb, t;
    const double djac = dPS2_dlnM2(xa, xb, t, x, jp);
    const double s = xa * xb * jp->sh;

    double dij = 1.;
    if (jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
        dij = .5;
    }

    double ga, qa[2][6];
    pdfX(ga, qa, xa, jp->mufs);
    double gb, qb[2][6];
    pdfX(gb, qb, xb, jp->mufs);

    double sig = 0.;
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            jp->in1 = i0;
            jp->in2 = i1;
            // Test charge conservation
            //old:     if (jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
            if ((jp->out1 < 10 && jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
                    || (jp->out1 >= 10 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 10) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out1 >= 16 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out2 >= 16 && jp->out2 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 13) / 3)) {
                sig += Born(s, t, jp) *
                       (qa[0][i0] * qb[1 - jp->ic][i1] + qb[jp->ic][i0] * qa[1][i1]);
            }
        }
    }
    // Flux, symmetry factor, spin and color average
    return dij / 24. / s * sig * djac;
}

double IV_dlnM2(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj; // conversion of type void* into IOS*
    double xa, xb, t;
    const double djac = dPS2_dlnM2(xa, xb, t, x, jp);
    const double s = xa * xb * jp->sh;

    double dij = 1.;
    if (jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
        dij = .5;
    }

    double ga, qa[2][6];
    pdfX(ga, qa, xa, jp->mufs);
    double gb, qb[2][6];
    pdfX(gb, qb, xb, jp->mufs);

    double sig = 0.;
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            jp->in1 = i0;
            jp->in2 = i1;
            // Test charge conservation
            //old:    if (jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
            if ((jp->out1 < 10 && jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
                    || (jp->out1 >= 10 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 10) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out1 >= 16 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out2 >= 16 && jp->out2 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 13) / 3)) {
                sig += (qa[0][i0] * qb[1 - jp->ic][i1] + qb[jp->ic][i0] * qa[1][i1])
                       * (Virt(s, t, jp) + DipI(s, t, jp));
            }
        }
    }
    // Flux, symmetry factor, spin and color average
    return dij / 24. / s * sig * djac * 4. / 3.*aS(jp->murs, jp->set);
}


double IC_dlnM2(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj; // conversion of type void* into IOS*
    double djacdelta, djacplus, xa, xb, xc, tc;
    const double djac = dPS2C_dlnM2(djacdelta, djacplus, xa, xb, xc, tc, x, jp);
    const double sc = jp->mis;
    const double xd = xb * xc;

    double dij = 1.;
    if (jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
        dij = .5;
    }

    double ga, qa[2][6];
    pdfX(ga, qa, xa, jp->mufs);
    double gb, qb[2][6];
    pdfX(gb, qb, xb, jp->mufs);
    double gd, qd[2][6];
    pdfX(gd, qd, xd, jp->mufs);

    double sig = 0.;
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            jp->in1 = i0;
            jp->in2 = i1;
            // Test charge conservation
            //old:   if (jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4) {
            if ((jp->out1 < 10 && jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
                    || (jp->out1 >= 10 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 10) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out1 >= 16 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out2 >= 16 && jp->out2 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 13) / 3)) {
                // Pqq and Kqq
                sig += (qa[0][i0] * qb[1 - jp->ic][i1] + qb[jp->ic][i0] * qa[1][i1] +
                        qb[0][i0] * qa[1 - jp->ic][i1] + qa[jp->ic][i0] * qb[1][i1]) *
                       (-(1. + std::pow(xc, 2)) / (1. - xc) * std::log(jp->mufs / sc)
                        + 2.*(2. / (1. - xc) - 1. - xc) * std::log(1. - xc) + 1. - xc
                       ) * Born(sc, tc, jp) * djac * 4. / 3.;

                sig += (qa[0][i0] * qd[1 - jp->ic][i1] + qd[jp->ic][i0] * qa[1][i1] +
                        qd[0][i0] * qa[1 - jp->ic][i1] + qa[jp->ic][i0] * qd[1][i1]) *
                       (((1. + std::pow(xc, 2)) / (1. - xc) * std::log(jp->mufs / sc)
                         - 4. / (1. - xc) * std::log(1. - xc))
                       ) * Born(sc, tc, jp) * djacplus * 4. / 3.;

                sig -= (qa[0][i0] * qd[1 - jp->ic][i1] + qd[jp->ic][i0] * qa[1][i1]) *
                       std::pow(M_PI, 2) * Born(sc, tc, jp) * djacdelta * 4. / 9.;

                // Pgq, Kgq, Pgqb, Kgqb
                sig += (qa[0][i0] * gb + qb[jp->ic][i0] * ga +
                        qa[1][i1] * gb + qb[1 - jp->ic][i1] * ga) *
                       (-(std::pow(xc, 2) + std::pow(1. - xc, 2)) *
                        std::log(jp->mufs / sc * xc / std::pow(1. - xc, 2))
                        + 2.*xc * (1. - xc)) * Born(sc, tc, jp) * djac * .5;
            }
        }
    }

    // Flux, symmetry factor, spin and color average
    return dij / 24.*sig * aS(jp->murs, jp->set) / sc;
}

// No dipoles -> cf. addition to joint res.
double IG_dlnM2j(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj; // conversion of type void* into IOS*
    double xa, xb, mi2, pt2, th, ph;
    const double djac = dPS3_dlnM2(xa, xb, mi2, pt2, th, ph, x, jp);
    const double s = xa * xb * jp->sh;
    if (djac == 0.) {
        return 0.;
    }

    double dij = 1.;
    if (jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
        dij = .5;
    }

    double ga, qa[2][6];
    pdfX(ga, qa, xa, jp->mufs);
    double gb, qb[2][6];
    pdfX(gb, qb, xb, jp->mufs);

    double sig = 0.;
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            jp->in1 = i0;
            jp->in2 = i1;
            // Test charge conservation
            //old:   if (jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
            if ((jp->out1 < 10 && jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
                    || (jp->out1 >= 10 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 10) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out1 >= 16 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out2 >= 16 && jp->out2 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 13) / 3)) {
                sig += (qa[0][i0] * qb[1 - jp->ic][i1] + qb[jp->ic][i0] * qa[1][i1]) *
                       (RealG(s, mi2, pt2, th, ph, 0, jp) +
                        RealG(s, mi2, pt2, th, ph, 1, jp));
            }
        }
    }
    // Flux, symmetry factor, spin and color average
    return dij * std::pow(M_PI, 2) / 4.5 / s * sig * djac * aS(jp->murs, jp->set);

}

double IG_dlnM2(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj; // conversion of type void* into IOS*
    double xa, xb, mi2, pt2, th, ph;
    const double djac = dPS3_dlnM2(xa, xb, mi2, pt2, th, ph, x, jp);
    const double s = xa * xb * jp->sh;
    if (djac == 0.) {
        return 0.;
    }

    double dij = 1.;
    if (jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
        dij = .5;
    }

    double ga, qa[2][6];
    pdfX(ga, qa, xa, jp->mufs);
    double gb, qb[2][6];
    pdfX(gb, qb, xb, jp->mufs);

    double sig = 0.;
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            jp->in1 = i0;
            jp->in2 = i1;
            // Test charge conservation
            //old:  if (jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
            if ((jp->out1 < 10 && jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
                    || (jp->out1 >= 10 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 10) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out1 >= 16 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out2 >= 16 && jp->out2 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 13) / 3)) {
                sig += (qa[0][i0] * qb[1 - jp->ic][i1] + qb[jp->ic][i0] * qa[1][i1]) *
                       (RealG(s, mi2, pt2, th, ph, 0, jp) +
                        RealG(s, mi2, pt2, th, ph, 1, jp) -
                        DipGA(s, mi2, pt2, th, ph, 0, jp) -
                        DipGA(s, mi2, pt2, th, ph, 1, jp) -
                        DipGB(s, mi2, pt2, th, ph, 0, jp) -
                        DipGB(s, mi2, pt2, th, ph, 1, jp)
                       );
            }
        }
    }
    // Flux, symmetry factor, spin and color average
    return dij * std::pow(M_PI, 2) / 4.5 / s * sig * djac * aS(jp->murs, jp->set);

}

double IQ_dlnM2(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj; // conversion of type void* into IOS*
    double xa, xb, mi2, pt2, th, ph;
    const double djac = dPS3_dlnM2(xa, xb, mi2, pt2, th, ph, x, jp);
    const double s = xa * xb * jp->sh;
    if (djac == 0.) {
        return 0.;
    }

    double dij = 1.;
    if (jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
        dij = .5;
    }

    double ga, qa[2][6];
    pdfX(ga, qa, xa, jp->mufs);
    double gb, qb[2][6];
    pdfX(gb, qb, xb, jp->mufs);

    double sig = 0.;
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            jp->in1 = i0;
            jp->in2 = i1;
            // Test charge conservation
            //old:    if (jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4) {
            if ((jp->out1 < 10 && jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
                    || (jp->out1 >= 10 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 10) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out1 >= 16 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out2 >= 16 && jp->out2 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 13) / 3)) {
                sig += (qa[0][i0] * gb + qb[jp->ic][i0] * ga) *
                       (RealQ(s, mi2, pt2, th, ph, 0, jp) +
                        RealQ(s, mi2, pt2, th, ph, 1, jp) -
                        DipQB(s, mi2, pt2, th, ph, 0, jp) -
                        DipQB(s, mi2, pt2, th, ph, 1, jp));
                sig += (qa[1][i1] * gb + qb[1 - jp->ic][i1] * ga) *
                       (RealQB(s, mi2, pt2, th, ph, 0, jp) +
                        RealQB(s, mi2, pt2, th, ph, 1, jp) -
                        DipQA(s, mi2, pt2, th, ph, 0, jp) -
                        DipQA(s, mi2, pt2, th, ph, 1, jp));
            }
        }
    }

    // Flux, symmetry factor, spin and color average
    return dij * std::pow(M_PI, 2) / 12. / s * sig * djac * aS(jp->murs, jp->set);
}

// No dipoles -> cf. addition to the joint
double IQ_dlnM2j(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj; // conversion of type void* into IOS*
    double xa, xb, mi2, pt2, th, ph;
    const double djac = dPS3_dlnM2(xa, xb, mi2, pt2, th, ph, x, jp);
    const double s = xa * xb * jp->sh;
    if (djac == 0.) {
        return 0.;
    }
    double dij = 1.;
    if (jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
        dij = .5;
    }

    double ga, qa[2][6];
    pdfX(ga, qa, xa, jp->mufs);
    double gb, qb[2][6];
    pdfX(gb, qb, xb, jp->mufs);

    double sig = 0.;
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            jp->in1 = i0;
            jp->in2 = i1;
            // Test charge conservation
            //old:  if (jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4) {
            if ((jp->out1 < 10 && jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
                    || (jp->out1 >= 10 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 10) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out1 >= 16 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out2 >= 16 && jp->out2 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 13) / 3)) {
                sig += (qa[0][i0] * gb + qb[jp->ic][i0] * ga) *
                       (RealQ(s, mi2, pt2, th, ph, 0, jp) +
                        RealQ(s, mi2, pt2, th, ph, 1, jp));
                sig += (qa[1][i1] * gb + qb[1 - jp->ic][i1] * ga) *
                       (RealQB(s, mi2, pt2, th, ph, 0, jp) +
                        RealQB(s, mi2, pt2, th, ph, 1, jp));
            }
        }
    }

    // Flux, symmetry factor, spin and color average
    return dij * std::pow(M_PI, 2) / 12. / s * sig * djac * aS(jp->murs, jp->set);
}

double IR_dlnM2(double *x, size_t dim, void *prm)
{
    Parameters *jp = (Parameters *)prm; // conversion of type void* into IOS*

    const double sh = jp->sh;
    //old: const double m1 = jp->mCH[jp->out1];
    //old: const double m2 = jp->mCH[jp->out2];
    double m1;
    double m2;

    if (jp->out1 < 10) {
        m1 = jp->mCH[jp->out1];
        m2 = jp->mCH[jp->out2];
    }

    else  if (jp->out1 >= 10 && jp->out1 < 20) {
        m1 = jp->mSL[jp->out1 - 10];
        m2 = jp->mSL[jp->out2 - 10];
    }

    const double m1s = std::pow(m1, 2);
    const double m2s = std::pow(m2, 2);
    double s = jp->mis;

    if (s < std::pow(m1 + m2, 2)) {
        std::cout << "IvMlln: M2 too small\n";
        exit(0);
    } else if (s > sh) {
        std::cout << "IvMlln: M2 too large\n";
        exit(1);
    }

    // Jacobian initialization
    std::complex<double> djac(389379304., 0.);

    // Inverse Mellin
    const std::complex<double> ephi(-M_SQRT1_2, M_SQRT1_2);
    const std::complex<double> nm = (s / sh / .09) + .09 - jp->a1min - ephi * std::log(x[0]);
    djac *= ephi / x[0];

    // t integration
    const double tmin = -.5 * (s - m1s - m2s + kln(s, m1s, m2s));
    const double tmax = tmin + kln(s, m1s, m2s);
    double t = (tmax - tmin) * x[1] + tmin;
    djac *= tmax - tmin;

    djac /= 8.*M_PI * s;

    return M_1_PI * std::imag(djac * std::pow(s / sh, -nm + 1.) * Thadronic_xs2(nm, s, t, jp));
}

void hadronic_xs_dlnM2(double &res, double &err, double &chi2,
                       int Flag, int Verb, Parameters *Params)
{

    // Definition of the integral
    size_t calls = 0;
    gsl_monte_function I;
    switch (Flag) {
    case 0:
        I.f = &IB_dlnM2;
        I.dim = 2;
        I.params = Params;
        calls = 10000;
        break;
    case 1:
        I.f = &IV_dlnM2;
        I.dim = 2;
        I.params = Params;
        calls = 1000;
        break;
    case 2:
        I.f = &IC_dlnM2;
        I.dim = 3;
        I.params = Params;
        calls = 15000;
        break;
    case 3:
        I.f = &IG_dlnM2;
        I.dim = 5;
        I.params = Params;
        calls = 25000;
        break;
    case 4:
        I.f = &IQ_dlnM2;
        I.dim = 5;
        I.params = Params;
        calls = 25000;
        break;
    case 5:
        I.f = &IR_dlnM2;
        I.dim = 2;
        I.params = Params;
        calls = 20000;
        break;
        //case 6: I.f=&IJ_dlnM2; I.dim=5; I.params=Params; calls=50000; break; // Bugged ?
    default:
        std::cout << "hadronic_xs: Flag=" << Flag << std::endl;
        exit(0);
    }
    const size_t dnum = I.dim;

    double xmin[dnum], xmax[dnum];
    for (size_t i0 = 0; i0 < dnum; i0++) {
        xmin[i0] = 0.;
        xmax[i0] = 1.;
    }

    // Initialization
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dnum);
    s->ostream = Params->fout;
    s->verbose = Verb;

    // Integration warm-up.
    s->stage = 0;
    s->iterations = 5;
    gsl_monte_vegas_integrate(&I, xmin, xmax, dnum, calls / 10, r, s,
                              &res, &err);

    // Integrates.
    s->stage = 1;
    s->iterations = 5;
    gsl_monte_vegas_integrate(&I, xmin, xmax, dnum, calls, r, s, &res, &err);

    // Convergence step.
    double prec;
    int counter = 1;

    do {

      if (abs(res) < 1.0e-15) {
        break;
      }
      
      s->iterations = 1;
      s->stage = 3;

        gsl_monte_vegas_integrate(&I, xmin, xmax, dnum, calls / 5, r, s,
                                  &res, &err);
        prec = abs(err / res);
        counter++;
    } while (prec > Params->precision && counter <= Params->max_iters);

    err = pow2(err);
    chi2 = s->chisq;

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);
}
