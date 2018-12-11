// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2012 David R. Lamprea.
// Copyright 2011-2012 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Computes the hadronic cross section at LO, NLO (collinear, virtual,
// gluon emission and quark emission) and NLL.

#include <cmath>
#include <complex>
#include <iostream>
#include <cstdlib>
#include "utils.h"
#include "gsl_all.h"
#include "prm.h"
#include "pxs.h"
#include "pdf.h"


using namespace std;

double* get_masses(Parameters *Params)
{
    const double sh = Params->sh;
    double m1;
    double m2;
    double *masses  = (double*)malloc(2 * sizeof(double));

    if (Params->out1 < 10) {
        m1 = Params->mCH[Params->out1];
        m2 = Params->mCH[Params->out2];
    }

    else  if (Params->out1 >= 10 && Params->out1 < 20) {
        m1 = Params->mSL[Params->out1 - 10];
        m2 = Params->mSL[Params->out2 - 10];
    }

    *masses = m1;
    *(masses + 1) = m2;

    return masses;
}

static inline double kln(double x, double y, double z)
{

    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)
                - 2.*x * y - 2.*y * z - 2.*z * x);
}

double dPS2(double &xa, double &xb, double &t, double *x, Parameters *Params)
{
    const double sh = Params->sh;
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


    const double m1s = pow2(m1);
    const double m2s = pow2(m2);

    // Conversion factor from GeV^-2 to pb.
    double djac = 389379304.0;

    // Integration variables xa and xb.
    double xamin = pow2(m1 + m2) / sh;
    double xamax = 1.0;
    double xbmin = xamin;
    double xbmax = 1.0;

    if (xamin == 0.0) {
        xa = (xamax - xamin) * x[0] + xamin;
        djac *= xamax - xamin;
        xbmin /= xa;
        xb = (xbmax - xbmin) * x[1] + xbmin;
        djac *= xbmax - xbmin;
    } else {
        xa = xamin * pow(xamax / xamin, x[0]);
        djac *= xa * log(xamax / xamin);
        xbmin /= xa;
        xb = xbmin * pow(xbmax / xbmin, x[1]);
        djac *= xb * log(xbmax / xbmin);
    }

    const double s = xa * xb * sh;

    // Integration variable t.
    const double tmin = -.5 * (s - m1s - m2s + kln(s, m1s, m2s));
    const double tmax = tmin + kln(s, m1s, m2s);
    t = (tmax - tmin) * x[2] + tmin;
    djac *= tmax - tmin;

    // Phase Space factor.
    djac *= 0.125 * M_1_PI / s;

    return djac;
}

double dPS2c(double &djacdelta, double &djacplus,
             double &xa, double &xb, double &xc,
             double &t, double &tc, double *x, Parameters *Params)
{
    const double sh = Params->sh;

    double m1;
    double m2;

    if (Params->out1 < 10) {
        m1 = Params->mCH[Params->out1];
        m2 = Params->mCH[Params->out2];
    } else  if (Params->out1 >= 10 && Params->out1 < 20) {
        m1 = Params->mSL[Params->out1 - 10];
        m2 = Params->mSL[Params->out2 - 10];
    }

    const double m1s = pow2(m1);
    const double m2s = pow2(m2);

    // Conversion factor from GeV^-2 to pb.
    double djac = 389379304.0;

    // Integration variablesxa and xb.
    double xamin = pow2(m1 + m2) / sh;
    double xamax = 1.;
    double xbmin = xamin;
    double xbmax = 1.;
    double xcmin = xbmin;
    double xcmax = 1.;


    if (xamin == 0.0) {
        xa = (xamax - xamin) * x[0] + xamin;
        djac *= xamax - xamin;
        xbmin /= xa;
        xb = (xbmax - xbmin) * x[1] + xbmin;
        djac *= xbmax - xbmin;
        djacdelta = djac;
        xcmin /= xa * xb;
        xc = (xcmax - xcmin) * x[2] + xcmin;
        djac *= xcmax - xcmin;
    } else {
        xa = xamin * pow(xamax / xamin, x[0]);
        djac *= xa * log(xamax / xamin);
        xbmin /= xa;
        xb = xbmin * pow(xbmax / xbmin, x[1]);
        djac *= xb * log(xbmax / xbmin);
        djacdelta = djac;
        xcmin /= xa * xb;
        xc = xcmin * pow(xcmax / xcmin, x[2]);
        djac *= xc * log(xcmax / xcmin);
    }

    djacplus = djac;
    const double s = xa * xb * sh;

    // Integration variable t.
    const double tmin = -0.5 * (s - m1s - m2s + kln(s, m1s, m2s));
    const double tmax = tmin + kln(s, m1s, m2s);
    t = (tmax - tmin) * x[3] + tmin;
    djacdelta *= tmax - tmin;
    djacplus *= tmax - tmin;

    const double sc = xc * s;

    // Integration variable tc.
    const double tcmin = -0.5 * (sc - m1s - m2s + kln(sc, m1s, m2s));
    const double tcmax = tcmin + kln(sc, m1s, m2s);
    tc = (tcmax - tcmin) * x[3] + tcmin;
    djac *= tcmax - tcmin;

    // Phase Space factor.
    djacdelta *= 0.125 * M_1_PI / s;
    djacplus *= 0.125 * M_1_PI / s;
    djac     *= 0.125 * M_1_PI / sc;

    return djac;
}

double dPS3(double &xa, double &xb, double &M2, double &pt2,
            double &Tp, double &Pp, double *x, Parameters *Params)
{
    const double sh = Params->sh;

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

    const double m1s = pow2(m1);
    const double m2s = pow2(m2);

    // Conversion factor from GeV^-2 to pb.
    double djac = 389379304.0;

    // Integration variable xa, xb and M2
    double xamin = pow2(m1 + m2) / sh;
    double xamax = 1.;
    double xbmin = xamin;
    double xbmax = 1.;
    double M2min = xamin * sh;
    double M2max = sh;

    if (xamin == 0.0) {
        xa = (xamax - xamin) * x[0] + xamin;
        djac *= xamax - xamin;
        xbmin /= xa;
        xb = (xbmax - xbmin) * x[1] + xbmin;
        djac *= xbmax - xbmin;
        M2max *= xa * xb;
        M2 = (M2max - M2min) * x[2] + M2min;
        djac *= M2max - M2min;
    } else {
        xa = xamin * pow(xamax / xamin, x[0]);
        djac *= xa * log(xamax / xamin);
        xbmin /= xa;
        xb = xbmin * pow(xbmax / xbmin, x[1]);
        djac *= xb * log(xbmax / xbmin);
        M2max *= xa * xb;
        M2 = M2min * pow(M2max / M2min, x[2]);
        djac *= M2 * log(M2max / M2min);
    }

    const double s = xa * xb * sh;

    double pt2min = 0.0;
    double pt2max = 0.25 * pow2(s - M2) / s;

    pt2 = (pt2max - pt2min) * x[3] + pt2min;

    djac *= pt2max - pt2min;

    // Integration variable Tp.
    const double Tpmin = 0.0;
    const double Tpmax = M_PI;
    Tp = (Tpmax - Tpmin) * x[4] + Tpmin;
    djac *= Tpmax - Tpmin;

    // Integration variable Pp.
    const double Ppmin = 0.0;
    const double Ppmax = M_PI;
    djac *= 2.0;
    Pp = (Ppmax - Ppmin) * x[5] + Ppmin;
    djac *= Ppmax - Ppmin;

    // Phase Space factor
    djac *= kln(M2, m1s, m2s) * sin(Tp)
            / sqrt(pow2(s - M2) - 4.0 * s * pt2) / M2 / pow4(4.0 * M_PI);

    return djac;
}

double IB(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj;
    double xa, xb, t;
    const double djac = dPS2(xa, xb, t, x, jp);
    const double s = xa * xb * jp->sh;

    double dij = 1.0;
    if (jp->out1 < 10 && jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
        dij = 0.5;
    }

    // PDFs of gluons and quarks for incoming particle a and b.
    double ga, qa[2][6];
    double gb, qb[2][6];

    pdfX(ga, qa, xa, jp->mufs);
    pdfX(gb, qb, xb, jp->mufs);

    double sig = 0.0;

    // Sums over all possible initial states.
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            // Set initial state
            jp->in1 = i0;
            jp->in2 = i1;

            // Test charge conservation.
            // down type: 0 1 2; up type: 3 4 5
            // neutralino: 0 1 2 3; charginos 4 5
            // sneutrinos: 10 11 12; sleptons: 13 14 15 16 17 18
            // 12 and 18 works, but 12 and 19 not

            if ((jp->out1 < 10 && jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
                    || (jp->out1 >= 10 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 10) / 3 - (jp->out2 - 10) / 3)
                    || (jp->out1 >= 16 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 10) / 3)
                    || (jp->out2 >= 16 && jp->out2 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 13) / 3)) {
                sig += (qa[0][i0] * qb[1 - jp->ic][i1] + qb[jp->ic][i0] * qa[1][i1])
                       * Born(s, t, jp);
            }
        }
    }
    // Flux, symmetry factor, spin and color average.
    return dij * sig * djac / (24.0 * s);
}

double IV(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj; // conversion of type void* into IOS*
    double xa, xb, t;
    const double djac = dPS2(xa, xb, t, x, jp);
    const double s = xa * xb * jp->sh;

    double dij = 1.0;
    if (jp->out1 < 10 && jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
        dij = 0.5;
    }

    // PDFs of the gluons and the quarks for incoming particles a and b.
    double ga, qa[2][6];
    pdfX(ga, qa, xa, jp->mufs);
    double gb, qb[2][6];
    pdfX(gb, qb, xb, jp->mufs);

    double sig = 0.;
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            jp->in1 = i0;
            jp->in2 = i1;

            // Charge conservation.
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
    // Flux, symmetry factor, spin and color average.
    return aS(jp->murs, jp->set) * dij * sig * djac * 4.0 / (24.0 * s * 3.0);
}


// Collinear emission
double IC(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj;
    double djacdelta, djacplus, xa, xb, xc, t, tc;
    const double djac = dPS2c(djacdelta, djacplus, xa, xb, xc, t, tc, x, jp);
    const double s = xa * xb * jp->sh;
    const double sc = xc * s;
    double z;

    if (jp->out1 >= 10 && jp->out1 < 20) {
        z = pow2(jp->mSL[jp->out1 - 10] + jp->mSL[jp->out2 - 10]) / s;
    } else {
        z = pow2(jp->mCH[jp->out1] + jp->mCH[jp->out2]) / s;
    }

    double dij = 1.0;
    if (jp->out1 < 10 && jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
        dij = 0.5;
    }

    double ga, qa[2][6];
    pdfX(ga, qa, xa, jp->mufs);
    double gb, qb[2][6];
    pdfX(gb, qb, xb, jp->mufs);

    double sig = 0.;

    // Sums over all possible initial states.
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            jp->in1 = i0;
            jp->in2 = i1;

            if ((jp->out1 < 10 && jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
                    || (jp->out1 >= 10 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 10) / 3 - (jp->out2 - 10) / 3)
                    || (jp->out1 >= 16 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 10) / 3)
                    || (jp->out2 >= 16 && jp->out2 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 13) / 3)) {

                // Splitting functions Pqq and Kqq.
                sig += (qa[0][i0] * qb[1 - jp->ic][i1] + qb[jp->ic][i0] * qa[1][i1])
                       * ((-(1. + pow2(xc)) / (1. - xc) * log(jp->mufs / s)
                           + (2. / (1. - xc) - 1. - xc) * 2.*log(1. - xc) + 1. - xc)
                          * Born(sc, tc, jp) * djac / xc
                          + (((1. + pow2(xc)) / (1. - xc) * log(jp->mufs / s)
                              - 4. / (1. - xc) * log(1. - xc))
                             * djacplus
                             + (-pow2(M_PI) / 6.0 - (.5 * z * z + z + 2.*log(1. - z)) * log(jp->mufs / s)
                                + 2.*pow2(log(1. - z))) * djacdelta)
                          * Born(s, t, jp)) * 8.0 / 3.0;

                // Splitting functions Pgq, Kgq, Pgqb and Kgqb.
                sig += (qa[0][i0] * gb + qb[jp->ic][i0] * ga
                        + qa[1][i1] * gb + qb[1 - jp->ic][i1] * ga)
                       * (-(pow(xc, 2) + pow(1. - xc, 2)) *
                          log(jp->mufs / s / pow(1. - xc, 2))
                          + 2.*xc * (1. - xc))
                       * Born(sc, tc, jp) * djac * .5 / xc;
            }
        }
    }

    // Flux, symmetry factor, spin and color average.
    return dij / 24. / s * sig * aS(jp->murs, jp->set);
}

double IG(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj;
    double xa, xb, mi2, pt2, th, ph;
    const double djac = dPS3(xa, xb, mi2, pt2, th, ph, x, jp);
    const double s = xa * xb * jp->sh;

    double dij = 1.0;
    if (jp->out1 < 10 && jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
        dij = 0.5;
    }

    double ga, qa[2][6];
    pdfX(ga, qa, xa, jp->mufs);
    double gb, qb[2][6];
    pdfX(gb, qb, xb, jp->mufs);

    double sig = 0.;
    // Sums over all possible initial states.
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            jp->in1 = i0;
            jp->in2 = i1;

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
    // Flux, symmetry factor, spin and color average.
    return dij * pow(M_PI, 2) / 4.5 / s * sig * djac * aS(jp->murs, jp->set);

}

double IQ(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj; // conversion of type void* into IOS*
    double xa, xb, mi2, pt2, th, ph;
    const double djac = dPS3(xa, xb, mi2, pt2, th, ph, x, jp);
    const double s = xa * xb * jp->sh;

    double dij = 1.0;
    if (jp->out1 < 10 && jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
        dij = 0.5;
    }

    double ga, qa[2][6];
    pdfX(ga, qa, xa, jp->mufs);
    double gb, qb[2][6];
    pdfX(gb, qb, xb, jp->mufs);

    double sig = 0.;
    // Sums over all possible initial states.
    for (int i0 = 0; i0 < 5; i0++) {
        for (int i1 = 0; i1 < 5; i1++) {
            jp->in1 = i0;
            jp->in2 = i1;

            // Tests charge conservation.
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

    // Flux, symmetry factor, spin and color average.
    return dij * pow(M_PI, 2) / 12. / s * sig * djac * aS(jp->murs, jp->set);
}

double IR(double *x, size_t dim, void *prm)
{
    Parameters *jp = (Parameters *)prm;

    const double sh = jp->sh;

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


    const double m1s = pow2(m1);
    const double m2s = pow2(m2);

    // Conversion factor from GeV^-2 to pb.
    complex<double> djac(389379304.0, 0.);

    // M2 integration.
    double M2min = pow2(m1 + m2);
    const double M2max = sh;

    double s = M2min * pow(M2max / M2min, x[2]);
    djac *= log(M2max / M2min);

    // Inverse Mellin.
    const complex<double> ephi(-M_SQRT1_2, M_SQRT1_2);
    const complex<double> nm
        = (M2min / sh / 0.09 + 0.09) - jp->a1min - ephi * log(x[0]);
    djac *= ephi / x[0];

    // t integration.
    const double tmin = -0.5 * (s - m1s - m2s + kln(s, m1s, m2s));
    const double tmax = tmin + kln(s, m1s, m2s);
    double t = (tmax - tmin) * x[1] + tmin;
    djac *= tmax - tmin;
    djac /= 8.*M_PI * s;
    return M_1_PI *
           imag(djac * pow(s / sh, -nm + 1.) * Thadronic_xs2(nm, s, t, jp));
}

void hadronic_xs(double &res, double &err, double &chi2, int Flag, int Verb, Parameters *Params)
{

    // Selects the integrand.
    size_t calls = 0;
    gsl_monte_function I;
    switch (Flag) {
    case 0:
        I.f = &IB;
        I.dim = 3;
        I.params = Params;
        calls = 15000;
        break;
    case 1:
        I.f = &IV;
        I.dim = 3;
        I.params = Params;
        calls = 1500;
        break;
    case 2:
        I.f = &IC;
        I.dim = 4;
        I.params = Params;
        calls = 20000;
        break;
    case 3:
        I.f = &IG;
        I.dim = 6;
        I.params = Params;
        calls = 30000;
        break;
    case 4:
        I.f = &IQ;
        I.dim = 6;
        I.params = Params;
        calls = 30000;
        break;
    case 5:
        I.f = &IR;
        I.dim = 3;
        I.params = Params;
        calls = 30000;
        break;
    default:
        cout << "hadronic_xs: Flag=" << Flag << endl;
        exit(0);
    }
    const size_t dnum = I.dim;

    double xmin[dnum], xmax[dnum];
    for (size_t i0 = 0; i0 < dnum; i0++) {
        xmin[i0] = 0.0;
        xmax[i0] = 1.0;
    }

    // Initializes the integration routine.
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
