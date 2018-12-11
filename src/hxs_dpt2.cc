// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Computes the hadronic cross section transverse-momentum distribution at LO,
// NLO (collinear, virtual, gluon emission and quark emission) and NLL.

#include <cmath>
#include <complex>
#include <iostream>
#include <cstdlib>
#include "gsl_all.h"
#include "prm.h"
#include "utils.h"
#include "pdf.h"
#include "pxs.h"

#define CMLLN 0.9 // C coeff. of the inverse Mellin transform
#define VBSSL 2.9 // V coeff. of the inverse Bessel transform

static inline double kln(double x, double y, double z)
{
    return std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)
                     - 2.*x * y - 2.*y * z - 2.*z * x);
}

double dPS3_dPT2(double &xa, double &xb, double &M2, double &pt2,
                 double &Tp, double &Pp, double *x, Parameters *Params)
{
    const double sh = Params->sh;
    // old: const double m1=Params->mCH[Params->out1];
    // old: const double m2=Params->mCH[Params->out2];

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
    pt2 = Params->pts;

    if (pt2 < 0) {
        std::cout << "IvMlln: PT too small\n";
        exit(0);
    } else if (pt2 > std::pow(sh - std::pow(m1 + m2, 2), 2) / sh * .25) {
        std::cout << "IvMlln: PT too large\n";
        exit(1);
    }

    // Jacobian initialization
    double djac = 389379304.;

    // Integration variable xa, xb and M2
    double M2min = std::pow(m1 + m2, 2);

    double M2max = sh;
    double xamin = (M2min + (pt2 + std::sqrt(pt2 * (pt2 + M2min))) * 2.) / sh;
    double xamax = 1.;
    double xbmin = xamin;
    double xbmax = 1.;
    if (xamin == 0.) {
        xa = (xamax - xamin) * x[0] + xamin;
        djac *= xamax - xamin;
        xbmin /= xa;
        xb = (xbmax - xbmin) * x[1] + xbmin;
        djac *= xbmax - xbmin;
        M2max *= xa * xb;
        M2max -= std::sqrt(M2max * pt2) * 2.;
        M2 = (M2max - M2min) * x[2] + M2min;
        djac *= M2max - M2min;
    } else {
        xa = xamin * std::pow(xamax / xamin, x[0]);
        djac *= xa * std::log(xamax / xamin);
        xbmin /= xa;
        xb = xbmin * std::pow(xbmax / xbmin, x[1]);
        djac *= xb * std::log(xbmax / xbmin);
        M2max *= xa * xb;
        M2max -= std::sqrt(M2max * pt2) * 2.;
        M2 = M2min * std::pow(M2max / M2min, x[2]);
        djac *= std::log(M2max / M2min);
    }

    const double s = xa * xb * sh;

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

std::complex<double> dPSN_dPT2(std::complex<double> &nm, double &s,
                               double &t, double *x, Parameters *Params)
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

    const double m1s = std::pow(m1, 2);
    const double m2s = std::pow(m2, 2);
    const double pt2 = Params->pts;

    if (pt2 < 0) {
        std::cout << "IvMlln: PT too small\n";
        exit(0);
    } else if (pt2 > std::pow(sh - std::pow(m1 + m2, 2), 2) / sh * .25) {
        std::cout << "IvMlln: PT too large\n";
        exit(1);
    }

    // Jacobian initialization
    std::complex<double> djac(389379304., 0.);

    // Inverse Mellin
    const std::complex<double> ephi(-M_SQRT1_2, M_SQRT1_2);
    nm = CMLLN - Params->a1min - ephi * std::log(x[1]);
    djac *= ephi / x[1];

    // M2 integration
    double M2min = std::pow(m1 + m2, 2);

    const double M2max = sh - 2.*std::sqrt(sh * pt2);
    s = M2min * std::pow(M2max / M2min, x[0]);
    djac *= std::log(M2max / M2min);

    // t integration
    const double tmin = -.5 * (s - m1s - m2s + kln(s, m1s, m2s));
    const double tmax = tmin + kln(s, m1s, m2s);
    t = (tmax - tmin) * x[2] + tmin;
    djac *= tmax - tmin;

    djac /= 8.*M_PI * s;

    return djac;
}

double IG_dPT2(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj; // conversion of type void* into IOS*
    double xa, xb, mi2, pt2, th, ph;
    const double djac = dPS3_dPT2(xa, xb, mi2, pt2, th, ph, x, jp);
    const double s = xa * xb * jp->sh;

    double dij = 1.;
    if (jp->out1 < 10 && jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
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
            //if(jp->in1/3-jp->in2/3==jp->out1/4-jp->out2/4) {
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

double IQ_dPT2(double *x, size_t dim, void *jj)
{
    Parameters *jp = (Parameters *)jj; // conversion of type void* into IOS*
    double xa, xb, mi2, pt2, th, ph;
    const double djac = dPS3_dPT2(xa, xb, mi2, pt2, th, ph, x, jp);
    const double s = xa * xb * jp->sh;

    double dij = 1.;
    if (jp->out1 < 10 && jp->out1 == jp->out2 && jp->out1 / 4 == 0) {
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
            if ((jp->out1 < 10 && jp->in1 / 3 - jp->in2 / 3 == jp->out1 / 4 - jp->out2 / 4)
                    || (jp->out1 >= 10 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 10) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out1 >= 16 && jp->out1 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 10) / 3)


                    || (jp->out2 >= 16 && jp->out2 < 20
                        && jp->in1 / 3 - jp->in2 / 3 == (jp->out1 - 13) / 3 - (jp->out2 - 13) / 3)) {
                sig += (qa[0][i0] * gb + qb[jp->ic][i0] * ga) *
                       (RealQ(s, mi2, pt2, th, ph, 0, jp) +
                        RealQ(s, mi2, pt2, th, ph, 1, jp)// -
                       );
                sig += (qa[1][i1] * gb + qb[1 - jp->ic][i1] * ga) *
                       (RealQB(s, mi2, pt2, th, ph, 0, jp) +
                        RealQB(s, mi2, pt2, th, ph, 1, jp)// -
                       );
            }
        }
    }

    // Flux, symmetry factor, spin and color average
    return dij * std::pow(M_PI, 2) / 12. / s * sig * djac * aS(jp->murs, jp->set);
}

double IR_dPT2(double *x, size_t dim, void *prm)
{
    Parameters *jp = (Parameters *)prm; // conversion of type void* into IOS*
    std::complex<double> nm;
    double s, t;
    const double sh = jp->sh;
    const double pt = std::sqrt(jp->pts);
    const std::complex<double> djac = dPSN_dPT2(nm, s, t, x, jp);

    // Inverse Bessel: First branch
    std::complex<double> ephi(M_SQRT1_2, M_SQRT1_2);
    std::complex<double> b = -ephi * std::log(x[3]) / pt;
    std::complex<double> bjac = ephi / x[3] / pt;
    // H function for Bessel transform
    const double v = VBSSL;
    const std::complex<double> II(0., 1.);
    std::complex<double> theta = M_PI * (x[4] * (2.*II * v - 1.) - II * v);
    std::complex<double> hbssl = -M_1_PI * std::exp(-II * b * pt * std::sin(theta));
    bjac *= M_PI * (2.*II * v - 1.);

    std::complex<double> result = PtXS(nm, b, s, t, jp) * bjac * b * hbssl;

    // Inverse Bessel: Second branch
    ephi = std::complex<double>(M_SQRT1_2, -M_SQRT1_2);
    b = -ephi * std::log(x[3]) / pt;
    bjac = ephi / x[3] / pt;
    // H function for Bessel transform
    theta = M_PI * (x[4] * (2.*II * v + 1.) - II * v);
    bjac *= -M_PI * (2.*II * v + 1.);
    hbssl = -M_1_PI * std::exp(-II * b * pt * std::sin(theta));

    result += PtXS(nm, b, s, t, jp) * bjac * b * hbssl;

    return .25 * M_1_PI * std::imag(result * std::pow(s / sh, -nm + 1.) * djac);
}

double IJ_dPT2(double *x, size_t dim, void *prm)
{
    Parameters *jp = (Parameters *)prm; // conversion of type void* into IOS*
    std::complex<double> nm;
    double s, t;
    const double sh = jp->sh;
    const double pt = std::sqrt(jp->pts);
    const std::complex<double> djac = dPSN_dPT2(nm, s, t, x, jp);

    // Inverse Bessel: First branch
    std::complex<double> ephi(M_SQRT1_2, M_SQRT1_2);
    std::complex<double> b = -ephi * std::log(x[3]) / pt;
    std::complex<double> bjac = ephi / x[3] / pt;
    // H function for Bessel transform
    const double v = VBSSL;
    const std::complex<double> II(0., 1.);
    std::complex<double> theta = M_PI * (x[4] * (2.*II * v - 1.) - II * v);
    std::complex<double> hbssl = -M_1_PI * std::exp(-II * b * pt * std::sin(theta));
    bjac *= M_PI * (2.*II * v - 1.);

    std::complex<double> result = JtXS(nm, b, s, t, jp) * bjac * b * hbssl;

    // Inverse Bessel: Second branch
    ephi = std::complex<double>(M_SQRT1_2, -M_SQRT1_2);
    b = -ephi * std::log(x[3]) / pt;
    bjac = ephi / x[3] / pt;
    // H function for Bessel transform
    theta = M_PI * (x[4] * (2.*II * v + 1.) - II * v);
    bjac *= -M_PI * (2.*II * v + 1.);
    hbssl = -M_1_PI * std::exp(-II * b * pt * std::sin(theta));

    result += JtXS(nm, b, s, t, jp) * bjac * b * hbssl;

    return .25 * M_1_PI * std::imag(result * std::pow(s / sh, -nm + 1.) * djac);
}

void hadronic_xs_dPT2(double &res, double &err, double &chi2,
                      int Flag, int Verb, Parameters *Params)
{

    size_t calls = 0;
    gsl_monte_function I;
    switch (Flag) {
    case 0:
        res = 0.;
        err = 0.;
        chi2 = 0.;
        return;
    case 1:
        res = 0.;
        err = 0.;
        chi2 = 0.;
        return;
    case 2:
        res = 0.;
        err = 0.;
        chi2 = 0.;
        return;
    case 3:
        I.f = &IG_dPT2;
        I.dim = 5;
        I.params = Params;
        calls = 50000;
        break;
    case 4:
        I.f = &IQ_dPT2;
        I.dim = 5;
        I.params = Params;
        calls = 50000;
        break;
    case 5:
        I.f = &IR_dPT2;
        I.dim = 5;
        I.params = Params;
        calls = 100000;
        break;
    case 6:
        I.f = &IJ_dPT2;
        I.dim = 5;
        I.params = Params;
        calls = 100000;
        break;
    default:
        std::cout << "hadronic_xs_dPT2: Flag=" << Flag << std::endl;
        exit(0);
    }
    const size_t dnum = I.dim;

    double xmin[dnum], xmax[dnum];
    for (size_t i0 = 0; i0 < dnum; i0++) {
        xmin[i0] = 0.;
        xmax[i0] = 1.;
    }

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

    s->iterations = 1;
    s->stage = 3;

    do {
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
