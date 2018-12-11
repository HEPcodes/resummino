// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Implements the kinematics module.

#include <cmath>
#include <complex>
#include <iostream>
#include <cstdlib>
#include "fi.h"
#include <stdio.h>
#include "utils.h"

using namespace std;

#define WIDTH 1.e-2

static inline double kln(double x, double y, double z)
{
    return sqrt(pow2(x) + pow2(y) + pow2(z)
                - 2.0 * x * y - 2.0 * y * z - 2.0 * z * x);
}

void FI::SetKinematic(const double mi, const double mj,
                      const double s,  const double t)
{
    m1 = mi; // Mass of outgoing antiparticle
    m2 = mj; // Mass of outgoing particle
    m1s = pow2(m1);
    m2s = pow2(m2);
    papb = 0.5 * s;
    p1p2 = 0.5 * (s - m1s - m2s);
    const double u = m1s + m2s - s - t;
    pap1 = 0.5 * (m1s - t);
    pbp2 = 0.5 * (m2s - t);
    pbp1 = 0.5 * (m1s - u);
    pap2 = 0.5 * (m2s - u);
}

void FI::SetKinematic(const double mi,  const double mj,  const double s,
                      const double invariant_mass2, const double pt2,
                      const double th, const double ph, const int ys)
{
    m1 = mi;
    m1s = pow2(m1);
    m2 = mj;
    m2s = pow2(m2);
    papb = 0.5 * s;
    p1p2 = 0.5 * (invariant_mass2 - m1s - m2s);

    if (ys == 0) {
        pap3 = 0.25 * (s - invariant_mass2 + sqrt(pow2(s - invariant_mass2) - 4.0 * s * pt2));
        pbp3 = 0.25 * (s - invariant_mass2 - sqrt(pow2(s - invariant_mass2) - 4.0 * s * pt2));
    } else if (ys == 1) {
        pap3 = 0.25 * (s - invariant_mass2 - sqrt(pow(s - invariant_mass2, 2) - 4.0 * s * pt2));
        pbp3 = 0.25 * (s - invariant_mass2 + sqrt(pow(s - invariant_mass2, 2) - 4.0 * s * pt2));
    } else {
        cout << "Rapidity flag = " << ys << "\n";
        exit(0);
    }

    double kM12, cA, sA;
    kM12 = kln(invariant_mass2, m1s, m2s);
    cA = 1. - 2.0 * s * invariant_mass2 / (invariant_mass2 + 2.0 * pap3) / (invariant_mass2 + 2.0 * pbp3);
    sA = sqrt(1. - cA * cA);

    pap1 = .25 * (invariant_mass2 + 2.0 * pbp3) / invariant_mass2 * (invariant_mass2 + m1s - m2s - kM12 * cos(th));
    pbp1 = .25 * (invariant_mass2 + 2.0 * pap3) / invariant_mass2 *
           (invariant_mass2 + m1s - m2s - kM12 * (cA * cos(th) + sA * sin(th) * cos(ph)));
    pap2 = papb - pap1 - pap3;
    pbp2 = papb - pbp1 - pbp3;
    p1p3 = pap1 + pbp1 - m1s - p1p2;
    p2p3 = pap3 + pbp3 - p1p3;
}

void FI::SetDipKinematicA(double &x, double &fact)
{
    x = 1. - (pap3 + pbp3) / papb;
    fact = 1. / x / pap3;

    papb *= x;
    const double L = 4.0 * papb + (1. - x) * pap3;
    pap1 = (pap3 * m1s + pap3 * p1p2 - 2.0 * papb * p1p3 + 2.0 * papb * pap1 + x * pap3 * m1s
            + x * pap3 * p1p3 + x * pap3 * p1p2 + 2.0 * x * papb * pap1) / L;
    pap2 = (pap3 * m2s + pap3 * p1p2 - 2.0 * papb * p2p3 + 2.0 * papb * pap2 + x * pap3 * m2s
            + x * pap3 * p2p3 + x * pap3 * p1p2 + 2.0 * x * papb * pap2) / L;
    pbp1 = (2.0 * papb * m1s + 2.0 * papb * p1p2 + 2.0 * papb * pbp1 - 2.0 * x * pap3 * m1s
            - x * pap3 * p1p3 - 2.0 * x * pap3 * p1p2 - 2.0 * x * papb * m1s - 2.0 * x * papb * p1p3
            - 2.0 * x * papb * p1p2 + 2.0 * x * papb * pbp1) / L;
    pbp2 = (2.0 * papb * m2s + 2.0 * papb * p1p2 + 2.0 * papb * pbp2 - 2.0 * x * pap3 * m2s
            - x * pap3 * p2p3 - 2.0 * x * pap3 * p1p2 - 2.0 * x * papb * m2s - 2.0 * x * papb * p2p3
            - 2.0 * x * papb * p1p2 + 2.0 * x * papb * pbp2) / L;
}

void FI::SetDipKinematicB(double & x, double & fact)
{
    x = 1. - (pap3 + pbp3) / papb;
    fact = 1. / x / pbp3;

    papb *= x;
    const double L = 4.0 * papb + (1. - x) * pbp3;
    pbp1 = (pbp3 * m1s + pbp3 * p1p2 - 2.0 * papb * p1p3 + 2.0 * papb * pbp1 + x * pbp3 * m1s
            + x * pbp3 * p1p3 + x * pbp3 * p1p2 + 2.0 * x * papb * pbp1) / L;
    pbp2 = (pbp3 * m2s + pbp3 * p1p2 - 2.0 * papb * p2p3 + 2.0 * papb * pbp2 + x * pbp3 * m2s
            + x * pbp3 * p2p3 + x * pbp3 * p1p2 + 2.0 * x * papb * pbp2) / L;
    pap1 = (2.0 * papb * m1s + 2.0 * papb * p1p2 + 2.0 * papb * pap1 - 2.0 * x * pbp3 * m1s
            - x * pbp3 * p1p3 - 2.0 * x * pbp3 * p1p2 - 2.0 * x * papb * m1s - 2.0 * x * papb * p1p3
            - 2.0 * x * papb * p1p2 + 2.0 * x * papb * pap1) / L;
    pap2 = (2.0 * papb * m2s + 2.0 * papb * p1p2 + 2.0 * papb * pap2 - 2.0 * x * pbp3 * m2s
            - x * pbp3 * p2p3 - 2.0 * x * pbp3 * p1p2 - 2.0 * x * papb * m2s - 2.0 * x * papb * p2p3
            - 2.0 * x * papb * p1p2 + 2.0 * x * papb * pap2) / L;
}


void FI::SetPropagator(const double mass1, const double mass2, const double width1, const double width2)
{

    ivs = 0.5 / papb;
    complex<double> II(0.0, WIDTH);
    complex<double> III(0.0, 1.0);

    // Electroweak s-channel propagators in t-channel diagrams with real quark
    // emissions.
    ivs1s1 = 1. / (m2s + 2.*p2p3 - std::pow(mass1, 2) -  mass1 * width1 * III);
    ivs2s1 = 1. / (m1s + 2.*p1p3 - std::pow(mass1, 2) -  mass1 * width1 * III);
    ivs1s2 = 1. / (m2s + 2.*p2p3 - std::pow(mass2, 2) +  mass2 * width2 * III);
    ivs2s2 = 1. / (m1s + 2.*p1p3 - std::pow(mass2, 2) +  mass2 * width2 * III);

    // Electroweak s-channel propagators.
    ivs3v1 = 1. / (m1s + 2.*p1p2 + m2s - std::pow(mass1, 2) - mass1 * III * width1);
    ivs3v2 = 1. / (m1s + 2.*p1p2 + m2s - std::pow(mass2, 2) + mass2 * III * width2);

    // Propagators without width.
    ivt1s1 = 1. / (m1s - 2.*pap1 - std::pow(mass1, 2));
    ivt1s2 = 1. / (m1s - 2.*pap1 - std::pow(mass2, 2));
    ivt2s1 = 1. / (m2s - 2.*pbp2 - std::pow(mass1, 2));
    ivt2s2 = 1. / (m2s - 2.*pbp2 - std::pow(mass2, 2));

    // Strong quark propagators for real corrections.
    ivt3 = -0.5 / pap3;
    ivu3 = -0.5 / pbp3;

    // Electroweak u-channel propagators.
    // Propagators without width.
    ivu1s1 = 1. / (m1s - 2.*pbp1 - std::pow(mass1, 2));
    ivu1s2 = 1. / (m1s - 2.*pbp1 - std::pow(mass2, 2));
    ivu2s1 = 1. / (m2s - 2.*pap2 - std::pow(mass1, 2));
    ivu2s2 = 1. / (m2s - 2.*pap2 - std::pow(mass2, 2));

}

void FI::SetLoopMass(const double mls[5])
{
    mul = mls[0]; // Renormalization scale
    m1l = pow(mls[1], 2); // Squared masses in loops
    m2l = pow(mls[2], 2);
    m3l = pow(mls[3], 2);
    m4l = pow(mls[4], 2);
    ml1 = sqrt(m1l); // Masses in loops
    ml2 = sqrt(m2l);
    ml3 = sqrt(m3l);
    ml4 = sqrt(m4l);
}

void FI::SetSCoupling(struct Coupling C[2])
{
    LL = C[0].L * C[1].L;
    LR = C[0].L * C[1].R;
    RR = C[0].R * C[1].R;
    RL = C[0].R * C[1].L;
}

void FI::SetWCoupling(struct Coupling C[4])
{
    LLLL = C[0].L * C[1].L * conj(C[2].L * C[3].L);
    LLLR = C[0].L * C[1].L * conj(C[2].L * C[3].R);
    LLRL = C[0].L * C[1].L * conj(C[2].R * C[3].L);
    LRLL = C[0].L * C[1].R * conj(C[2].L * C[3].L);
    RLLL = C[0].R * C[1].L * conj(C[2].L * C[3].L);
    LLRR = C[0].L * C[1].L * conj(C[2].R * C[3].R);
    LRLR = C[0].L * C[1].R * conj(C[2].L * C[3].R);
    LRRL = C[0].L * C[1].R * conj(C[2].R * C[3].L);

    RRRR = C[0].R * C[1].R * conj(C[2].R * C[3].R);
    RRRL = C[0].R * C[1].R * conj(C[2].R * C[3].L);
    RRLR = C[0].R * C[1].R * conj(C[2].L * C[3].R);
    RLRR = C[0].R * C[1].L * conj(C[2].R * C[3].R);
    LRRR = C[0].L * C[1].R * conj(C[2].R * C[3].R);
    RRLL = C[0].R * C[1].R * conj(C[2].L * C[3].L);
    RLRL = C[0].R * C[1].L * conj(C[2].R * C[3].L);
    RLLR = C[0].R * C[1].L * conj(C[2].L * C[3].R);
}
