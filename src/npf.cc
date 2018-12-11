// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Implements the n-point functions.

#include <cmath>
#include <complex>
#include <iostream>
#include <cstdlib>
#include "npf.h"

void Npf::SetIeps(int Index)
{
    ieps = Index;
}

void Npf::GetNpf(std::complex<double> Array[7], int Index)
{
    // Initialization
    for (size_t i0 = 0; i0 < 7; i0++) {
        Array[i0] = std::complex<double>(0., 0.);
    }

    // Npoint
    const int ipt = npt + Index;
    switch (ipt) {
    case 10:
        GetA0(Array[0]);
        break;
    case 20:
        GetB0(Array[0]);
        break;
    case 21:
        GetB1(Array);
        break;
    case 22:
        GetB2(Array);
        break;
    case 30:
        GetC0(Array[0]);
        break;
    case 31:
        GetC1(Array);
        break;
    case 32:
        GetC2(Array);
        break;
    case 40:
        GetD0(Array[0]);
        break;
    case 41:
        GetD1(Array);
        break;
    case 42:
        GetD2(Array);
        break;
    default:
        std::cout << "npt=" << ipt << "\n";
        exit(0);
    }
}

void Npf::SetArgument(const double m1)
{
    m1a = m1;
}

void Npf::SetArgument(const double p1, const double m1, const double m2)
{
    p1b = p1;
    m1b = m1;
    m2b = m2;
    f1b = m2b - m1b - p1b;
}

void Npf::SetArgument(const double p1, const double p2, const double p3,
                      const double m1, const double m2, const double m3)
{
    p1c = p1;
    p2c = p2;
    p3c = p3;
    m1c = m1;
    m2c = m2;
    m3c = m3;
    f1c = m2c - m1c - p1c;
    f2c = m3c - m2c - p3c + p1c;
    double detXc = p1c * p2c - .25 * (p3c - p1c - p2c) * (p3c - p1c - p2c);
    Xc[0][0] = p2c / detXc;
    Xc[0][1] = -.5 * (p3c - p1c - p2c) / detXc;
    Xc[1][0] = Xc[0][1];
    Xc[1][1] = p1c / detXc;
}

void Npf::SetArgument(const double p1, const double p2, const double p3,
                      const double p4, const double s12, const double s23,
                      const double m1, const double m2, const double m3,
                      const double m4)
{
    p1d = p1;
    p2d = p2;
    p3d = p3;
    p4d = p4;
    s12d = s12;
    s23d = s23;
    m1d = m1;
    m2d = m2;
    m3d = m3;
    m4d = m4;
    f1d = m2d - m1d - p1d;
    f2d = m3d - m2d - s12d + p1d;
    f3d = m4d - m3d - p4d + s12d;
    double detXd = (4.*p1d * p2d * p3d - p3d * std::pow(p1d + p2d - s12d, 2) -
                    p1d * std::pow(p2d + p3d - s23d, 2) + (p1d + p2d - s12d) *
                    (p2d + p3d - s23d) * (p2d + p4d - s12d - s23d) -
                    p2d * std::pow(p2d + p4d - s12d - s23d, 2)) * .25;

    Xd[0][0] = (p2d * p3d - std::pow(p2d + p3d - s23d, 2) * .25) / detXd;
    Xd[0][1] = (2.*p3d * (p1d + p2d - s12d) - (p2d + p3d - s23d)
                * (p2d + p4d - s12d - s23d)) * .25 / detXd;
    Xd[0][2] = (-std::pow(p2d, 2) + p1d * (p2d + p3d - s23d) + s12d * (-p3d + s23d) +
                p2d * (p3d - 2.*p4d + s12d + s23d)) * .25 / detXd;
    Xd[1][1] = (p1d * p3d - std::pow(p2d + p4d - s12d - s23d, 2) * .25) / detXd;
    Xd[1][2] = (2.*p1d * (p2d + p3d - s23d) - (p1d + p2d - s12d)
                * (p2d + p4d - s12d - s23d)) * .25 / detXd;
    Xd[2][2] = (p1d * p2d - std::pow(p1d + p2d - s12d, 2) * .25) / detXd;
    Xd[1][0] = Xd[0][1];
    Xd[2][0] = Xd[0][2];
    Xd[2][1] = Xd[1][2];
}

void Npf::GetA0(std::complex<double> &a0)
{
    // A0
    if (ieps < 0 || ieps > 2) {
        a0 = std::complex<double>(0., 0.);
    } else {
        int qleps = -ieps;
        a0 = qli1_(m1a, mu, qleps);
    }
}

void Npf::GetB0(std::complex<double> &b0)
{
    // B0
    if (ieps < 0 || ieps > 2) {
        b0 = std::complex<double>(0., 0.);
    } else {
        int qleps = -ieps;
        b0 = qli2_(p1b, m1b, m2b, mu, qleps);
    }
}

void Npf::GetB1(std::complex<double> *b1)
{
    // B1
    std::complex<double> a0, b0;
    if (p1b != 0.) {
        GetB0(b0);
        b1[0] = f1b * b0;
        if (m1b != m2b) {
            SetArgument(m1b);
            GetA0(a0);
            b1[0] += a0;
            SetArgument(m2b);
            GetA0(a0);
            b1[0] -= a0;
        }
        b1[0] *= .5 / p1b;
    } else if (f1b != 0.) {
        ieps += 2;
        SetArgument(m2b);
        GetA0(a0);
        GetB0(b0);
        b1[0] = m1b * b0 + a0;
        ieps -= 1;
        SetArgument(m2b);
        GetA0(a0);
        GetB0(b0);
        b1[0] = .5 * b1[0] + m1b * b0 + a0;
        ieps -= 1;
        SetArgument(m2b);
        GetA0(a0);
        GetB0(b0);
        b1[0] = .5 * b1[0] + m1b * b0 + a0;
        b1[0] = .5 * b1[0] - a0;
        b1[0] /= f1b;
    } else {
        GetB0(b0);
        b1[0] = -.5 * b0;
    }
}

void Npf::GetB2(std::complex<double> *b2)
{
    // B22 ---------------------------------------------------------------------
    std::complex<double> a0, b0, b1[1];
    ieps += 2;
    SetArgument(m2b);
    GetA0(a0);
    GetB0(b0);
    GetB1(b1);
    b2[1] = 2.*m1b * b0 - f1b * b1[0] + a0;
    ieps -= 1;
    SetArgument(m2b);
    GetA0(a0);
    GetB0(b0);
    GetB1(b1);
    b2[1] = 2. / 3.*b2[1] + 2.*m1b * b0 - f1b * b1[0] + a0;
    ieps -= 1;
    SetArgument(m2b);
    GetA0(a0);
    GetB0(b0);
    GetB1(b1);
    b2[1] = 2. / 3.*b2[1] + 2.*m1b * b0 - f1b * b1[0] + a0;
    b2[1] /= 6.;
    // B21 ---------------------------------------------------------------------
    if (p1b != 0.) {
        b2[0] = (f1b * b1[0] + a0 - 2.*b2[1]) * .5 / p1b;
    } else if (f1b != 0.) {
        std::cout << "B21(0,m1,m2) not yet implemented!\n";
    } else {
        GetB0(b0);
        b2[0] = b0 / 3.;
    }
}

void Npf::GetC0(std::complex<double> &c0)
{
    // C0
    if (ieps < 0 || ieps > 2) {
        c0 = std::complex<double>(0., 0.);
    } else {
        int qleps = -ieps;
        c0 = qli3_(p1c, p2c, p3c, m1c, m2c, m3c, mu, qleps);
    }
}

void Npf::GetC1(std::complex<double> *c1)
{
    // C11, C12 ----------------------------------------------------------------
    std::complex<double> b0, c0;
    GetC0(c0);
    std::complex<double> rr[2] = { f1c * c0, f2c * c0 };
    SetArgument(p1c, m1c, m2c);
    GetB0(b0);
    rr[1] += b0;
    SetArgument(p2c, m2c, m3c);
    GetB0(b0);
    rr[0] -= b0;
    SetArgument(p3c, m1c, m3c);
    GetB0(b0);
    rr[0] += b0;
    rr[1] -= b0;
    c1[0] = .5 * (Xc[0][0] * rr[0] + Xc[0][1] * rr[1]);
    c1[1] = .5 * (Xc[1][0] * rr[0] + Xc[1][1] * rr[1]);
}

void Npf::GetC2(std::complex<double> *c2)
{
    // C24 ---------------------------------------------------------------------
    std::complex<double> b0, b1[1], c0, c1[2];
    ieps += 2;
    SetArgument(p2c, m2c, m3c);
    GetB0(b0);
    GetC0(c0);
    GetC1(c1);
    c2[3] = 2.*m1c * c0 - f1c * c1[0] - f2c * c1[1] + b0;
    ieps -= 1;
    SetArgument(p2c, m2c, m3c);
    GetB0(b0);
    GetC0(c0);
    GetC1(c1);
    c2[3] += 2.*m1c * c0 - f1c * c1[0] - f2c * c1[1] + b0;
    ieps -= 1;
    SetArgument(p2c, m2c, m3c);
    GetB0(b0);
    GetC0(c0);
    GetC1(c1);
    c2[3] += 2.*m1c * c0 - f1c * c1[0] - f2c * c1[1] + b0;
    c2[3] *= .5;
    // C21, C22, C23 -----------------------------------------------------------
    std::complex<double> rr[4] = { f1c*c1[0], f2c*c1[0], f1c*c1[1], f2c*c1[1] };
    SetArgument(p1c, m1c, m2c);
    GetB1(b1);
    rr[1] += b1[0];
    SetArgument(p2c, m2c, m3c);
    GetB0(b0);
    GetB1(b1);
    rr[0] += b0;
    rr[2] -= b1[0];
    SetArgument(p3c, m1c, m3c);
    GetB1(b1);
    rr[0] += b1[0];
    rr[1] -= b1[0];
    rr[2] += b1[0];
    rr[3] -= b1[0];
    c2[0] = Xc[0][0] * (rr[0] - c2[3]) + Xc[0][1] * rr[1];
    c2[2] = Xc[1][0] * (rr[0] - c2[3]) + Xc[1][1] * rr[1];
//   c2[2]=Xc[0][1]*(rr[3]-c2[3])+Xc[0][0]*rr[2];
    c2[1] = Xc[1][1] * (rr[3] - c2[3]) + Xc[1][0] * rr[2];
    for (size_t i0 = 0; i0 < 4; i0++) {
        c2[i0] *= .5;
    }
}

void Npf::GetD0(std::complex<double> &d0)
{
    // D0
    if (ieps < 0 || ieps > 2) {
        d0 = std::complex<double>(0., 0.);
    } else {
        int qleps = -ieps;
        d0 = qli4_(p1d, p2d, p3d, p4d, s12d, s23d, m1d, m2d, m3d, m4d, mu, qleps);
    }
}

void Npf::GetD1(std::complex<double> *d1)
{
    // D11, D12, D13
    std::complex<double> c0, d0;
    GetD0(d0);
    std::complex<double> rr[3] = { f1d * d0, f2d * d0, f3d * d0 };
    SetArgument(p1d, p2d, s12d, m1d, m2d, m3d);
    GetC0(c0);
    rr[2] += c0;
    SetArgument(p1d, s23d, p4d, m1d, m2d, m4d);
    GetC0(c0);
    rr[1] += c0;
    rr[2] -= c0;
    SetArgument(s12d, p3d, p4d, m1d, m3d, m4d);
    GetC0(c0);
    rr[0] += c0;
    rr[1] -= c0;
    SetArgument(p2d, p3d, s23d, m2d, m3d, m4d);
    GetC0(c0);
    rr[0] -= c0;
    d1[0] = .5 * (Xd[0][0] * rr[0] + Xd[0][1] * rr[1] + Xd[0][2] * rr[2]);
    d1[1] = .5 * (Xd[1][0] * rr[0] + Xd[1][1] * rr[1] + Xd[1][2] * rr[2]);
    d1[2] = .5 * (Xd[2][0] * rr[0] + Xd[2][1] * rr[1] + Xd[2][2] * rr[2]);
}

void Npf::GetD2(std::complex<double> *d2)
{
    // D27
    std::complex<double> c0, c1[2], d0, d1[3];
    ieps += 2;
    SetArgument(p2d, p3d, s23d, m2d, m3d, m4d);
    GetC0(c0);
    GetD0(d0);
    GetD1(d1);
    d2[6] = 2.*m1d * d0 - f1d * d1[0] - f2d * d1[1] - f3d * d1[2] + c0;
    ieps -= 1;
    SetArgument(p2d, p3d, s23d, m2d, m3d, m4d);
    GetC0(c0);
    GetD0(d0);
    GetD1(d1);
    d2[6] = 2.*d2[6] + 2.*m1d * d0 - f1d * d1[0] - f2d * d1[1] - f3d * d1[2] + c0;
    ieps -= 1;
    SetArgument(p2d, p3d, s23d, m2d, m3d, m4d);
    GetC0(c0);
    GetD0(d0);
    GetD1(d1);
    d2[6] = 2.*d2[6] + 2.*m1d * d0 - f1d * d1[0] - f2d * d1[1] - f3d * d1[2] + c0;

    // D21, D22, D23, D24, D25, D26
    std::complex<double> rr[9] = { f1d*d1[0], f2d*d1[0], f3d*d1[0],
                                   f1d*d1[1], f2d*d1[1], f3d*d1[1],
                                   f1d*d1[2], f2d*d1[2], f3d*d1[2]
                                 };
    SetArgument(p1d, p2d, s12d, m1d, m2d, m3d);
    GetC1(c1);
    rr[2] += c1[0];
    rr[5] += c1[1];
    SetArgument(p1d, s23d, p4d, m1d, m2d, m4d);
    GetC1(c1);
    rr[1] += c1[0];
    rr[2] -= c1[0];
    rr[4] += c1[1];
    rr[5] -= c1[1];
    rr[7] += c1[1];
    rr[8] -= c1[1];
    SetArgument(s12d, p3d, p4d, m1d, m3d, m4d);
    GetC1(c1);
    rr[0] += c1[0];
    rr[1] -= c1[0];
    rr[3] += c1[0];
    rr[4] -= c1[0];
    rr[6] += c1[1];
    rr[7] -= c1[1];
    SetArgument(p2d, p3d, s23d, m2d, m3d, m4d);
    GetC0(c0);
    GetC1(c1);
    rr[0] += c0;
    rr[3] -= c1[0];
    rr[6] -= c1[1];
    d2[0] = Xd[0][0] * (rr[0] - d2[6]) + Xd[0][1] * rr[1] + Xd[0][2] * rr[2];
    d2[3] = Xd[1][0] * (rr[0] - d2[6]) + Xd[1][1] * rr[1] + Xd[1][2] * rr[2];
    d2[4] = Xd[2][0] * (rr[0] - d2[6]) + Xd[2][1] * rr[1] + Xd[2][2] * rr[2];
//   d2[3]=Xd[0][0]*rr[3]+Xd[0][1]*(rr[4]-d2[6])+Xd[0][2]*rr[5];
    d2[1] = Xd[1][0] * rr[3] + Xd[1][1] * (rr[4] - d2[6]) + Xd[1][2] * rr[5];
    d2[5] = Xd[2][0] * rr[3] + Xd[2][1] * (rr[4] - d2[6]) + Xd[2][2] * rr[5];
//   d2[4]=Xd[0][0]*rr[6]+Xd[0][1]*rr[7]+Xd[0][2]*(rr[8]-d2[6]);
//   d2[5]=Xd[1][0]*rr[6]+Xd[1][1]*rr[7]+Xd[1][2]*(rr[8]-d2[6]);
    d2[2] = Xd[2][0] * rr[6] + Xd[2][1] * rr[7] + Xd[2][2] * (rr[8] - d2[6]);
    for (size_t i0 = 0; i0 < 7; i0++) {
        d2[i0] *= .5;
    }
}
