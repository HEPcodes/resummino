// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Bubble matrix elements.
//
// Notes on notation:
// - MVsbuns:
//   - M = squared Matrix element
//   - V = virtual correction
//   - s = channel 1 is s-channel
//   - bu = channel 1 includes self energy bubble n
//   - s = channel 2 is s-channel (Born level)
//
// See mbn.cc for more details on notation.

#include <complex>
#include "fi.h"
#include "npf.h"

using namespace std;

// Process: q + \bar{q} -> f + \bar{f}.

double FI::MVsbu1s(int ieps)
{

    //Squared matrix elements of O(epsilon^0).
    complex<double> me0 = complex<double> (0., 0.);
    //Squared matrix elements of O(epsilon^1).
    complex<double> me1 = complex<double> (0., 0.);
    //Squared matrix elements of O(epsilon^2).
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];

    // Constructs general two point function:
    // B(incoming momentum squared, loop mass1, loop mass2, renormalization scale,
    // term proportional to a certain order of 1/epsilon ).
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);

    // General idea: M_final = M B  (some squared matrix element times a complicated divergent integral)
    // M = M0 + eps M1 + eps^2 M2 + O(eps^3)
    // B = 1/eps^2 B2 + 1/eps B1 + B0 + O(eps)
    // -> M_final = B0 M0 + B1 M1 + B2 M2 + terms going to zero + divergent terms which we absorb
    // Don't get confused. Number refers to the order of epsilon and not to a certain B-function!

    // Two-point function prop. to O(1/epsilon^0).
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];

    // Part of scalar integral of O(1/eps^0) times squared matrix element of O(eps^0).
    me0 = bb1 * ivs3v1 * ivs3v2 * (16.*pap2 * pbp1 * RL * LLLL + 16.*pap2 * pbp1 * LR *
                                   RRRR + 16.*pap1 * pbp2 * RL * RLRL + 16.*pap1 * pbp2 * LR * LRLR + 8.*m1 * m2 *
                                   papb * RL * RLLL + 8.*m1 * m2 * papb * RL * LLRL + 8.*m1 * m2 * papb * LR * LRRR + 8
                                   * m1 * m2 * papb * LR * RRLR);

    bb->SetIeps(ieps + 1); // get two-point function prop. to O(1/epsilon^1).
    bb->GetNpf(la, 1);
    bb1 = la[0];

    // Part of scalar integral of O(1/eps^1) times squared matrix element of O(eps^1).
    me1 = bb1 * ivs3v1 * ivs3v2 * (24.*pap2 * pbp1 * RL * RLRL - 24.*pap2 * pbp1 * RL *
                                   LLLL - 24.*pap2 * pbp1 * LR * RRRR + 24.*pap2 * pbp1 * LR * LRLR - 24.*pap1 *
                                   pbp2 * RL * RLRL + 24.*pap1 * pbp2 * RL * LLLL + 24.*pap1 * pbp2 * LR * RRRR -
                                   24.*pap1 * pbp2 * LR * LRLR - 8.*papb * p1p2 * RL * RLRL - 8.*papb * p1p2 * RL *
                                   LLLL - 8.*papb * p1p2 * LR * RRRR - 8.*papb * p1p2 * LR * LRLR - 8.*m1 * m2 *
                                   papb * RL * RLLL - 8.*m1 * m2 * papb * RL * LLRL - 8.*m1 * m2 * papb * LR * LRRR - 8
                                   * m1 * m2 * papb * LR * RRLR);

    bb->SetIeps(ieps + 2);// get two-point function prop. to O(1/epsilon^2).
    bb->GetNpf(la, 1);
    bb1 = la[0];
    delete bb;

    // Part of scalar integral of O(1/eps^1) times squared matrix element of O(eps^2).
    me2 = bb1 * ivs3v1 * ivs3v2 * (- 16.*pap2 * pbp1 * RL * RLRL + 16.*pap2 * pbp1 *
                                   RL * LLLL + 16.*pap2 * pbp1 * LR * RRRR - 16.*pap2 * pbp1 * LR * LRLR + 16.*
                                   pap1 * pbp2 * RL * RLRL - 16.*pap1 * pbp2 * RL * LLLL - 16.*pap1 * pbp2 * LR *
                                   RRRR + 16.*pap1 * pbp2 * LR * LRLR);

    return real(me0 + me1 + me2); // returns the finite result of the whole squared matrix element.
}


// All the other virtual corrections are computed in the same way!


double FI::MVsbu1t(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    me0 =   + bb0 * ivs3v1 * ivt1s2 * (- 4.*m2 * pbp1 * ml1 * RR * LLRR - 4.*m2 * pbp1 * ml1
                                       * LL * RRLL - 8.*m1 * pbp2 * ml1 * RR * RLRR - 8.*m1 * pbp2 * ml1 * LL * LRLL);
    me0 +=  + bb1 * ivs3v1 * ivt1s2 * (8.*pap1 * pbp2 * RL * RLRL + 8.*pap1 * pbp2
                                       * LR * LRLR + 4.*m1 * m2 * papb * RL * LLRL + 4.*m1 * m2 * papb * LR * RRLR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;

    me1 = + bb0 * ivs3v1 * ivt1s2 * (4.*m2 * pbp1 * ml1 * RR * LLRR + 4.*m2 * pbp1 * ml1 * LL
                                     * RRLL + 4.*m1 * pbp2 * ml1 * RR * RLRR + 4.*m1 * pbp2 * ml1 * LL * LRLL);
    me1 +=  + bb1 * ivs3v1 * ivt1s2 * (4.*pap2 * pbp1 * RL * RLRL + 4.*pap2 * pbp1
                                       * LR * LRLR - 4.*pap1 * pbp2 * RL * RLRL - 4.*pap1 * pbp2 * LR * LRLR - 4.*papb *
                                       p1p2 * RL * RLRL - 4.*papb * p1p2 * LR * LRLR - 4.*m1 * m2 * papb * RL * LLRL - 4.*
                                       m1 * m2 * papb * LR * RRLR);

    return real(me0 + me1 + me2);
}

double FI::MVsbu1u(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    me0 = + bb0 * ivs3v1 * ivu2s2 * (8.*m2 * pbp1 * ml1 * RR * LLRR + 8.*m2 * pbp1 * ml1 * LL
                                     * RRLL + 4.*m1 * pbp2 * ml1 * RR * RLRR + 4.*m1 * pbp2 * ml1 * LL * LRLL);
    me0 +=  + bb1 * ivs3v1 * ivu2s2 * (- 8.*pap2 * pbp1 * RL * LLRL - 8.*pap2 *
                                       pbp1 * LR * RRLR - 4.*m1 * m2 * papb * RL * RLRL - 4.*m1 * m2 * papb * LR * LRLR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;

    me1 = + bb0 * ivs3v1 * ivu2s2 * (- 4.*m2 * pbp1 * ml1 * RR * LLRR - 4.*m2 * pbp1 * ml1
                                     * LL * RRLL - 4.*m1 * pbp2 * ml1 * RR * RLRR - 4.*m1 * pbp2 * ml1 * LL * LRLL);
    me1 +=  + bb1 * ivs3v1 * ivu2s2 * (4.*pap2 * pbp1 * RL * LLRL + 4.*pap2 * pbp1
                                       * LR * RRLR - 4.*pap1 * pbp2 * RL * LLRL - 4.*pap1 * pbp2 * LR * RRLR + 4.*papb *
                                       p1p2 * RL * LLRL + 4.*papb * p1p2 * LR * RRLR + 4.*m1 * m2 * papb * RL * RLRL + 4.*
                                       m1 * m2 * papb * LR * LRLR);

    return real(me0 + me1 + me2);
}

double FI::MVsbu2s(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];

    me0 = + bb1 * ivs3v1 * ivs3v2 * (16.*pap2 * pbp1 * RL * LLLL + 16.*pap2 * pbp1 * LR *
                                     RRRR + 16.*pap1 * pbp2 * RL * RLRL + 16.*pap1 * pbp2 * LR * LRLR + 8.*m1 * m2 *
                                     papb * RL * RLLL + 8.*m1 * m2 * papb * RL * LLRL + 8.*m1 * m2 * papb * LR * LRRR + 8
                                     * m1 * m2 * papb * LR * RRLR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 1);
    bb1 = la[0];

    me1 = + bb1 * ivs3v1 * ivs3v2 * (24.*pap2 * pbp1 * RL * RLRL - 24.*pap2 * pbp1 * RL *
                                     LLLL - 24.*pap2 * pbp1 * LR * RRRR + 24.*pap2 * pbp1 * LR * LRLR - 24.*pap1 *
                                     pbp2 * RL * RLRL + 24.*pap1 * pbp2 * RL * LLLL + 24.*pap1 * pbp2 * LR * RRRR -
                                     24.*pap1 * pbp2 * LR * LRLR - 8.*papb * p1p2 * RL * RLRL - 8.*papb * p1p2 * RL *
                                     LLLL - 8.*papb * p1p2 * LR * RRRR - 8.*papb * p1p2 * LR * LRLR - 8.*m1 * m2 *
                                     papb * RL * RLLL - 8.*m1 * m2 * papb * RL * LLRL - 8.*m1 * m2 * papb * LR * LRRR - 8
                                     * m1 * m2 * papb * LR * RRLR);

    bb->SetIeps(ieps + 2);
    bb->GetNpf(la, 1);
    bb1 = la[0];
    delete bb;

    me2 = + bb1 * ivs3v1 * ivs3v2 * (- 16.*pap2 * pbp1 * RL * RLRL + 16.*pap2 * pbp1 *
                                     RL * LLLL + 16.*pap2 * pbp1 * LR * RRRR - 16.*pap2 * pbp1 * LR * LRLR + 16.*
                                     pap1 * pbp2 * RL * RLRL - 16.*pap1 * pbp2 * RL * LLLL - 16.*pap1 * pbp2 * LR *
                                     RRRR + 16.*pap1 * pbp2 * LR * LRLR);
    return real(me0 + me1 + me2);
}


double FI::MVsbu2t(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    me0 =   + bb0 * ivs3v1 * ivt1s2 * (- 8.*m2 * pap1 * ml1 * RR * LRRR - 8.*m2 * pap1 * ml1
                                       * LL * RLLL - 4.*m1 * pap2 * ml1 * RR * RRRR - 4.*m1 * pap2 * ml1 * LL * LLLL);
    me0 +=  + bb1 * ivs3v1 * ivt1s2 * (8.*pap1 * pbp2 * RL * RLRL + 8.*pap1 * pbp2
                                       * LR * LRLR + 4.*m1 * m2 * papb * RL * LLRL + 4.*m1 * m2 * papb * LR * RRLR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;

    me1 =   + bb0 * ivs3v1 * ivt1s2 * (4.*m2 * pap1 * ml1 * RR * LRRR + 4.*m2 * pap1 * ml1 * LL
                                       * RLLL + 4.*m1 * pap2 * ml1 * RR * RRRR + 4.*m1 * pap2 * ml1 * LL * LLLL);
    me1 +=  + bb1 * ivs3v1 * ivt1s2 * (4.*pap2 * pbp1 * RL * RLRL + 4.*pap2 * pbp1
                                       * LR * LRLR - 4.*pap1 * pbp2 * RL * RLRL - 4.*pap1 * pbp2 * LR * LRLR - 4.*papb *
                                       p1p2 * RL * RLRL - 4.*papb * p1p2 * LR * LRLR - 4.*m1 * m2 * papb * RL * LLRL - 4.*
                                       m1 * m2 * papb * LR * RRLR);
    return real(me0 + me1 + me2);
}

double FI::MVsbu2u(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    me0 =   + bb0 * ivs3v1 * ivu2s2 * (4.*m2 * pap1 * ml1 * RR * LRRR + 4.*m2 * pap1 * ml1 * LL
                                       * RLLL + 8.*m1 * pap2 * ml1 * RR * RRRR + 8.*m1 * pap2 * ml1 * LL * LLLL);
    me0 +=  + bb1 * ivs3v1 * ivu2s2 * (- 8.*pap2 * pbp1 * RL * LLRL - 8.*pap2 *
                                       pbp1 * LR * RRLR - 4.*m1 * m2 * papb * RL * RLRL - 4.*m1 * m2 * papb * LR * LRLR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;

    me1 =   + bb0 * ivs3v1 * ivu2s2 * (- 4.*m2 * pap1 * ml1 * RR * LRRR - 4.*m2 * pap1 * ml1
                                       * LL * RLLL - 4.*m1 * pap2 * ml1 * RR * RRRR - 4.*m1 * pap2 * ml1 * LL * LLLL);
    me1 +=  + bb1 * ivs3v1 * ivu2s2 * (4.*pap2 * pbp1 * RL * LLRL + 4.*pap2 * pbp1
                                       * LR * RRLR - 4.*pap1 * pbp2 * RL * LLRL - 4.*pap1 * pbp2 * LR * RRLR + 4.*papb *
                                       p1p2 * RL * LLRL + 4.*papb * p1p2 * LR * RRLR + 4.*m1 * m2 * papb * RL * RLRL + 4.*
                                       m1 * m2 * papb * LR * LRLR);
    return real(me0 + me1 + me2);
}

double FI::MVtbu1s(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    me0 =   + bb0 * ivs3v2 * ivt1s1 * (- 4.*m2 * pbp1 * ml1 * RR * LLRR - 4.*m2 * pbp1 * ml1
                                       * LL * RRLL - 8.*m1 * pbp2 * ml1 * RR * LLLR - 8.*m1 * pbp2 * ml1 * LL * RRRL);
    me0 +=  + bb1 * ivs3v2 * ivt1s1 * (8.*pap1 * pbp2 * RL * RLRL + 8.*pap1 * pbp2
                                       * LR * LRLR + 4.*m1 * m2 * papb * RL * RLLL + 4.*m1 * m2 * papb * LR * LRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;

    me1 =   + bb0 * ivs3v2 * ivt1s1 * (4.*m2 * pbp1 * ml1 * RR * LLRR + 4.*m2 * pbp1 * ml1 * LL
                                       * RRLL + 4.*m1 * pbp2 * ml1 * RR * LLLR + 4.*m1 * pbp2 * ml1 * LL * RRRL);
    me1 +=  + bb1 * ivs3v2 * ivt1s1 * (4.*pap2 * pbp1 * RL * RLRL + 4.*pap2 * pbp1
                                       * LR * LRLR - 4.*pap1 * pbp2 * RL * RLRL - 4.*pap1 * pbp2 * LR * LRLR - 4.*papb *
                                       p1p2 * RL * RLRL - 4.*papb * p1p2 * LR * LRLR - 4.*m1 * m2 * papb * RL * RLLL - 4.*
                                       m1 * m2 * papb * LR * LRRR);
    return real(me0 + me1 + me2);
}

double FI::MVtbu1t(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];
    delete bb;

    me0 =   + bb0 * ivt1s1 * ivt1s2 * (- 4.*m1 * pbp2 * ml1 * RR * RLRR - 4.*m1 * pbp2 * ml1
                                       * RR * LLLR - 4.*m1 * pbp2 * ml1 * LL * RRRL - 4.*m1 * pbp2 * ml1 * LL * LRLL);
    me0 +=  + bb1 * ivt1s1 * ivt1s2 * (4.*pap1 * pbp2 * RL * RLRL + 4.*pap1 * pbp2
                                       * RL * LLLL + 4.*pap1 * pbp2 * LR * RRRR + 4.*pap1 * pbp2 * LR * LRLR);

    return real(me0 + me1 + me2);
}

double FI::MVtbu1u(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    me0 =   + bb0 * ivt1s1 * ivu2s2 * (2.*m2 * pbp1 * ml1 * RR * LLLR + 2.*m2 * pbp1 * ml1 * LL
                                       * RRRL + 2.*m1 * pbp2 * ml1 * RR * RLRR + 2.*m1 * pbp2 * ml1 * LL * LRLL);
    me0 +=  + bb1 * ivt1s1 * ivu2s2 * (- 2.*pap2 * pbp1 * RL * LLLL - 2.*pap2 *
                                       pbp1 * LR * RRRR - 2.*pap1 * pbp2 * RL * LLLL - 2.*pap1 * pbp2 * LR * RRRR + 2.*
                                       papb * p1p2 * RL * LLLL + 2.*papb * p1p2 * LR * RRRR - 2.*m1 * m2 * papb * RL * RLRL
                                       - 2.*m1 * m2 * papb * LR * LRLR);
    delete bb;

    return real(me0 + me1 + me2);
}

double FI::MVtbu2s(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    me0 =   + bb0 * ivs3v2 * ivt1s1 * (- 8.*m2 * pap1 * ml1 * RR * LLRL - 8.*m2 * pap1 * ml1
                                       * LL * RRLR - 4.*m1 * pap2 * ml1 * RR * LLLL - 4.*m1 * pap2 * ml1 * LL * RRRR);
    me0 +=  + bb1 * ivs3v2 * ivt1s1 * (8.*pap1 * pbp2 * RL * RLRL + 8.*pap1 * pbp2
                                       * LR * LRLR + 4.*m1 * m2 * papb * RL * RLLL + 4.*m1 * m2 * papb * LR * LRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;

    me1 =   + bb0 * ivs3v2 * ivt1s1 * (4.*m2 * pap1 * ml1 * RR * LLRL + 4.*m2 * pap1 * ml1 * LL
                                       * RRLR + 4.*m1 * pap2 * ml1 * RR * LLLL + 4.*m1 * pap2 * ml1 * LL * RRRR);
    me1 +=  + bb1 * ivs3v2 * ivt1s1 * (4.*pap2 * pbp1 * RL * RLRL + 4.*pap2 * pbp1
                                       * LR * LRLR - 4.*pap1 * pbp2 * RL * RLRL - 4.*pap1 * pbp2 * LR * LRLR - 4.*papb *
                                       p1p2 * RL * RLRL - 4.*papb * p1p2 * LR * LRLR - 4.*m1 * m2 * papb * RL * RLLL - 4.*
                                       m1 * m2 * papb * LR * LRRR);

    return real(me0 + me1 + me2);
}

double FI::MVtbu2t(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];
    delete bb;

    me0 =   + bb0 * ivt1s1 * ivt1s2 * (- 4.*m2 * pap1 * ml1 * RR * LRRR - 4.*m2 * pap1 * ml1
                                       * RR * LLRL - 4.*m2 * pap1 * ml1 * LL * RRLR - 4.*m2 * pap1 * ml1 * LL * RLLL);
    me0 +=  + bb1 * ivt1s1 * ivt1s2 * (4.*pap1 * pbp2 * RL * RLRL + 4.*pap1 * pbp2
                                       * RL * RRRR + 4.*pap1 * pbp2 * LR * LRLR + 4.*pap1 * pbp2 * LR * LLLL);

    return real(me0 + me1 + me2);
}

double FI::MVtbu2u(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];
    delete bb;

    me0 =  + bb0 * ivt1s1 * ivu2s2 * (2.*m2 * pap1 * ml1 * RR * LRRR + 2.*m2 * pap1 * ml1 * LL
                                      * RLLL + 2.*m1 * pap2 * ml1 * RR * LLRL + 2.*m1 * pap2 * ml1 * LL * RRLR);
    me0 +=  + bb1 * ivt1s1 * ivu2s2 * (- 2.*pap2 * pbp1 * RL * RRRR - 2.*pap2 *
                                       pbp1 * LR * LLLL - 2.*pap1 * pbp2 * RL * RRRR - 2.*pap1 * pbp2 * LR * LLLL + 2.*
                                       papb * p1p2 * RL * RRRR + 2.*papb * p1p2 * LR * LLLL - 2.*m1 * m2 * papb * RL * RLRL
                                       - 2.*m1 * m2 * papb * LR * LRLR);

    return real(me0 + me1 + me2);
}

double FI::MVubu1s(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    me0 =   + bb0 * ivs3v2 * ivu2s1 * (8.*m2 * pbp1 * ml1 * RR * LLRR + 8.*m2 * pbp1 * ml1 * LL
                                       * RRLL + 4.*m1 * pbp2 * ml1 * RR * LLLR + 4.*m1 * pbp2 * ml1 * LL * RRRL);
    me0 +=  + bb1 * ivs3v2 * ivu2s1 * (- 8.*pap2 * pbp1 * RL * RLLL - 8.*pap2 *
                                       pbp1 * LR * LRRR - 4.*m1 * m2 * papb * RL * RLRL - 4.*m1 * m2 * papb * LR * LRLR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;

    me1 =   + bb0 * ivs3v2 * ivu2s1 * (- 4.*m2 * pbp1 * ml1 * RR * LLRR - 4.*m2 * pbp1 * ml1
                                       * LL * RRLL - 4.*m1 * pbp2 * ml1 * RR * LLLR - 4.*m1 * pbp2 * ml1 * LL * RRRL);
    me1 +=  + bb1 * ivs3v2 * ivu2s1 * (4.*pap2 * pbp1 * RL * RLLL + 4.*pap2 * pbp1
                                       * LR * LRRR - 4.*pap1 * pbp2 * RL * RLLL - 4.*pap1 * pbp2 * LR * LRRR + 4.*papb *
                                       p1p2 * RL * RLLL + 4.*papb * p1p2 * LR * LRRR + 4.*m1 * m2 * papb * RL * RLRL + 4.*
                                       m1 * m2 * papb * LR * LRLR);

    return real(me0 + me1 + me2);
}

double FI::MVubu1t(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];
    delete bb;

    me0 =   + bb0 * ivt1s2 * ivu2s1 * (2.*m2 * pbp1 * ml1 * RR * RLRR + 2.*m2 * pbp1 * ml1 * LL
                                       * LRLL + 2.*m1 * pbp2 * ml1 * RR * LLLR + 2.*m1 * pbp2 * ml1 * LL * RRRL);
    me0 +=  + bb1 * ivt1s2 * ivu2s1 * (- 2.*pap2 * pbp1 * RL * LLLL - 2.*pap2 *
                                       pbp1 * LR * RRRR - 2.*pap1 * pbp2 * RL * LLLL - 2.*pap1 * pbp2 * LR * RRRR + 2.*
                                       papb * p1p2 * RL * LLLL + 2.*papb * p1p2 * LR * RRRR - 2.*m1 * m2 * papb * RL * RLRL
                                       - 2.*m1 * m2 * papb * LR * LRLR);

    return real(me0 + me1 + me2);
}

double FI::MVubu1u(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];
    delete bb;

    me0 =  + bb0 * ivu2s1 * ivu2s2 * (- 4.*m2 * pbp1 * ml1 * RR * RLRR - 4.*m2 * pbp1 * ml1
                                      * RR * LLLR - 4.*m2 * pbp1 * ml1 * LL * RRRL - 4.*m2 * pbp1 * ml1 * LL * LRLL);
    me0 +=  + bb1 * ivu2s1 * ivu2s2 * (4.*pap2 * pbp1 * RL * RLRL + 4.*pap2 * pbp1
                                       * RL * LLLL + 4.*pap2 * pbp1 * LR * RRRR + 4.*pap2 * pbp1 * LR * LRLR);

    return real(me0 + me1 + me2);
}

double FI::MVubu2s(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    me0 =  + bb0 * ivs3v2 * ivu2s1 * (4.*m2 * pap1 * ml1 * RR * LLRL + 4.*m2 * pap1 * ml1 * LL
                                      * RRLR + 8.*m1 * pap2 * ml1 * RR * LLLL + 8.*m1 * pap2 * ml1 * LL * RRRR);
    me0 +=  + bb1 * ivs3v2 * ivu2s1 * (- 8.*pap2 * pbp1 * RL * RLLL - 8.*pap2 *
                                       pbp1 * LR * LRRR - 4.*m1 * m2 * papb * RL * RLRL - 4.*m1 * m2 * papb * LR * LRLR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;

    me1 =   + bb0 * ivs3v2 * ivu2s1 * (- 4.*m2 * pap1 * ml1 * RR * LLRL - 4.*m2 * pap1 * ml1
                                       * LL * RRLR - 4.*m1 * pap2 * ml1 * RR * LLLL - 4.*m1 * pap2 * ml1 * LL * RRRR);
    me1 +=  + bb1 * ivs3v2 * ivu2s1 * (4.*pap2 * pbp1 * RL * RLLL + 4.*pap2 * pbp1
                                       * LR * LRRR - 4.*pap1 * pbp2 * RL * RLLL - 4.*pap1 * pbp2 * LR * LRRR + 4.*papb *
                                       p1p2 * RL * RLLL + 4.*papb * p1p2 * LR * LRRR + 4.*m1 * m2 * papb * RL * RLRL + 4.*
                                       m1 * m2 * papb * LR * LRLR);
    return real(me0 + me1 + me2);
}

double FI::MVubu2t(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];
    delete bb;

    me0 =   + bb0 * ivt1s2 * ivu2s1 * (2.*m2 * pap1 * ml1 * RR * LLRL + 2.*m2 * pap1 * ml1 * LL
                                       * RRLR + 2.*m1 * pap2 * ml1 * RR * LRRR + 2.*m1 * pap2 * ml1 * LL * RLLL);
    me0 +=  + bb1 * ivt1s2 * ivu2s1 * (- 2.*pap2 * pbp1 * RL * RRRR - 2.*pap2 *
                                       pbp1 * LR * LLLL - 2.*pap1 * pbp2 * RL * RRRR - 2.*pap1 * pbp2 * LR * LLLL + 2.*
                                       papb * p1p2 * RL * RRRR + 2.*papb * p1p2 * LR * LLLL - 2.*m1 * m2 * papb * RL * RLRL
                                       - 2.*m1 * m2 * papb * LR * LRLR);

    return real(me0 + me1 + me2);
}

double FI::MVubu2u(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];
    delete bb;

    me0 =  + bb0 * ivu2s1 * ivu2s2 * (- 4.*m1 * pap2 * ml1 * RR * LRRR - 4.*m1 * pap2 * ml1
                                      * RR * LLRL - 4.*m1 * pap2 * ml1 * LL * RRLR - 4.*m1 * pap2 * ml1 * LL * RLLL);
    me0 +=  + bb1 * ivu2s1 * ivu2s2 * (4.*pap2 * pbp1 * RL * RLRL + 4.*pap2 * pbp1
                                       * RL * RRRR + 4.*pap2 * pbp1 * LR * LRLR + 4.*pap2 * pbp1 * LR * LLLL);

    return real(me0 + me1 + me2);
}

double FI::MVtbu3s(int ieps)
{
    Npf * bb = new Npf(m1s - 2.*pap1, m1l, m2l, mul, ieps);
    complex<double> la[7];
    bb->GetNpf(la, 2);
    complex<double> bb21 = la[0];
    complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    complex<double> me0 = 4.*ivs3v2 * pow(ivt1s1, 2) *
                          (4.*bb22 + (bb0 - 2.*bb1 + bb21) * (m1s - 2.*pap1)) *
                          (m1 * m2 * papb * (LRRR + RLLL) + 2.*pap1 * pbp2 * (LRLR + RLRL)) * RR;

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb21 = la[0];
    bb22 = la[1];
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];

    complex<double> me1 = -4.*ivs3v2 * pow(ivt1s1, 2) *
                          ((bb0 - 2.*bb1 + bb21) * (m1s - 2.*pap1) *
                           (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - LRLR * pap2 * pbp1 + LRLR * pap1 * pbp2 +
                            m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL) +
                           bb22 * (6.*LRRR * m1 * m2 * papb + 4.*LRLR * p1p2 * papb - 4.*LRLR * pap2 * pbp1 +
                                   8.*LRLR * pap1 * pbp2 + 6.*m1 * m2 * papb * RLLL + 4.*p1p2 * papb * RLRL -
                                   4.*pap2 * pbp1 * RLRL + 8.*pap1 * pbp2 * RLRL)) * RR;

    bb->SetIeps(ieps + 2);
    bb->GetNpf(la, 2);
    bb21 = la[0];
    bb22 = la[1];
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;

    complex<double> me2 = 8.*bb22 * ivs3v2 * pow(ivt1s1, 2) *
                          (LRRR * m1 * m2 * papb +
                           LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLLL +
                           p1p2 * papb * RLRL - pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL) * RR;

    return real(me0 + me1 + me2);
}

double FI::MVtbu3t(int ieps)
{
    Npf * bb = new Npf(m1s - 2.*pap1, m1l, m2l, mul, ieps);
    complex<double> la[7];
    bb->GetNpf(la, 2);
    complex<double> bb21 = la[0];
    complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    complex<double> me0 = 4.*pow(ivt1s1, 2) * ivt1s2 *
                          (4.*bb22 + (bb0 - 2.*bb1 + bb21) * (m1s - 2.*pap1)) * pap1 * pbp2 * RR *
                          (LLLL + LRLR + RLRL + RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    complex<double> me1 = -8.*bb22 * pow(ivt1s1, 2) * ivt1s2 * pap1 * pbp2
                          * RR * (LLLL + LRLR + RLRL + RRRR);

    return real(me0 + me1);
}

double FI::MVtbu3u(int ieps)
{
    Npf * bb = new Npf(m1s - 2.*pap1, m1l, m2l, mul, ieps);
    complex<double> la[7];
    bb->GetNpf(la, 2);
    complex<double> bb21 = la[0];
    complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    complex<double> me0 = -2.*pow(ivt1s1, 2) * ivu2s2 *
                          (4.*bb22 + (bb0 - 2.*bb1 + bb21) * (m1s - 2.*pap1)) * RR *
                          (LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                           m1 * m2 * papb * RLRL - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    complex<double> me1 = 4.*bb22 * pow(ivt1s1, 2) * ivu2s2 *
                          RR * (LRLR * m1 * m2 * papb +
                                LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL -
                                p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    return real(me0 + me1);
}

double FI::MVubu3s(int ieps)
{
    Npf * bb = new Npf(m2s - 2.*pap2, m1l, m2l, mul, ieps);
    complex<double> la[7];
    bb->GetNpf(la, 2);
    complex<double> bb21 = la[0];
    complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    complex<double> me0 = -4.*ivs3v2 * pow(ivu2s1, 2) *
                          (4.*bb22 + (bb0 - 2.*bb1 + bb21) * (m2s - 2.*pap2)) *
                          (2.*pap2 * pbp1 * (LRRR + RLLL) + m1 * m2 * papb * (LRLR + RLRL)) * RR;

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb21 = la[0];
    bb22 = la[1];
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];

    complex<double> me1 = 4.*ivs3v2 * pow(ivu2s1, 2) *
                          ((bb0 - 2.*bb1 + bb21) * (m2s - 2.*pap2) *
                           (LRLR * m1 * m2 * papb + (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) +
                            m1 * m2 * papb * RLRL) + bb22 * (6.*LRLR * m1 * m2 * papb +
                                    4.*(p1p2 * papb + 2.*pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) +
                                    6.*m1 * m2 * papb * RLRL)) * RR;

    bb->SetIeps(ieps + 2);
    bb->GetNpf(la, 2);
    bb21 = la[0];
    bb22 = la[1];
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;

    complex<double> me2 = -8.*bb22 * ivs3v2 * pow(ivu2s1, 2) *
                          (LRLR * m1 * m2 * papb + (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) *
                           (LRRR + RLLL) + m1 * m2 * papb * RLRL) * RR;

    return real(me0 + me1 + me2);
}

double FI::MVubu3t(int ieps)
{
    Npf * bb = new Npf(m2s - 2.*pap2, m1l, m2l, mul, ieps);
    complex<double> la[7];
    bb->GetNpf(la, 2);
    complex<double> bb21 = la[0];
    complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    complex<double> me0 = -2.*ivt1s2 * pow(ivu2s1, 2) *
                          (4.*bb22 + (bb0 - 2.*bb1 + bb21) * (m2s - 2.*pap2)) * RR *
                          (LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                           m1 * m2 * papb * RLRL - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    complex<double> me1 = 4.*bb22 * ivt1s2 * pow(ivu2s1, 2) *
                          RR * (LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                                m1 * m2 * papb * RLRL - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    return real(me0 + me1);
}

double FI::MVubu3u(int ieps)
{
    Npf * bb = new Npf(m2s - 2.*pap2, m1l, m2l, mul, ieps);
    complex<double> la[7];
    bb->GetNpf(la, 2);
    complex<double> bb21 = la[0];
    complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    complex<double> me0 = 4.*pow(ivu2s1, 2) * ivu2s2 *
                          (4.*bb22 + (bb0 - 2.*bb1 + bb21) * (m2s - 2.*pap2)) * pap2 * pbp1 * RR *
                          (LLLL + LRLR + RLRL + RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    complex<double> me1 = -8.*bb22 * pow(ivu2s1, 2) * ivu2s2 * pap2 * pbp1 *
                          RR * (LLLL + LRLR + RLRL + RRRR);

    return real(me0 + me1);
}

double FI::MVtbu4s(int ieps)
{
    ivt1s1 = 1. / (m1s - 2.*pap1 - m3l);
    ivt2s1 = 1. / (m1s - 2.*pap1 - m4l);

    Npf * bb = new Npf(m1s - 2.*pap1, m1l, m2l, mul, ieps);
    complex<double> la[7];
    bb->GetNpf(la, 2);
    complex<double> bb21 = la[0];
    complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    complex<double> me0 = 8.*ivs3v2 * ivt1s1 * ivt2s1 *
                          (m1 * m2 * papb * (LRRR + RLLL) + 2.*pap1 * pbp2 * (LRLR + RLRL)) *
                          ((4.*bb22 + (bb1 + bb21) * (m1s - 2.*pap1)) * (LR + RL) +
                           bb0 * ml1 * ml2 * (LL + RR));

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb21 = la[0];
    bb22 = la[1];
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];

    complex<double> me1 = -8.*ivs3v2 * ivt1s1 * ivt2s1 *
                          (2.*bb22 * (LR + RL) * (3.*LRRR * m1 * m2 * papb + 2.*LRLR * p1p2 * papb -
                                  2.*LRLR * pap2 * pbp1 + 4.*LRLR * pap1 * pbp2 + 3.*m1 * m2 * papb * RLLL +
                                  2.*p1p2 * papb * RLRL - 2.*pap2 * pbp1 * RLRL + 4.*pap1 * pbp2 * RLRL) +
                           (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - LRLR * pap2 * pbp1 + LRLR * pap1 * pbp2 +
                            m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL) *
                           ((bb1 + bb21) * (m1s - 2.*pap1) * (LR + RL) + bb0 * ml1 * ml2 * (LL + RR)));

    bb->SetIeps(ieps + 2);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    complex<double> me2 = 16.*bb22 * ivs3v2 * ivt1s1 * ivt2s1 *
                          (LR + RL) * (LRRR * m1 * m2 * papb + LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)
                                       + m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL);

    return real(me0 + me1 + me2);
}

double FI::MVtbu4t(int ieps)
{
    ivt1s1 = 1. / (m1s - 2.*pap1 - m3l);
    ivt2s1 = 1. / (m1s - 2.*pap1 - m4l);

    Npf * bb = new Npf(m1s - 2.*pap1, m1l, m2l, mul, ieps);
    complex<double> la[7];
    bb->GetNpf(la, 2);
    complex<double> bb21 = la[0];
    complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    complex<double> me0 = 8.*ivt1s1 * ivt2s1 * ivt1s2 * pap1 * pbp2 *
                          ((4.*bb22 + (bb1 + bb21) * (m1s - 2.*pap1)) *
                           (LR + RL) + bb0 * ml1 * ml2 * (LL + RR)) * (LLLL + LRLR + RLRL + RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    complex<double> me1 = -16.*bb22 * ivt1s1 * ivt2s1 * ivt1s2 *
                          pap1 * pbp2 * (LR + RL) * (LLLL + LRLR + RLRL + RRRR);

    return real(me0 + me1);
}

double FI::MVtbu4u(int ieps)
{
    ivt1s1 = 1. / (m1s - 2.*pap1 - m3l);
    ivt2s1 = 1. / (m1s - 2.*pap1 - m4l);

    Npf * bb = new Npf(m1s - 2.*pap1, m1l, m2l, mul, ieps);
    complex<double> la[7];
    bb->GetNpf(la, 2);
    complex<double> bb21 = la[0];
    complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    complex<double> me0 = -4.*ivt1s1 * ivt2s1 * ivu2s2 *
                          ((4.*bb22 + (bb1 + bb21) * (m1s - 2.*pap1)) * (LR + RL) +
                           bb0 * ml1 * ml2 * (LL + RR)) * (LRLR * m1 * m2 * papb +
                                   LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL -
                                   p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    complex<double> me1 = 8.*bb22 * ivt1s1 * ivt2s1 * ivu2s2 *
                          (LR + RL) * (LRLR * m1 * m2 * papb +
                                       LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL -
                                       p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    return real(me0 + me1);
}

double FI::MVubu4s(int ieps)
{
    ivu1s1 = 1. / (m2s - 2.*pap2 - m3l);
    ivu2s1 = 1. / (m2s - 2.*pap2 - m4l);

    Npf * bb = new Npf(m2s - 2.*pap2, m1l, m2l, mul, ieps);
    complex<double> la[7];
    bb->GetNpf(la, 2);
    complex<double> bb21 = la[0];
    complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    complex<double> me0 = -8.*ivs3v2 * ivu1s1 * ivu2s1 *
                          (2.*pap2 * pbp1 * (LRRR + RLLL) + m1 * m2 * papb * (LRLR + RLRL)) *
                          ((4.*bb22 + (bb1 + bb21) * (m2s - 2.*pap2)) * (LR + RL) + bb0 * ml1 * ml2 * (LL + RR));

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb21 = la[0];
    bb22 = la[1];
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];

    complex<double> me1 = 8.*ivs3v2 * ivu1s1 * ivu2s1 *
                          (2.*bb22 * (LR + RL) * (3.*LRLR * m1 * m2 * papb +
                                  2.*(p1p2 * papb + 2.*pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) +
                                  3.*m1 * m2 * papb * RLRL) +
                           (LRLR * m1 * m2 * papb + (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) +
                            m1 * m2 * papb * RLRL) * ((bb1 + bb21) * (m2s - 2.*pap2) * (LR + RL) +
                                    bb0 * ml1 * ml2 * (LL + RR)));

    bb->SetIeps(ieps + 2);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    complex<double> me2 = -16.*bb22 * ivs3v2 * ivu1s1 * ivu2s1 *
                          (LR + RL) * (LRLR * m1 * m2 * papb +
                                       (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) + m1 * m2 * papb * RLRL);

    return real(me0 + me1 + me2);
}

double FI::MVubu4t(int ieps)
{
    ivu1s1 = 1. / (m2s - 2.*pap2 - m3l);
    ivu2s1 = 1. / (m2s - 2.*pap2 - m4l);

    Npf * bb = new Npf(m2s - 2.*pap2, m1l, m2l, mul, ieps);
    complex<double> la[7];
    bb->GetNpf(la, 2);
    complex<double> bb21 = la[0];
    complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    complex<double> me0 = -4.*ivt1s2 * ivu1s1 * ivu2s1 *
                          ((4.*bb22 + (bb1 + bb21) * (m2s - 2.*pap2)) * (LR + RL) +
                           bb0 * ml1 * ml2 * (LL + RR)) * (LRLR * m1 * m2 * papb +
                                   LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL -
                                   p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    complex<double> me1 = 8.*bb22 * ivt1s2 * ivu1s1 * ivu2s1 *
                          (LR + RL) * (LRLR * m1 * m2 * papb +
                                       LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL -
                                       p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    return real(me0 + me1);
}

double FI::MVubu4u(int ieps)
{
    ivu1s1 = 1. / (m2s - 2.*pap2 - m3l);
    ivu2s1 = 1. / (m2s - 2.*pap2 - m4l);

    Npf * bb = new Npf(m2s - 2.*pap2, m1l, m2l, mul, ieps);
    complex<double> la[7];
    bb->GetNpf(la, 2);
    complex<double> bb21 = la[0];
    complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    complex<double> bb0 = la[0];

    complex<double> me0 = 8.*ivu1s1 * ivu2s1 * ivu2s2 *
                          pap2 * pbp1 * ((4.*bb22 + (bb1 + bb21) * (m2s - 2.*pap2)) *
                                         (LR + RL) + bb0 * ml1 * ml2 * (LL + RR)) * (LLLL + LRLR + RLRL + RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    complex<double> me1 = -16.*bb22 * ivu1s1 * ivu2s1 * ivu2s2 *
                          pap2 * pbp1 * (LR + RL) * (LLLL + LRLR + RLRL + RRRR);

    return real(me0 + me1);
}

double FI::MVtbu5s(int ieps)
{
    ivt1s1 = 1. / (m1s - 2.*pap1 - m3l);
    ivt2s1 = 1. / (m1s - 2.*pap1 - m4l);

    Npf * aa = new Npf(m1l, mul, ieps);
    complex<double> la[7];
    aa->GetNpf(la, 0);
    complex<double> aa0 = la[0];

    complex<double> me0 = aa0 * ivs3v2 * ivt1s1 * ivt2s1 *
                          (8.*RR * pap1 * pbp2 * RLRL + 8.*RR * pap1 * pbp2 * LRLR +
                           4.*RR * m1 * m2 * papb * LRRR + 4.*RR * m1 * m2 * papb * RLLL);

    aa->SetIeps(ieps + 1);
    aa->GetNpf(la, 0);
    aa0 = la[0];
    delete aa;

    complex<double> me1 = aa0 * ivs3v2 * ivt1s1 * ivt2s1 *
                          (4.*RR * pap2 * pbp1 * RLRL + 4.*RR * pap2 * pbp1 * LRLR - 4.*RR * pap1 * pbp2 * RLRL -
                           4.*RR * pap1 * pbp2 * LRLR - 4.*RR * papb * p1p2 * RLRL - 4.*RR * papb * p1p2 * LRLR -
                           4.*RR * m1 * m2 * papb * LRRR - 4.*RR * m1 * m2 * papb * RLLL);

    return real(me0 + me1);
}

double FI::MVtbu5t(int ieps)
{
    ivt1s1 = 1. / (m1s - 2.*pap1 - m3l);
    ivt2s1 = 1. / (m1s - 2.*pap1 - m4l);

    Npf * aa = new Npf(m1l, mul, ieps);
    complex<double> la[7];
    aa->GetNpf(la, 0);
    complex<double> aa0 = la[0];
    delete aa;

    complex<double> me0 = aa0 * ivt1s1 * ivt2s1 * ivt1s2 *
                          (4.*RR * pap1 * pbp2 * RLRL + 4.*RR * pap1 *
                           pbp2 * RRRR + 4.*RR * pap1 * pbp2 * LRLR + 4.*RR * pap1 * pbp2 * LLLL);

    return real(me0);
}

double FI::MVtbu5u(int ieps)
{
    ivt1s1 = 1. / (m1s - 2.*pap1 - m3l);
    ivt2s1 = 1. / (m1s - 2.*pap1 - m4l);

    Npf * aa = new Npf(m1l, mul, ieps);
    complex<double> la[7];
    aa->GetNpf(la, 0);
    complex<double> aa0 = la[0];
    delete aa;

    complex<double> me0 = aa0 * ivt1s1 * ivt2s1 * ivu2s2 *
                          (- 2.*RR * pap2 * pbp1 * RRRR - 2.*RR * pap2 * pbp1 * LLLL
                           - 2.*RR * pap1 * pbp2 * RRRR - 2.*RR * pap1 * pbp2 * LLLL
                           + 2.*RR * papb * p1p2 * RRRR + 2.*RR * papb * p1p2 * LLLL
                           - 2.*RR * m1 * m2 * papb * RLRL - 2.*RR * m1 * m2 * papb * LRLR);

    return real(me0);
}

double FI::MVubu5s(int ieps)
{
    ivu1s1 = 1. / (m2s - 2.*pap2 - m3l);
    ivu2s1 = 1. / (m2s - 2.*pap2 - m4l);

    Npf * aa = new Npf(m1l, mul, ieps);
    complex<double> la[7];
    aa->GetNpf(la, 0);
    complex<double> aa0 = la[0];

    complex<double> me0 = aa0 * ivs3v2 * ivu1s1 * ivu2s1 *
                          (- 8.*RR * pap2 * pbp1 * LRRR - 8.*RR * pap2 * pbp1 * RLLL -
                           4.*RR * m1 * m2 * papb * RLRL - 4.*RR * m1 * m2 * papb * LRLR);

    aa->SetIeps(ieps + 1);
    aa->GetNpf(la, 0);
    aa0 = la[0];
    delete aa;

    complex<double> me1 = aa0 * ivs3v2 * ivu1s1 * ivu2s1 *
                          (4.*RR * pap2 * pbp1 * LRRR + 4.*RR * pap2 * pbp1 * RLLL -
                           4.*RR * pap1 * pbp2 * LRRR - 4.*RR * pap1 * pbp2 * RLLL +
                           4.*RR * papb * p1p2 * LRRR + 4.*RR * papb * p1p2 * RLLL +
                           4.*RR * m1 * m2 * papb * RLRL + 4.*RR * m1 * m2 * papb * LRLR);

    return real(me0 + me1);
}

double FI::MVubu5t(int ieps)
{
    ivu1s1 = 1. / (m2s - 2.*pap2 - m3l);
    ivu2s1 = 1. / (m2s - 2.*pap2 - m4l);

    Npf * aa = new Npf(m1l, mul, ieps);
    complex<double> la[7];
    aa->GetNpf(la, 0);
    complex<double> aa0 = la[0];
    delete aa;

    complex<double> me0 = aa0 * ivt1s2 * ivu1s1 * ivu2s1 *
                          (- 2.*RR * pap2 * pbp1 * RRRR - 2.*RR *
                           pap2 * pbp1 * LLLL - 2.*RR * pap1 * pbp2 * RRRR - 2.*RR * pap1 * pbp2 * LLLL +
                           2.*RR * papb * p1p2 * RRRR + 2.*RR * papb * p1p2 * LLLL -
                           2.*RR * m1 * m2 * papb * RLRL - 2.*RR * m1 * m2 * papb * LRLR);

    return real(me0);
}

double FI::MVubu5u(int ieps)
{
    ivu1s1 = 1. / (m2s - 2.*pap2 - m3l);
    ivu2s1 = 1. / (m2s - 2.*pap2 - m4l);

    Npf * aa = new Npf(m1l, mul, ieps);
    complex<double> la[7];
    aa->GetNpf(la, 0);
    complex<double> aa0 = la[0];
    delete aa;

    complex<double> me0 = aa0 * ivu1s1 * ivu2s1 * ivu2s2 *
                          (4.*RR * pap2 * pbp1 * RLRL + 4.*RR * pap2 *
                           pbp1 * RRRR + 4.*RR * pap2 * pbp1 * LRLR + 4.*RR * pap2 * pbp1 * LLLL);


    return real(me0);
}


// Process: q + \bar{q} -> sl + \bar{sl}.

// For the squared matrix element in D = 4 - epsilon dimensions there are only
// terms of order O(eps^0).
double FI::MVsbu1sSL(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    delete bb;

    me0 = bb1 * (8.*ivs3v1 * ivs3v2 * (LRLR * LR + LLLL * RL) * (2.*pap1 * pap2 - pap1 * m2s - pap2 * m1s));

    return real(me0);
}

double FI::MVsbu2sSL(int ieps)
{
    complex<double> me0 = complex<double> (0., 0.);
    complex<double> me1 = complex<double> (0., 0.);
    complex<double> me2 = complex<double> (0., 0.);

    complex<double> la[7];
    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    bb->GetNpf(la, 1);
    complex<double> bb1 = la[0];
    delete bb;

    me0 = bb1 * (8.*ivs3v1 * ivs3v2 * (LRLR * LR + LLLL * RL) * (2.*pap1 * pap2 - pap1 * m2s - pap2 * m1s));

    return real(me0);
}
