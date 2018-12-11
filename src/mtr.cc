// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Triangle matrix elements.
//
// See mbu.cc and mbn.cc for notes on notation.

#include <complex>
#include "fi.h"
#include "npf.h"

using namespace std;

// Process: q + \bar{q} -> f + \bar{f}.

double FI::MVstr1s(int ieps)
{
    // Constructs general three point function:
    // C(pa^2, pb^2,(pa+pb)^2,  loop mass1, loop mass2, loop mass3, renormalization scale,
    // term proportional to a certain order of 1/epsilon ).
    Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    // Part of scalar integral of O(1/eps^0) times squared matrix element of O(eps^0).
    complex<double> me0 = 32.*ivs3v1 * ivs3v2 * (cc24 + (cc11 + cc23) * papb) *
                          (LL * (m1 * m2 * papb * (LLRL + RLLL) + 2.*(LLLL * pap2 * pbp1 + pap1 * pbp2 * RLRL)) +
                           RR * (m1 * m2 * papb * (LRRR + RRLR) + 2.*(LRLR * pap1 * pbp2 + pap2 * pbp1 * RRRR)));

    cc->SetIeps(ieps + 1);  // get three-point function prop. to O(1/epsilon^1).
    cc->GetNpf(la, 2);
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
    cc11 = la[0];
    cc12 = la[1];

    // Part of scalar integral of O(1/eps^1) times squared matrix element of O(eps^1).
    complex<double> me1 = -32.*ivs3v1 * ivs3v2 *
                          (cc24 * (LL * (3.*LLRL * m1 * m2 * papb + LLLL * p1p2 * papb + 7.*LLLL * pap2 * pbp1 -
                                         3.*LLLL * pap1 * pbp2 + 3.*m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                                         3.*pap2 * pbp1 * RLRL + 7.*pap1 * pbp2 * RLRL) +
                                   RR * (3.*LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - 3.*LRLR * pap2 * pbp1 +
                                         7.*LRLR * pap1 * pbp2 + 3.*m1 * m2 * papb * RRLR + p1p2 * papb * RRRR +
                                         7.*pap2 * pbp1 * RRRR - 3.*pap1 * pbp2 * RRRR)) +
                           papb * (cc11 * (LL * (LLRL * m1 * m2 * papb + LLLL * p1p2 * papb + 3.*LLLL * pap2 * pbp1 -
                                           3.*LLLL * pap1 * pbp2 + m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                                           3.*pap2 * pbp1 * RLRL + 3.*pap1 * pbp2 * RLRL) +
                                           RR * (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - 3.*LRLR * pap2 * pbp1 +
                                                   3.*LRLR * pap1 * pbp2 + m1 * m2 * papb * RRLR + p1p2 * papb * RRRR +
                                                   3.*pap2 * pbp1 * RRRR - 3.*pap1 * pbp2 * RRRR)) +
                                   cc23 * (LL * (2.*LLRL * m1 * m2 * papb + LLLL * p1p2 * papb + 5.*LLLL * pap2 * pbp1 -
                                           3.*LLLL * pap1 * pbp2 + 2.*m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                                           3.*pap2 * pbp1 * RLRL + 5.*pap1 * pbp2 * RLRL) +
                                           RR * (2.*LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - 3.*LRLR * pap2 * pbp1 +
                                                   5.*LRLR * pap1 * pbp2 + 2.*m1 * m2 * papb * RRLR + p1p2 * papb * RRRR +
                                                   5.*pap2 * pbp1 * RRRR - 3.*pap1 * pbp2 * RRRR)) +
                                   cc12 * (LL * (m1 * m2 * papb * (LLRL + RLLL) + 2.*(LLLL * pap2 * pbp1 +
                                           pap1 * pbp2 * RLRL)) + RR * (m1 * m2 * papb * (LRRR + RRLR) +
                                                   2.*(LRLR * pap1 * pbp2 + pap2 * pbp1 * RRRR)))));

    cc->SetIeps(ieps + 2); // Gets three-point function prop. to O(1/epsilon^2).
    cc->GetNpf(la, 2);
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
    cc11 = la[0];
    cc12 = la[1];
    delete cc;

    // Part of scalar integral of O(1/eps^2) times squared matrix element of O(eps^2).
    complex<double> me2 = 32.*ivs3v1 * ivs3v2 *
                          (cc24 * (LL * (3.*LLRL * m1 * m2 * papb + 2.*LLLL * (p1p2 * papb + 5.*pap2 * pbp1 -
                                         4.*pap1 * pbp2) + 3.*m1 * m2 * papb * RLLL +
                                         2.*(p1p2 * papb - 4.*pap2 * pbp1 + 5.*pap1 * pbp2) * RLRL) +
                                   RR * (3.*LRRR * m1 * m2 * papb + 2.*LRLR * (p1p2 * papb - 4.*pap2 * pbp1 + 5.*pap1 * pbp2) +
                                         3.*m1 * m2 * papb * RRLR + 2.*(p1p2 * papb + 5.*pap2 * pbp1 - 4.*pap1 * pbp2) * RRRR)) +
                           papb * (2.*cc11 * (pap2 * pbp1 - pap1 * pbp2) * (LL * (LLLL - RLRL) +
                                   RR * (-LRLR + RRRR)) +
                                   cc23 * (LL * (LLRL * m1 * m2 * papb + LLLL * p1p2 * papb + 5.*LLLL * pap2 * pbp1 -
                                           5.*LLLL * pap1 * pbp2 + m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                                           5.*pap2 * pbp1 * RLRL + 5.*pap1 * pbp2 * RLRL) +
                                           RR * (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - 5.*LRLR * pap2 * pbp1 +
                                                   5.*LRLR * pap1 * pbp2 + m1 * m2 * papb * RRLR + p1p2 * papb * RRRR +
                                                   5.*pap2 * pbp1 * RRRR - 5.*pap1 * pbp2 * RRRR)) +
                                   cc12 * (LL * (LLRL * m1 * m2 * papb + LLLL * p1p2 * papb + 3.*LLLL * pap2 * pbp1 -
                                           3.*LLLL * pap1 * pbp2 + m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                                           3.*pap2 * pbp1 * RLRL + 3.*pap1 * pbp2 * RLRL) +
                                           RR * (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - 3.*LRLR * pap2 * pbp1 +
                                                   3.*LRLR * pap1 * pbp2 + m1 * m2 * papb * RRLR + p1p2 * papb * RRRR +
                                                   3.*pap2 * pbp1 * RRRR - 3.*pap1 * pbp2 * RRRR))));

    // Returns the finite result of the whole squared matrix element.
    return real(me0 + me1 + me2);
}

double FI::MVstr1t(int ieps)
{
    Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = 16.*ivs3v1 * ivt1s2 * (cc24 + (cc11 + cc23) * papb) *
                          (LL * LLRL * m1 * m2 * papb + 2.*LL * pap1 * pbp2 * RLRL + 2.*LRLR * pap1 * pbp2 * RR +
                           m1 * m2 * papb * RR * RRLR);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
    cc11 = la[0];
    cc12 = la[1];

    complex<double> me1 = -16.*ivs3v1 * ivt1s2 *
                          (cc24 * (LL * (3.*LLRL * m1 * m2 * papb + (p1p2 * papb - pap2 * pbp1 + 5.*pap1 * pbp2) * RLRL) +
                                   RR * (LRLR * (p1p2 * papb - pap2 * pbp1 + 5.*pap1 * pbp2) + 3.*m1 * m2 * papb * RRLR)) +
                           papb * (cc12 * (LL * LLRL * m1 * m2 * papb + 2.*LL * pap1 * pbp2 * RLRL +
                                           2.*LRLR * pap1 * pbp2 * RR + m1 * m2 * papb * RR * RRLR) +
                                   cc11 * (LL * (LLRL * m1 * m2 * papb + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL +
                                           pap1 * pbp2 * RLRL) + RR * (LRLR * p1p2 * papb - LRLR * pap2 * pbp1 +
                                                   LRLR * pap1 * pbp2 + m1 * m2 * papb * RRLR)) +
                                   cc23 * (LL * (2.*LLRL * m1 * m2 * papb + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL +
                                           3.*pap1 * pbp2 * RLRL) + RR * (LRLR * p1p2 * papb - LRLR * pap2 * pbp1 +
                                                   3.*LRLR * pap1 * pbp2 + 2.*m1 * m2 * papb * RRLR))));

    cc->SetIeps(ieps + 2);
    cc->GetNpf(la, 2);
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
    cc11 = la[0];
    cc12 = la[1];
    delete cc;

    complex<double> me2 = 16.*ivs3v1 * ivt1s2 * ((cc12 + cc23) * papb *
                          (LL * (LLRL * m1 * m2 * papb + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL) +
                           RR * (LRLR * p1p2 * papb - LRLR * pap2 * pbp1 + LRLR * pap1 * pbp2 +
                                 m1 * m2 * papb * RRLR)) +
                          cc24 * (LL * (3.*LLRL * m1 * m2 * papb + 2.*p1p2 * papb * RLRL - 2.*pap2 * pbp1 * RLRL +
                                        4.*pap1 * pbp2 * RLRL) + RR * (2.*LRLR * p1p2 * papb - 2.*LRLR * pap2 * pbp1 +
                                                4.*LRLR * pap1 * pbp2 + 3.*m1 * m2 * papb * RRLR)));

    return real(me0 + me1 + me2);
}

double FI::MVstr1u(int ieps)
{
    Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = -16.*ivs3v1 * ivu2s2 * (cc24 + (cc11 + cc23) * papb) *
                          (2.*LL * LLRL * pap2 * pbp1 + LL * m1 * m2 * papb * RLRL + LRLR * m1 * m2 * papb * RR +
                           2.*pap2 * pbp1 * RR * RRLR);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
    cc11 = la[0];
    cc12 = la[1];

    complex<double> me1 = 16.*ivs3v1 * ivu2s2 *
                          (cc24 * (LL * (LLRL * (p1p2 * papb + 5.*pap2 * pbp1 - pap1 * pbp2) +
                                         3.*m1 * m2 * papb * RLRL) + RR * (3.*LRLR * m1 * m2 * papb +
                                                 (p1p2 * papb + 5.*pap2 * pbp1 - pap1 * pbp2) * RRLR)) +
                           papb * (cc12 * (2.*LL * LLRL * pap2 * pbp1 + LL * m1 * m2 * papb * RLRL + LRLR * m1 * m2 * papb * RR +
                                           2.*pap2 * pbp1 * RR * RRLR) + cc23 * (LL * LLRL * p1p2 * papb + 3.*LL * LLRL * pap2 * pbp1 -
                                                   LL * LLRL * pap1 * pbp2 + 2.*LL * m1 * m2 * papb * RLRL + 2.*LRLR * m1 * m2 * papb * RR +
                                                   (p1p2 * papb + 3.*pap2 * pbp1 - pap1 * pbp2) * RR * RRLR) +
                                   cc11 * (LL * (LLRL * p1p2 * papb + LLRL * pap2 * pbp1 - LLRL * pap1 * pbp2 +
                                           m1 * m2 * papb * RLRL) + RR * (LRLR * m1 * m2 * papb + p1p2 * papb * RRLR +
                                                   pap2 * pbp1 * RRLR - pap1 * pbp2 * RRLR))));

    cc->SetIeps(ieps + 2);
    cc->GetNpf(la, 2);
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
    cc11 = la[0];
    cc12 = la[1];
    delete cc;

    complex<double> me2 = -16.*ivs3v1 * ivu2s2 *
                          (cc24 * (LL * (2.*LLRL * p1p2 * papb + 4.*LLRL * pap2 * pbp1 - 2.*LLRL * pap1 * pbp2 +
                                         3.*m1 * m2 * papb * RLRL) + RR * (3.*LRLR * m1 * m2 * papb + 2.*p1p2 * papb * RRLR +
                                                 4.*pap2 * pbp1 * RRLR - 2.*pap1 * pbp2 * RRLR)) +
                           (cc12 + cc23) * papb * (LL * (LLRL * p1p2 * papb + LLRL * pap2 * pbp1 - LLRL * pap1 * pbp2 +
                                   m1 * m2 * papb * RLRL) + RR * (LRLR * m1 * m2 * papb + p1p2 * papb * RRLR +
                                           pap2 * pbp1 * RRLR - pap1 * pbp2 * RRLR)));

    return real(me0 + me1 + me2);
}

double FI::MVstr2s(int ieps)
{
    Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc24 = la[3];

    complex<double> me0 = 16.*cc24 * ivs3v1 * ivs3v2 *
                          (RL * (m1 * m2 * papb * (LRRL + RRLL) + 2.*(LRLL * pap2 * pbp1 + pap1 * pbp2 * RRRL)) +
                           LR * (m1 * m2 * papb * (LRRR + RRLR) + 2.*(LRLR * pap1 * pbp2 + pap2 * pbp1 * RRRR)));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];

    complex<double> me1 = -16.*cc24 * ivs3v1 * ivs3v2 *
                          (RL * (LRRL * m1 * m2 * papb + LRLL * p1p2 * papb +
                                 3.*LRLL * pap2 * pbp1 - 3.*LRLL * pap1 * pbp2 + m1 * m2 * papb * RRLL + p1p2 * papb * RRRL -
                                 3.*pap2 * pbp1 * RRRL + 3.*pap1 * pbp2 * RRRL) +
                           LR * (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - 3.*LRLR * pap2 * pbp1 +
                                 3.*LRLR * pap1 * pbp2 + m1 * m2 * papb * RRLR + p1p2 * papb * RRRR + 3.*pap2 * pbp1 * RRRR -
                                 3.*pap1 * pbp2 * RRRR));

    cc->SetIeps(ieps + 2);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me2 = 32.*cc24 * ivs3v1 * ivs3v2 * (pap2 * pbp1 - pap1 * pbp2) *
                          (RL * (LRLL - RRRL) + LR * (-LRLR + RRRR));

    return real(me0 + me1 + me2);
}


double FI::MVstr2t(int ieps)
{
    Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    complex<double> me0 = 4.*ivs3v1 * ivt1s2 *
                          (2.*cc24 * (2.*LR * LRLR * pap1 * pbp2 + LRRL * m1 * m2 * papb * RL +
                                      LR * m1 * m2 * papb * RRLR + 2.*pap1 * pbp2 * RL * RRRL) +
                           ml1 * papb * (cc0 * LL * LRLL * m1 * pbp2 + 2.*cc11 * LL * LRLL * m1 * pbp2 -
                                         2.*cc12 * LL * LRLL * m1 * pbp2 + cc0 * LRRR * m2 * pap1 * RR + 2.*cc12 * LRRR * m2 * pap1 * RR +
                                         cc0 * LL * m2 * pap1 * RRLL + 2.*cc12 * LL * m2 * pap1 * RRLL +
                                         (cc0 + 2.*cc11 - 2.*cc12) * m1 * pbp2 * RR * RRRR));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 = -8.*cc24 * ivs3v1 * ivt1s2 *
                          (LR * (LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                                 m1 * m2 * papb * RRLR) + RL * (LRRL * m1 * m2 * papb +
                                         (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RRRL));

    return real(me0 + me1);
}


double FI::MVstr2u(int ieps)
{
    Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    complex<double> me0 = -4.*ivs3v1 * ivu2s2 *
                          (2.*cc24 * (LR * LRLR * m1 * m2 * papb + 2.*LRRL * pap2 * pbp1 * RL +
                                      2.*LR * pap2 * pbp1 * RRLR + m1 * m2 * papb * RL * RRRL) +
                           ml1 * papb * (cc0 * (LL * LRLL * m1 * pap2 + LRRR * m2 * pbp1 * RR + LL * m2 * pbp1 * RRLL +
                                         m1 * pap2 * RR * RRRR) + 2.*(cc11 * m2 * pbp1 * (LRRR * RR + LL * RRLL) +
                                                 cc12 * (LL * LRLL * m1 * pap2 - LRRR * m2 * pbp1 * RR - LL * m2 * pbp1 * RRLL +
                                                         m1 * pap2 * RR * RRRR))));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 = 8.*cc24 * ivs3v1 * ivu2s2 *
                          (LR * (LRLR * m1 * m2 * papb + (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RRLR) +
                           RL * (LRRL * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) + m1 * m2 * papb * RRRL));

    return real(me0 + me1);
}

double FI::MVstr3t(int ieps)
{
    Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];
    delete cc;

    complex<double> me0 = cc0 * ml1 * ivs3v1 * ivt1s2 *
                          (- 2.*pap2 * pbp1 * RR * RRRR - 2.*pap2 * pbp1 * LL * LRLL + 2.*pap1 * pbp2 * RR * RRRR +
                           2.*pap1 * pbp2 * LL * LRLL + 2.*papb * p1p2 * RR * RRRR + 2.*papb * p1p2 * LL * LRLL -
                           2.*m1 * m2 * papb * RR * LRRR - 2.*m1 * m2 * papb * LL * RRLL);

    return real(me0);
}

double FI::MVstr3u(int ieps)
{
    Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];
    delete cc;

    complex<double> me0 = cc0 * ml1 * ivs3v1 * ivu2s2 *
                          (2.*pap2 * pbp1 * RR * RRRR + 2.*pap2 * pbp1 * LL * LRLL - 2.*pap1 * pbp2 * RR * RRRR -
                           2.*pap1 * pbp2 * LL * LRLL + 2.*papb * p1p2 * RR * RRRR + 2.*papb * p1p2 * LL * LRLL -
                           2.*m1 * m2 * papb * RR * LRRR - 2.*m1 * m2 * papb * LL * RRLL);

    return real(me0);
}

double FI::MVttr1s(int ieps)
{
    Npf *cc = new Npf(0., m1s - 2.*pap1, m1s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = 4.*ivs3v2 * ivt1s1 * (4.*cc24 + (2.*cc12 + cc22) * m1s +
                          2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pap1) * (LR * m1 * m2 * papb * RLLL +
                                  2.*LR * pap1 * pbp2 * RLRL + LRRR * m1 * m2 * papb * RR + 2.*LRLR * pap1 * pbp2 * RR);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc22 = la[1];
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
    cc11 = la[0];
    cc12 = la[1];

    complex<double> me1 = -4.*ivs3v2 * ivt1s1 *
                          (((2.*cc12 + cc22) * m1s + 2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) *
                            pap1) * (LR * (m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL +
                                           pap1 * pbp2 * RLRL) + (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - LRLR * pap2 * pbp1 +
                                                   LRLR * pap1 * pbp2) * RR) +
                           cc24 * (LR * (6.*m1 * m2 * papb * RLLL + 4.*p1p2 * papb * RLRL - 4.*pap2 * pbp1 * RLRL +
                                         8.*pap1 * pbp2 * RLRL) + 2.*(3.*LRRR * m1 * m2 * papb + 2.*LRLR * p1p2 * papb -
                                                 2.*LRLR * pap2 * pbp1 + 4.*LRLR * pap1 * pbp2) * RR));

    cc->SetIeps(ieps + 2);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me2 = 8.*cc24 * ivs3v2 * ivt1s1 *
                          (LR * (m1 * m2 * papb * RLLL + (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RLRL) +
                           (LRRR * m1 * m2 * papb + LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) * RR);

    return real(me0 + me1 + me2);
}

double FI::MVttr1t(int ieps)
{
    Npf *cc = new Npf(0., m1s - 2.*pap1, m1s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = 4.*ivt1s1 * ivt1s2 * pap1 *
                          (4.*cc24 + (2.*cc12 + cc22) * m1s + 2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pap1) *
                          pbp2 * (LR * (LLLL + RLRL) + RR * (LRLR + RRRR));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 = -8.*cc24 * ivt1s1 * ivt1s2 * pap1 * pbp2 *
                          (LR * (LLLL + RLRL) + RR * (LRLR + RRRR));

    return real(me0 + me1);
}

double FI::MVttr1u(int ieps)
{
    Npf *cc = new Npf(0., m1s - 2.*pap1, m1s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = 2.*ivt1s1 * ivu2s2 * (4.*cc24 + (2.*cc12 + cc22) * m1s +
                          2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pap1) *
                          ((RRRR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL) * RR +
                           LR * (LRLR * m1 * m2 * papb + (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * LLLL));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 = -4.*cc24 * ivt1s1 * ivu2s2 *
                          ((RRRR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                            m1 * m2 * papb * RLRL) * RR + LR * (LRLR * m1 * m2 * papb +
                                    (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * LLLL));

    return -real(me0 + me1);
}

double FI::MVttr2s(int ieps)
{
    Npf *cc = new Npf(0., m1s - 2.*pap1, m2s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = 4.*ivs3v2 * ivt1s1 * (4.*cc24 + (2.*cc12 + cc22) * m2s +
                          2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pbp2) * (LR * m1 * m2 * papb * RLLL +
                                  2.*LR * pap1 * pbp2 * RLRL + LRRR * m1 * m2 * papb * RR + 2.*LRLR * pap1 * pbp2 * RR);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc22 = la[1];
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
    cc11 = la[0];
    cc12 = la[1];

    complex<double> me1 = -4.*ivs3v2 * ivt1s1 *
                          (((2.*cc12 + cc22) * m2s + 2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) *
                            pbp2) * (LR * (m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL +
                                           pap1 * pbp2 * RLRL) + (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - LRLR * pap2 * pbp1 +
                                                   LRLR * pap1 * pbp2) * RR) +
                           cc24 * (LR * (6.*m1 * m2 * papb * RLLL + 4.*p1p2 * papb * RLRL - 4.*pap2 * pbp1 * RLRL +
                                         8.*pap1 * pbp2 * RLRL) + 2.*(3.*LRRR * m1 * m2 * papb + 2.*LRLR * p1p2 * papb -
                                                 2.*LRLR * pap2 * pbp1 + 4.*LRLR * pap1 * pbp2) * RR));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me2 = 8.*cc24 * ivs3v2 * ivt1s1 *
                          (LR * (m1 * m2 * papb * RLLL + (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RLRL) +
                           (LRRR * m1 * m2 * papb + LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) * RR);

    return real(me0 + me1 + me2);
}

double FI::MVttr2t(int ieps)
{
    Npf *cc = new Npf(0., m1s - 2.*pap1, m2s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = 4.*ivt1s1 * ivt1s2 * pap1 * pbp2 *
                          (4.*cc24 + (2.*cc12 + cc22) * m2s + 2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pbp2)
                          * ((LLLL + LRLR) * RR + LR * (RLRL + RRRR));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 = -8.*cc24 * ivt1s1 * ivt1s2 * pap1 * pbp2 *
                          ((LLLL + LRLR) * RR + LR * (RLRL + RRRR));

    return real(me0 + me1);
}

double FI::MVttr2u(int ieps)
{
    Npf *cc = new Npf(0., m1s - 2.*pap1, m2s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = -2.*ivt1s1 * ivu2s2 * (4.*cc24 + (2.*cc12 + cc22) * m2s +
                          2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pbp2) *
                          ((LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) * RR +
                           LR * (m1 * m2 * papb * RLRL + (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 = 4.*cc24 * ivt1s1 * ivu2s2 *
                          ((LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) * RR +
                           LR * (m1 * m2 * papb * RLRL + (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR));

    return real(me0 + me1);
}

double FI::MVutr1s(int ieps)
{
    Npf *cc = new Npf(0., m2s - 2.*pap2, m2s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = -4.*ivs3v2 * ivu2s1 * (4.*cc24 + (2.*cc12 + cc22) * m2s +
                          2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pap2) * (2.*LR * pap2 * pbp1 * RLLL +
                                  LR * m1 * m2 * papb * RLRL + LRLR * m1 * m2 * papb * RR + 2.*LRRR * pap2 * pbp1 * RR);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc22 = la[1];
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
    cc11 = la[0];
    cc12 = la[1];

    complex<double> me1 = 4.*ivs3v2 * ivu2s1 *
                          (cc24 * (LR * (4.*p1p2 * papb * RLLL + 8.*pap2 * pbp1 * RLLL - 4.*pap1 * pbp2 * RLLL +
                                         6.*m1 * m2 * papb * RLRL) + 2.*(3.*LRLR * m1 * m2 * papb + 2.*LRRR * p1p2 * papb +
                                                 4.*LRRR * pap2 * pbp1 - 2.*LRRR * pap1 * pbp2) * RR) +
                           ((2.*cc12 + cc22) * m2s + 2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pap2) *
                           (LR * (p1p2 * papb * RLLL + pap2 * pbp1 * RLLL - pap1 * pbp2 * RLLL + m1 * m2 * papb * RLRL) +
                            (LRLR * m1 * m2 * papb + LRRR * p1p2 * papb + LRRR * pap2 * pbp1 - LRRR * pap1 * pbp2) * RR));

    cc->SetIeps(ieps + 2);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me2 = -8.*cc24 * ivs3v2 * ivu2s1 *
                          (LR * (p1p2 * papb * RLLL + pap2 * pbp1 * RLLL - pap1 * pbp2 * RLLL + m1 * m2 * papb * RLRL) +
                           (LRLR * m1 * m2 * papb + LRRR * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) * RR);

    return real(me0 + me1 + me2);
}

double FI::MVutr1t(int ieps)
{
    Npf *cc = new Npf(0., m2s - 2.*pap2, m2s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = -2.*ivt1s2 * ivu2s1 *
                          (4.*cc24 + (2.*cc12 + cc22) * m2s + 2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pap2) *
                          (LLLL * LR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + LR * m1 * m2 * papb * RLRL +
                           RR * (LRLR * m1 * m2 * papb - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 = 4.*cc24 * ivt1s2 * ivu2s1 *
                          (LLLL * LR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                           LR * m1 * m2 * papb * RLRL + RR * (LRLR * m1 * m2 * papb - p1p2 * papb * RRRR +
                                   pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR));

    return real(me0 + me1);
}

double FI::MVutr1u(int ieps)
{
    Npf *cc = new Npf(0., m2s - 2.*pap2, m2s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = 4.*ivu2s1 * ivu2s2 * pap2 *
                          (4.*cc24 + (2.*cc12 + cc22) * m2s + 2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pap2) *
                          pbp1 * (LR * (LLLL + RLRL) + RR * (LRLR + RRRR));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 = -8.*cc24 * ivu2s1 * ivu2s2 * pap2 * pbp1 *
                          (LR * (LLLL + RLRL) + RR * (LRLR + RRRR));

    return real(me0 + me1);
}

double FI::MVutr2s(int ieps)
{
    Npf *cc = new Npf(0., m2s - 2.*pap2, m1s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = -4.*ivs3v2 * ivu2s1 * (4.*cc24 + (2.*cc12 + cc22) * m1s +
                          2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pbp1) * (2.*LR * pap2 * pbp1 * RLLL +
                                  LR * m1 * m2 * papb * RLRL + LRLR * m1 * m2 * papb * RR + 2.*LRRR * pap2 * pbp1 * RR);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc22 = la[1];
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
    cc11 = la[0];
    cc12 = la[1];

    complex<double> me1 = 4.*ivs3v2 * ivu2s1 *
                          (cc24 * (LR * (4.*p1p2 * papb * RLLL + 8.*pap2 * pbp1 * RLLL - 4.*pap1 * pbp2 * RLLL +
                                         6.*m1 * m2 * papb * RLRL) + 2.*(3.*LRLR * m1 * m2 * papb + 2.*LRRR * p1p2 * papb +
                                                 4.*LRRR * pap2 * pbp1 - 2.*LRRR * pap1 * pbp2) * RR) +
                           ((2.*cc12 + cc22) * m1s + 2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pbp1) *
                           (LR * (p1p2 * papb * RLLL + pap2 * pbp1 * RLLL - pap1 * pbp2 * RLLL + m1 * m2 * papb * RLRL) +
                            (LRLR * m1 * m2 * papb + LRRR * p1p2 * papb + LRRR * pap2 * pbp1 - LRRR * pap1 * pbp2) * RR));

    cc->SetIeps(ieps + 2);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me2 = -8.*cc24 * ivs3v2 * ivu2s1 *
                          (LR * (p1p2 * papb * RLLL + pap2 * pbp1 * RLLL - pap1 * pbp2 * RLLL + m1 * m2 * papb * RLRL) +
                           (LRLR * m1 * m2 * papb + LRRR * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) * RR);

    return real(me0 + me1 + me2);
}

double FI::MVutr2t(int ieps)
{
    Npf *cc = new Npf(0., m2s - 2.*pap2, m1s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = -2.*ivt1s2 * ivu2s1 * (4.*cc24 + (2.*cc12 + cc22) * m1s +
                          2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pbp1) *
                          ((LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) * RR +
                           LR * (m1 * m2 * papb * RLRL + (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 = 4.*cc24 * ivt1s2 * ivu2s1 *
                          ((LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) * RR +
                           LR * (m1 * m2 * papb * RLRL + (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR));

    return real(me0 + me1);
}

double FI::MVutr2u(int ieps)
{
    Npf *cc = new Npf(0., m2s - 2.*pap2, m1s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = 4.*ivu2s1 * ivu2s2 * pap2 * pbp1 *
                          (4.*cc24 + (2.*cc12 + cc22) * m1s + 2.*(2.*cc11 - 2.*cc12 - cc22 + cc23) * pbp1) *
                          ((LLLL + LRLR) * RR + LR * (RLRL + RRRR));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 = -8.*cc24 * ivu2s1 * ivu2s2 * pap2 * pbp1 *
                          ((LLLL + LRLR) * RR + LR * (RLRL + RRRR));

    return real(me0 + me1);
}

double FI::MVttr3s(int ieps)
{
    Npf *cc = new Npf(0., m1s - 2.*pap1, m1s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
//   complex<double> cc11=la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    complex<double> me0 = -4.*ivs3v2 * ivt1s1 * (cc22 * LR * LRRR * m1 * m1s * m2 * papb -
                          2.*cc23 * LR * LRRR * m1 * m2 * pap1 * papb + 2.*cc22 * LR * LRLR * m1s * pap1 * pbp2 -
                          4.*cc23 * LR * LRLR * pow(pap1, 2) * pbp2 + cc0 * LL * m1 * m2 * ml1 * ml3 * papb * RLLL +
                          cc22 * m1 * m1s * m2 * papb * RL * RLLL - 2.*cc23 * m1 * m2 * pap1 * papb * RL * RLLL +
                          2.*cc0 * LL * ml1 * ml3 * pap1 * pbp2 * RLRL + 2.*cc22 * m1s * pap1 * pbp2 * RL * RLRL -
                          4.*cc23 * pow(pap1, 2) * pbp2 * RL * RLRL + 4.*cc24 * (LR * LRRR * m1 * m2 * papb +
                                  2.*LR * LRLR * pap1 * pbp2 + m1 * m2 * papb * RL * RLLL + 2.*pap1 * pbp2 * RL * RLRL) +
                          cc0 * LLRR * m1s * m2 * ml1 * papb * RR + cc0 * LRRR * m1 * m2 * ml1 * ml3 * papb * RR +
                          2.*cc0 * LLLR * m1 * ml1 * pap1 * pbp2 * RR + 2.*cc0 * LRLR * ml1 * ml3 * pap1 * pbp2 * RR +
                          cc0 * LL * m1s * m2 * ml1 * papb * RRLL + 2.*cc0 * LL * m1 * ml1 * pap1 * pbp2 * RRRL +
                          cc12 * (LR * m2 * (LLRR * m1s * ml3 + LRRR * m1 * (m1s - 2.*pap1)) * papb +
                                  2.*LR * (LLLR * m1 * ml3 + LRLR * (m1s - 2.*pap1)) * pap1 * pbp2 +
                                  2.*m1s * pap1 * pbp2 * RL * RLRL - 4.*pow(pap1, 2) * pbp2 * RL * RLRL +
                                  LLRR * m1s * m2 * ml1 * papb * RR + m1s * m2 * papb * (LL * ml1 + ml3 * RL) * RRLL +
                                  m1 * (m1s * m2 * papb * RL * RLLL + 2.*pap1 * (-(m2 * papb * RL * RLLL) + LLLR * ml1 * pbp2 * RR +
                                          LL * ml1 * pbp2 * RRRL + ml3 * pbp2 * RL * RRRL))));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc22 = la[1];
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
//   cc11=la[0];
    cc12 = la[1];
    cc->GetNpf(la, 0);
    cc0 = la[0];

    complex<double> me1 = 4.*ivs3v2 * ivt1s1 *
                          (cc22 * LR * LRRR * m1 * m1s * m2 * papb + cc22 * LR * LRLR * m1s * p1p2 * papb -
                           2.*cc23 * LR * LRRR * m1 * m2 * pap1 * papb - 2.*cc23 * LR * LRLR * p1p2 * pap1 * papb -
                           cc22 * LR * LRLR * m1s * pap2 * pbp1 + 2.*cc23 * LR * LRLR * pap1 * pap2 * pbp1 +
                           cc22 * LR * LRLR * m1s * pap1 * pbp2 - 2.*cc23 * LR * LRLR * pow(pap1, 2) * pbp2 +
                           cc0 * LL * m1 * m2 * ml1 * ml3 * papb * RLLL + cc22 * m1 * m1s * m2 * papb * RL * RLLL -
                           2.*cc23 * m1 * m2 * pap1 * papb * RL * RLLL + cc0 * LL * ml1 * ml3 * p1p2 * papb * RLRL -
                           cc0 * LL * ml1 * ml3 * pap2 * pbp1 * RLRL + cc0 * LL * ml1 * ml3 * pap1 * pbp2 * RLRL +
                           cc22 * m1s * p1p2 * papb * RL * RLRL - 2.*cc23 * p1p2 * pap1 * papb * RL * RLRL -
                           cc22 * m1s * pap2 * pbp1 * RL * RLRL + 2.*cc23 * pap1 * pap2 * pbp1 * RL * RLRL +
                           cc22 * m1s * pap1 * pbp2 * RL * RLRL - 2.*cc23 * pow(pap1, 2) * pbp2 * RL * RLRL +
                           cc24 * (LR * (6.*LRRR * m1 * m2 * papb + 4.*LRLR * (p1p2 * papb - pap2 * pbp1 +
                                         2.*pap1 * pbp2)) + 6.*m1 * m2 * papb * RL * RLLL +
                                   4.*(p1p2 * papb - pap2 * pbp1 + 2.*pap1 * pbp2) * RL * RLRL) +
                           cc0 * LLRR * m1s * m2 * ml1 * papb * RR + cc0 * LRRR * m1 * m2 * ml1 * ml3 * papb * RR +
                           cc0 * LLLR * m1 * ml1 * p1p2 * papb * RR + cc0 * LRLR * ml1 * ml3 * p1p2 * papb * RR -
                           cc0 * LLLR * m1 * ml1 * pap2 * pbp1 * RR - cc0 * LRLR * ml1 * ml3 * pap2 * pbp1 * RR +
                           cc0 * LLLR * m1 * ml1 * pap1 * pbp2 * RR + cc0 * LRLR * ml1 * ml3 * pap1 * pbp2 * RR +
                           cc0 * LL * m1s * m2 * ml1 * papb * RRLL + cc0 * LL * m1 * ml1 * (p1p2 * papb - pap2 * pbp1 +
                                   pap1 * pbp2) * RRRL +
                           cc12 * (LR * (LLRR * m1s * m2 * ml3 * papb + LRRR * m1 * m2 * (m1s - 2.*pap1) * papb +
                                         (LRLR * m1s + LLLR * m1 * ml3 - 2.*LRLR * pap1) * (p1p2 * papb - pap2 * pbp1 +
                                                 pap1 * pbp2)) + m1s * p1p2 * papb * RL * RLRL - 2.*p1p2 * pap1 * papb * RL * RLRL -
                                   m1s * pap2 * pbp1 * RL * RLRL + 2.*pap1 * pap2 * pbp1 * RL * RLRL +
                                   m1s * pap1 * pbp2 * RL * RLRL - 2.*pow(pap1, 2) * pbp2 * RL * RLRL + LLRR * m1s * m2 * ml1 * papb * RR +
                                   m1s * m2 * papb * (LL * ml1 + ml3 * RL) * RRLL +
                                   m1 * (m1s * m2 * papb * RL * RLLL - 2.*m2 * pap1 * papb * RL * RLLL +
                                         (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * (LLLR * ml1 * RR + LL * ml1 * RRRL +
                                                 ml3 * RL * RRRL))));

    cc->SetIeps(ieps + 2);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me2 = -8.*cc24 * ivs3v2 * ivt1s1 *
                          (LR * (LRRR * m1 * m2 * papb + LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) +
                           RL * (m1 * m2 * papb * RLLL + (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RLRL));

    return real(me0 + me1 + me2);
}

double FI::MVttr3t(int ieps)
{
    Npf *cc = new Npf(0., m1s - 2.*pap1, m1s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
//   complex<double> cc11=la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    complex<double> me0 = -4.*ivt1s1 * ivt1s2 * pap1 * pbp2 *
                          (cc22 * LR * LRLR * m1s - 2.*cc23 * LR * LRLR * pap1 +
                           cc22 * LLLL * m1s * RL - 2.*cc23 * LLLL * pap1 * RL + cc22 * m1s * RL * RLRL -
                           2.*cc23 * pap1 * RL * RLRL + cc0 * LRLR * ml1 * ml3 * RR + cc0 * m1 * ml1 * RLRR * RR +
                           cc0 * ml1 * (LL * (LRLL * m1 + ml3 * (LLLL + RLRL)) + LLLR * m1 * RR) +
                           cc0 * LL * m1 * ml1 * RRRL + cc22 * LR * m1s * RRRR - 2.*cc23 * LR * pap1 * RRRR +
                           cc0 * ml1 * ml3 * RR * RRRR + 4.*cc24 * (RL * (LLLL + RLRL) + LR * (LRLR + RRRR)) +
                           cc12 * (LLLL * m1s * RL + LRLL * m1 * ml3 * RL - 2.*LLLL * pap1 * RL + m1s * RL * RLRL -
                                   2.*pap1 * RL * RLRL + LLLR * m1 * ml1 * RR + m1 * ml1 * RLRR * RR + m1 * ml3 * RL * RRRL +
                                   LL * m1 * ml1 * (LRLL + RRRL) + LR * (m1 * ml3 * (LLLR + RLRR) +
                                           (m1s - 2.*pap1) * (LRLR + RRRR))));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 = 8.*cc24 * ivt1s1 * ivt1s2 *
                          pap1 * pbp2 * (RL * (LLLL + RLRL) + LR * (LRLR + RRRR));

    return real(me0 + me1);
}

double FI::MVttr3u(int ieps)
{
    Npf *cc = new Npf(0., m1s - 2.*pap1, m1s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
//   complex<double> cc11=la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    complex<double> me0 = 2.*ivt1s1 * ivu2s2 *
                          (4.*cc24 * LR * LRLR * m1 * m2 * papb + cc22 * LR * LRLR * m1 * m1s * m2 * papb -
                           cc0 * LL * LRLL * m1 * ml1 * p1p2 * papb - cc12 * LL * LRLL * m1 * ml1 * p1p2 * papb -
                           cc0 * LL * LLLL * ml1 * ml3 * p1p2 * papb + cc12 * LR * LRLR * m1 * m2 * (m1s - 2.*pap1) * papb -
                           2.*cc23 * LR * LRLR * m1 * m2 * pap1 * papb + cc0 * LL * LRLL * m1 * ml1 * pap2 * pbp1 +
                           cc12 * LL * LRLL * m1 * ml1 * pap2 * pbp1 + cc0 * LL * LLLL * ml1 * ml3 * pap2 * pbp1 +
                           cc0 * LL * LRLL * m1 * ml1 * pap1 * pbp2 + cc12 * LL * LRLL * m1 * ml1 * pap1 * pbp2 +
                           cc0 * LL * LLLL * ml1 * ml3 * pap1 * pbp2 - cc12 * LLLL * m1s * p1p2 * papb * RL -
                           cc22 * LLLL * m1s * p1p2 * papb * RL - cc12 * LRLL * m1 * ml3 * p1p2 * papb * RL +
                           2.*cc12 * LLLL * p1p2 * pap1 * papb * RL + 2.*cc23 * LLLL * p1p2 * pap1 * papb * RL +
                           cc12 * LLLL * m1s * pap2 * pbp1 * RL + cc22 * LLLL * m1s * pap2 * pbp1 * RL +
                           cc12 * LRLL * m1 * ml3 * pap2 * pbp1 * RL - 2.*cc12 * LLLL * pap1 * pap2 * pbp1 * RL -
                           2.*cc23 * LLLL * pap1 * pap2 * pbp1 * RL + cc12 * LLLL * m1s * pap1 * pbp2 * RL +
                           cc22 * LLLL * m1s * pap1 * pbp2 * RL + cc12 * LRLL * m1 * ml3 * pap1 * pbp2 * RL -
                           2.*cc12 * LLLL * pow(pap1, 2) * pbp2 * RL -
                           2.*cc23 * LLLL * pow(pap1, 2) * pbp2 * RL +
                           cc0 * LL * m1 * m2 * ml1 * ml3 * papb * RLRL + cc12 * m1 * m1s * m2 * papb * RL * RLRL +
                           cc22 * m1 * m1s * m2 * papb * RL * RLRL - 2.*cc12 * m1 * m2 * pap1 * papb * RL * RLRL -
                           2.*cc23 * m1 * m2 * pap1 * papb * RL * RLRL +
                           4.*cc24 * RL * (LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL) +
                           cc12 * LR * ml3 * (LLLR * m1s * m2 * papb + m1 * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) *
                                   RLRR) + cc0 * LLLR * m1s * m2 * ml1 * papb * RR + cc12 * LLLR * m1s * m2 * ml1 * papb * RR +
                           cc0 * LRLR * m1 * m2 * ml1 * ml3 * papb * RR - cc0 * m1 * ml1 * p1p2 * papb * RLRR * RR -
                           cc12 * m1 * ml1 * p1p2 * papb * RLRR * RR + cc0 * m1 * ml1 * pap2 * pbp1 * RLRR * RR +
                           cc12 * m1 * ml1 * pap2 * pbp1 * RLRR * RR + cc0 * m1 * ml1 * pap1 * pbp2 * RLRR * RR +
                           cc12 * m1 * ml1 * pap1 * pbp2 * RLRR * RR + cc0 * LL * m1s * m2 * ml1 * papb * RRRL +
                           cc12 * m1s * m2 * papb * (LL * ml1 + ml3 * RL) * RRRL +
                           (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) * (-4.*cc24 * LR - (cc12 + cc22) * LR * m1s +
                                   2.*(cc12 + cc23) * LR * pap1 - cc0 * ml1 * ml3 * RR) * RRRR);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 = -4.*cc24 * ivt1s1 * ivu2s2 *
                          (RL * (LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                                 m1 * m2 * papb * RLRL) + LR * (LRLR * m1 * m2 * papb +
                                         (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR));

    return real(me0 + me1);
}

double FI::MVttr4s(int ieps)
{
    Npf *cc = new Npf(0., m1s - 2.*pap1, m2s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
//   complex<double> cc11=la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    complex<double> me0 = -4.*ivs3v2 * ivt1s1 *
                          ((4.*cc24 + cc22 * m2s - 2.*cc23 * pbp2) *
                           (LRRR * m1 * m2 * papb * RL + LR * m1 * m2 * papb * RLLL +
                            2.*pap1 * pbp2 * (LRLR * RL + LR * RLRL)) +
                           cc12 * (LRRR * m1 * m2 * m2s * papb * RL + 2.*LRLR * m2s * pap1 * pbp2 * RL -
                                   2.*LRRR * m1 * m2 * papb * pbp2 * RL - 4.*LRLR * pap1 * pow(pbp2, 2) * RL +
                                   LR * m1 * m2 * m2s * papb * RLLL - 2.*LR * m1 * m2 * papb * pbp2 * RLLL +
                                   2.*LR * m2s * pap1 * pbp2 * RLRL - 4.*LR * pap1 * pow(pbp2, 2) * RLRL +
                                   LLLL * m1 * m2s * papb * (LR * ml3 + ml1 * RR) + 2.*LLRL * m2 * pap1 * pbp2 *
                                   (LR * ml3 + ml1 * RR) + 2.*LL * m2 * ml1 * pap1 * pbp2 * RRLR +
                                   2.*m2 * ml3 * pap1 * pbp2 * RL * RRLR + m1 * m2s * papb * (LL * ml1 + ml3 * RL) * RRRR) +
                           cc0 * ml1 * ((LLLL * m1 * m2s * papb + 2.*LLRL * m2 * pap1 * pbp2 + m1 * m2 * ml3 * papb * RLLL +
                                         2.*ml3 * pap1 * pbp2 * RLRL) * RR + LL * (LRRR * m1 * m2 * ml3 * papb +
                                                 2.*pap1 * pbp2 * (LRLR * ml3 + m2 * RRLR) + m1 * m2s * papb * RRRR)));

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc22 = la[1];
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
//   cc11=la[0];
    cc12 = la[1];
    cc->GetNpf(la, 0);
    cc0 = la[0];

    complex<double> me1 = 4.*ivs3v2 * ivt1s1 *
                          (.6 * cc24 * LRRR * m1 * m2 * papb * RL + cc22 * LRRR * m1 * m2 * m2s * papb * RL +
                           4.*cc24 * LRLR * p1p2 * papb * RL + cc22 * LRLR * m2s * p1p2 * papb * RL -
                           4.*cc24 * LRLR * pap2 * pbp1 * RL - cc22 * LRLR * m2s * pap2 * pbp1 * RL +
                           8.*cc24 * LRLR * pap1 * pbp2 * RL + cc22 * LRLR * m2s * pap1 * pbp2 * RL -
                           2.*cc23 * LRRR * m1 * m2 * papb * pbp2 * RL - 2.*cc23 * LRLR * p1p2 * papb * pbp2 * RL +
                           2.*cc23 * LRLR * pap2 * pbp1 * pbp2 * RL - 2.*cc23 * LRLR * pap1 * pow(pbp2, 2) * RL +
                           6.*cc24 * LR * m1 * m2 * papb * RLLL + cc22 * LR * m1 * m2 * m2s * papb * RLLL -
                           2.*cc23 * LR * m1 * m2 * papb * pbp2 * RLLL +
                           LR * ((cc22 * m2s - 2.*cc23 * pbp2) * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                                 4.*cc24 * (p1p2 * papb - pap2 * pbp1 + 2.*pap1 * pbp2)) * RLRL +
                           cc12 * (LRRR * m1 * m2 * m2s * papb * RL + LRLR * m2s * p1p2 * papb * RL -
                                   LRLR * m2s * pap2 * pbp1 * RL + LRLR * m2s * pap1 * pbp2 * RL -
                                   2.*LRRR * m1 * m2 * papb * pbp2 * RL - 2.*LRLR * p1p2 * papb * pbp2 * RL +
                                   2.*LRLR * pap2 * pbp1 * pbp2 * RL - 2.*LRLR * pap1 * pow(pbp2, 2) * RL +
                                   LR * m1 * m2 * m2s * papb * RLLL - 2.*LR * m1 * m2 * papb * pbp2 * RLLL +
                                   LR * m2s * p1p2 * papb * RLRL - LR * m2s * pap2 * pbp1 * RLRL + LR * m2s * pap1 * pbp2 * RLRL -
                                   2.*LR * p1p2 * papb * pbp2 * RLRL + 2.*LR * pap2 * pbp1 * pbp2 * RLRL -
                                   2.*LR * pap1 * pow(pbp2, 2) * RLRL + LLLL * m1 * m2s * papb * (LR * ml3 + ml1 * RR) +
                                   LLRL * m2 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * (LR * ml3 + ml1 * RR) +
                                   LL * m2 * ml1 * p1p2 * papb * RRLR - LL * m2 * ml1 * pap2 * pbp1 * RRLR +
                                   LL * m2 * ml1 * pap1 * pbp2 * RRLR + m2 * ml3 * p1p2 * papb * RL * RRLR -
                                   m2 * ml3 * pap2 * pbp1 * RL * RRLR + m2 * ml3 * pap1 * pbp2 * RL * RRLR +
                                   m1 * m2s * papb * (LL * ml1 + ml3 * RL) * RRRR) +
                           cc0 * ml1 * ((LLLL * m1 * m2s * papb + m2 * (LLRL * p1p2 * papb - LLRL * pap2 * pbp1 +
                                         LLRL * pap1 * pbp2 + m1 * ml3 * papb * RLLL) +
                                         ml3 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RLRL) * RR +
                                        LL * (LRRR * m1 * m2 * ml3 * papb + (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) *
                                                (LRLR * ml3 + m2 * RRLR) + m1 * m2s * papb * RRRR)));

    cc->SetIeps(ieps + 2);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me2 = -8.*cc24 * ivs3v2 * ivt1s1 *
                          (LRRR * m1 * m2 * papb * RL + LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RL +
                           LR * (m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL));

    return real(me0 + me1 + me2);
}

double FI::MVttr4t(int ieps)
{
    Npf *cc = new Npf(0., m1s - 2.*pap1, m2s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
//   complex<double> cc11=la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    complex<double> me0 =
        + cc0 * ivt1s1 * ivt1s2 * (- 4.*pap1 * pbp2 * ml1 * ml3 * RR * RLRL - 4.*pap1 *
                                   pbp2 * ml1 * ml3 * RR * RRRR - 4.*pap1 * pbp2 * ml1 * ml3 * LL * LRLR - 4.*pap1 *
                                   pbp2 * ml1 * ml3 * LL * LLLL - 4.*m2 * pap1 * pbp2 * ml1 * RR * LRRR - 4.*m2 * pap1 *
                                   pbp2 * ml1 * RR * LLRL - 4.*m2 * pap1 * pbp2 * ml1 * LL * RRLR - 4.*m2 * pap1 * pbp2
                                   * ml1 * LL * RLLL);
    me0 +=  + cc12 * ivt1s1 * ivt1s2 * (- 4.*pap1 * pbp2 * RL * LRLR * m2s - 4.*
                                        pap1 * pbp2 * RL * LLLL * m2s - 4.*pap1 * pbp2 * LR * RLRL * m2s - 4.*pap1 * pbp2 *
                                        LR * RRRR * m2s + 8.*pap1 * pow(pbp2, 2) * RL * LRLR + 8.*pap1 * pow(pbp2, 2) *
                                        RL * LLLL + 8.*pap1 * pow(pbp2, 2) * LR * RLRL + 8.*pap1 * pow(pbp2, 2) * LR *
                                        RRRR - 4.*m2 * pap1 * pbp2 * ml3 * RL * RRLR - 4.*m2 * pap1 * pbp2 * ml3 * RL * RLLL
                                        - 4.*m2 * pap1 * pbp2 * ml3 * LR * LRRR - 4.*m2 * pap1 * pbp2 * ml3 * LR * LLRL - 4.
                                        *m2 * pap1 * pbp2 * ml1 * RR * LRRR - 4.*m2 * pap1 * pbp2 * ml1 * RR * LLRL - 4.*m2 *
                                        pap1 * pbp2 * ml1 * LL * RRLR - 4.*m2 * pap1 * pbp2 * ml1 * LL * RLLL);
    me0 +=  + cc22 * ivt1s1 * ivt1s2 * (- 4.*pap1 * pbp2 * RL * LRLR * m2s - 4.*
                                        pap1 * pbp2 * RL * LLLL * m2s - 4.*pap1 * pbp2 * LR * RLRL * m2s - 4.*pap1 * pbp2 *
                                        LR * RRRR * m2s);
    me0 +=  + cc23 * ivt1s1 * ivt1s2 * (8.*pap1 * pow(pbp2, 2) * RL * LRLR + 8.*
                                        pap1 * pow(pbp2, 2) * RL * LLLL + 8.*pap1 * pow(pbp2, 2) * LR * RLRL + 8.*pap1
                                        * pow(pbp2, 2) * LR * RRRR);
    me0 +=  + cc24 * ivt1s1 * ivt1s2 * (- 16.*pap1 * pbp2 * RL * LRLR - 16.*
                                        pap1 * pbp2 * RL * LLLL - 16.*pap1 * pbp2 * LR * RLRL - 16.*pap1 * pbp2 * LR *
                                        RRRR);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 =
        + cc24 * ivt1s1 * ivt1s2 * (8.*pap1 * pbp2 * RL * LRLR + 8.*pap1 * pbp2 * RL *
                                    LLLL + 8.*pap1 * pbp2 * LR * RLRL + 8.*pap1 * pbp2 * LR * RRRR);

    return real(me0 + me1);
}

double FI::MVttr4u(int ieps)
{
    Npf *cc = new Npf(0., m1s - 2.*pap1, m2s, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
//   complex<double> cc11=la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    complex<double> me0 =
        + cc0 * ivt1s1 * ivu2s2 * (2.*pap2 * pbp1 * ml1 * ml3 * RR * RRRR + 2.*pap2 *
                                   pbp1 * ml1 * ml3 * LL * LLLL + 2.*pap1 * pbp2 * ml1 * ml3 * RR * RRRR + 2.*pap1 *
                                   pbp2 * ml1 * ml3 * LL * LLLL - 2.*papb * p1p2 * ml1 * ml3 * RR * RRRR - 2.*papb *
                                   p1p2 * ml1 * ml3 * LL * LLLL + 2.*m2 * pap2 * pbp1 * ml1 * RR * LRRR + 2.*m2 * pap2 *
                                   pbp1 * ml1 * LL * RLLL + 2.*m2 * pap1 * pbp2 * ml1 * RR * LRRR + 2.*m2 * pap1 * pbp2
                                   * ml1 * LL * RLLL - 2.*m2 * papb * p1p2 * ml1 * RR * LRRR - 2.*m2 * papb * p1p2 * ml1
                                   * LL * RLLL + 2.*m1 * papb * ml1 * RR * LLRL * m2s + 2.*m1 * papb * ml1 * LL * RRLR *
                                   m2s + 2.*m1 * m2 * papb * ml1 * ml3 * RR * RLRL + 2.*m1 * m2 * papb * ml1 * ml3 * LL *
                                   LRLR);
    me0 +=  + cc12 * ivt1s1 * ivu2s2 * (2.*pap2 * pbp1 * RL * LLLL * m2s + 2.*pap2
                                        * pbp1 * LR * RRRR * m2s - 4.*pap2 * pbp1 * pbp2 * RL * LLLL - 4.*pap2 * pbp1 *
                                        pbp2 * LR * RRRR + 2.*pap1 * pbp2 * RL * LLLL * m2s + 2.*pap1 * pbp2 * LR * RRRR *
                                        m2s - 4.*pap1 * pow(pbp2, 2) * RL * LLLL - 4.*pap1 * pow(pbp2, 2) * LR * RRRR
                                        - 2.*papb * p1p2 * RL * LLLL * m2s - 2.*papb * p1p2 * LR * RRRR * m2s + 4.*papb *
                                        pbp2 * p1p2 * RL * LLLL + 4.*papb * pbp2 * p1p2 * LR * RRRR + 2.*m2 * pap2 * pbp1 *
                                        ml3 * RL * RLLL + 2.*m2 * pap2 * pbp1 * ml3 * LR * LRRR + 2.*m2 * pap2 * pbp1 * ml1 *
                                        RR * LRRR + 2.*m2 * pap2 * pbp1 * ml1 * LL * RLLL + 2.*m2 * pap1 * pbp2 * ml3 * RL *
                                        RLLL + 2.*m2 * pap1 * pbp2 * ml3 * LR * LRRR + 2.*m2 * pap1 * pbp2 * ml1 * RR * LRRR
                                        + 2.*m2 * pap1 * pbp2 * ml1 * LL * RLLL - 2.*m2 * papb * p1p2 * ml3 * RL * RLLL - 2
                                        * m2 * papb * p1p2 * ml3 * LR * LRRR - 2.*m2 * papb * p1p2 * ml1 * RR * LRRR - 2.*m2 *
                                        papb * p1p2 * ml1 * LL * RLLL + 2.*m1 * papb * ml3 * RL * RRLR * m2s + 2.*m1 * papb *
                                        ml3 * LR * LLRL * m2s + 2.*m1 * papb * ml1 * RR * LLRL * m2s + 2.*m1 * papb * ml1 * LL
                                        * RRLR * m2s + 2.*m1 * m2 * papb * RL * LRLR * m2s + 2.*m1 * m2 * papb * LR * RLRL *
                                        m2s - 4.*m1 * m2 * papb * pbp2 * RL * LRLR - 4.*m1 * m2 * papb * pbp2 * LR * RLRL);
    me0 +=  + cc22 * ivt1s1 * ivu2s2 * (2.*pap2 * pbp1 * RL * LLLL * m2s + 2.*pap2
                                        * pbp1 * LR * RRRR * m2s + 2.*pap1 * pbp2 * RL * LLLL * m2s + 2.*pap1 * pbp2 * LR *
                                        RRRR * m2s - 2.*papb * p1p2 * RL * LLLL * m2s - 2.*papb * p1p2 * LR * RRRR * m2s
                                        + 2.*m1 * m2 * papb * RL * LRLR * m2s + 2.*m1 * m2 * papb * LR * RLRL * m2s);
    me0 +=  + cc23 * ivt1s1 * ivu2s2 * (- 4.*pap2 * pbp1 * pbp2 * RL * LLLL - 4.*
                                        pap2 * pbp1 * pbp2 * LR * RRRR - 4.*pap1 * pow(pbp2, 2) * RL * LLLL - 4.*pap1 *
                                        pow(pbp2, 2) * LR * RRRR + 4.*papb * pbp2 * p1p2 * RL * LLLL + 4.*papb * pbp2 *
                                        p1p2 * LR * RRRR - 4.*m1 * m2 * papb * pbp2 * RL * LRLR - 4.*m1 * m2 * papb * pbp2 *
                                        LR * RLRL);
    me0 +=  + cc24 * ivt1s1 * ivu2s2 * (8.*pap2 * pbp1 * RL * LLLL + 8.*pap2 *
                                        pbp1 * LR * RRRR + 8.*pap1 * pbp2 * RL * LLLL + 8.*pap1 * pbp2 * LR * RRRR - 8.*
                                        papb * p1p2 * RL * LLLL - 8.*papb * p1p2 * LR * RRRR + 8.*m1 * m2 * papb * RL * LRLR
                                        + 8.*m1 * m2 * papb * LR * RLRL);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    complex<double> me1 =
        + cc24 * ivt1s1 * ivu2s2 * (- 4.*pap2 * pbp1 * RL * LLLL - 4.*pap2 * pbp1 * LR
                                    * RRRR - 4.*pap1 * pbp2 * RL * LLLL - 4.*pap1 * pbp2 * LR * RRRR + 4.*papb *
                                    p1p2 * RL * LLLL + 4.*papb * p1p2 * LR * RRRR - 4.*m1 * m2 * papb * RL * LRLR - 4.*
                                    m1 * m2 * papb * LR * RLRL);

    return real(me0 + me1);
}

double FI::MVutr3s(int ieps)
{
    complex<double> me0 = complex<double>(0., 0.);
    complex<double> me1 = complex<double>(0., 0.);
    complex<double> me2 = complex<double>(0., 0.);

    complex<double> la[7];
    Npf *cc = new Npf(0., m2s - 2.*pap2, m2s, m1l, m2l, m3l, mul, ieps);

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
//   complex<double> cc11=la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    me0 =  + ivs3v2 * ivu2s1 * (8.*cc0 * pap2 * pbp1 * ml1 * ml3 * RR * LRRR + 8.*cc0 * pap2
                                * pbp1 * ml1 * ml3 * LL * RLLL + 8.*cc0 * m2 * pap2 * pbp1 * ml1 * RR * LLRR + 8.*cc0
                                * m2 * pap2 * pbp1 * ml1 * LL * RRLL + 4.*cc0 * m1 * m2 * papb * ml1 * ml3 * RR * LRLR
                                + 4.*cc0 * m1 * m2 * papb * ml1 * ml3 * LL * RLRL + 4.*cc0 * m1 * m2s * papb *
                                ml1 * RR * LLLR + 4.*cc0 * m1 * m2s * papb * ml1 * LL * RRRL - 16.*cc12 *
                                pow(pap2, 2) * pbp1 * RL * RLLL - 16.*cc12 * pow(pap2, 2) * pbp1 * LR * LRRR +
                                8.*cc12 * m2 * pap2 * pbp1 * ml3 * RL * RRLL + 8.*cc12 * m2 * pap2 * pbp1 * ml3 * LR *
                                LLRR + 8.*cc12 * m2 * pap2 * pbp1 * ml1 * RR * LLRR + 8.*cc12 * m2 * pap2 * pbp1 *
                                ml1 * LL * RRLL + 8.*cc12 * m2s * pap2 * pbp1 * RL * RLLL + 8.*cc12 * pow(
                                    m2, 2) * pap2 * pbp1 * LR * LRRR - 8.*cc12 * m1 * m2 * papb * pap2 * RL * RLRL - 8.*
                                cc12 * m1 * m2 * papb * pap2 * LR * LRLR + 4.*cc12 * m1 * m2s * papb * ml3 * RL
                                * RRRL + 4.*cc12 * m1 * m2s * papb * ml3 * LR * LLLR + 4.*cc12 * m1 * pow(
                                    m2, 2) * papb * ml1 * RR * LLLR + 4.*cc12 * m1 * m2s * papb * ml1 * LL * RRRL
                                + 4.*cc12 * m1 * pow(m2, 3) * papb * RL * RLRL + 4.*cc12 * m1 * pow(m2, 3) * papb
                                * LR * LRLR + 8.*cc22 * m2s * pap2 * pbp1 * RL * RLLL + 8.*cc22 * pow(
                                    m2, 2) * pap2 * pbp1 * LR * LRRR);
    me0 +=  + ivs3v2 * ivu2s1 * (4.*cc22 * m1 * pow(m2, 3) * papb * RL * RLRL + 4.*
                                 cc22 * m1 * pow(m2, 3) * papb * LR * LRLR - 16.*cc23 * pow(pap2, 2) * pbp1 * RL *
                                 RLLL - 16.*cc23 * pow(pap2, 2) * pbp1 * LR * LRRR - 8.*cc23 * m1 * m2 * papb *
                                 pap2 * RL * RLRL - 8.*cc23 * m1 * m2 * papb * pap2 * LR * LRLR + 32.*cc24 * pap2 *
                                 pbp1 * RL * RLLL + 32.*cc24 * pap2 * pbp1 * LR * LRRR + 16.*cc24 * m1 * m2 * papb *
                                 RL * RLRL + 16.*cc24 * m1 * m2 * papb * LR * LRLR);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc22 = la[1];
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
//   cc11=la[0];
    cc12 = la[1];
    cc->GetNpf(la, 0);
    cc0 = la[0];

    me1 =  + ivs3v2 * ivu2s1 * (- 4.*cc0 * pap2 * pbp1 * ml1 * ml3 * RR * LRRR - 4.*cc0 *
                                pap2 * pbp1 * ml1 * ml3 * LL * RLLL + 4.*cc0 * pap1 * pbp2 * ml1 * ml3 * RR * LRRR +
                                4.*cc0 * pap1 * pbp2 * ml1 * ml3 * LL * RLLL - 4.*cc0 * papb * p1p2 * ml1 * ml3 * RR *
                                LRRR - 4.*cc0 * papb * p1p2 * ml1 * ml3 * LL * RLLL - 4.*cc0 * m2 * pap2 * pbp1 *
                                ml1 * RR * LLRR - 4.*cc0 * m2 * pap2 * pbp1 * ml1 * LL * RRLL + 4.*cc0 * m2 * pap1 *
                                pbp2 * ml1 * RR * LLRR + 4.*cc0 * m2 * pap1 * pbp2 * ml1 * LL * RRLL - 4.*cc0 * m2 *
                                papb * p1p2 * ml1 * RR * LLRR - 4.*cc0 * m2 * papb * p1p2 * ml1 * LL * RRLL - 4.*cc0
                                * m1 * m2 * papb * ml1 * ml3 * RR * LRLR - 4.*cc0 * m1 * m2 * papb * ml1 * ml3 * LL * RLRL
                                - 4.*cc0 * m1 * m2s * papb * ml1 * RR * LLLR - 4.*cc0 * m1 * m2s *
                                papb * ml1 * LL * RRRL + 8.*cc12 * pow(pap2, 2) * pbp1 * RL * RLLL + 8.*cc12 *
                                pow(pap2, 2) * pbp1 * LR * LRRR - 8.*cc12 * pap1 * pap2 * pbp2 * RL * RLLL - 8.*
                                cc12 * pap1 * pap2 * pbp2 * LR * LRRR + 8.*cc12 * papb * pap2 * p1p2 * RL * RLLL +
                                8.*cc12 * papb * pap2 * p1p2 * LR * LRRR - 4.*cc12 * m2 * pap2 * pbp1 * ml3 * RL *
                                RRLL - 4.*cc12 * m2 * pap2 * pbp1 * ml3 * LR * LLRR - 4.*cc12 * m2 * pap2 * pbp1 *
                                ml1 * RR * LLRR - 4.*cc12 * m2 * pap2 * pbp1 * ml1 * LL * RRLL + 4.*cc12 * m2 * pap1
                                * pbp2 * ml3 * RL * RRLL);
    me1 +=  + ivs3v2 * ivu2s1 * (4.*cc12 * m2 * pap1 * pbp2 * ml3 * LR * LLRR + 4.*
                                 cc12 * m2 * pap1 * pbp2 * ml1 * RR * LLRR + 4.*cc12 * m2 * pap1 * pbp2 * ml1 * LL *
                                 RRLL - 4.*cc12 * m2 * papb * p1p2 * ml3 * RL * RRLL - 4.*cc12 * m2 * papb * p1p2 *
                                 ml3 * LR * LLRR - 4.*cc12 * m2 * papb * p1p2 * ml1 * RR * LLRR - 4.*cc12 * m2 * papb
                                 * p1p2 * ml1 * LL * RRLL - 4.*cc12 * m2s * pap2 * pbp1 * RL * RLLL - 4.*
                                 cc12 * m2s * pap2 * pbp1 * LR * LRRR + 4.*cc12 * m2s * pap1 * pbp2 *
                                 RL * RLLL + 4.*cc12 * m2s * pap1 * pbp2 * LR * LRRR - 4.*cc12 * pow(
                                     m2, 2) * papb * p1p2 * RL * RLLL - 4.*cc12 * m2s * papb * p1p2 * LR * LRRR
                                 + 8.*cc12 * m1 * m2 * papb * pap2 * RL * RLRL + 8.*cc12 * m1 * m2 * papb * pap2 * LR *
                                 LRLR - 4.*cc12 * m1 * m2s * papb * ml3 * RL * RRRL - 4.*cc12 * m1 * pow(
                                     m2, 2) * papb * ml3 * LR * LLLR - 4.*cc12 * m1 * m2s * papb * ml1 * RR * LLLR
                                 - 4.*cc12 * m1 * m2s * papb * ml1 * LL * RRRL - 4.*cc12 * m1 * pow(m2, 3) *
                                 papb * RL * RLRL - 4.*cc12 * m1 * pow(m2, 3) * papb * LR * LRLR - 4.*cc22 * pow(
                                     m2, 2) * pap2 * pbp1 * RL * RLLL - 4.*cc22 * m2s * pap2 * pbp1 * LR * LRRR
                                 + 4.*cc22 * m2s * pap1 * pbp2 * RL * RLLL + 4.*cc22 * m2s * pap1 *
                                 pbp2 * LR * LRRR);
    me1 +=  + ivs3v2 * ivu2s1 * (- 4.*cc22 * m2s * papb * p1p2 * RL * RLLL
                                 - 4.*cc22 * m2s * papb * p1p2 * LR * LRRR - 4.*cc22 * m1 * pow(m2, 3) *
                                 papb * RL * RLRL - 4.*cc22 * m1 * pow(m2, 3) * papb * LR * LRLR + 8.*cc23 * pow(
                                     pap2, 2) * pbp1 * RL * RLLL + 8.*cc23 * pow(pap2, 2) * pbp1 * LR * LRRR - 8.*
                                 cc23 * pap1 * pap2 * pbp2 * RL * RLLL - 8.*cc23 * pap1 * pap2 * pbp2 * LR * LRRR +
                                 8.*cc23 * papb * pap2 * p1p2 * RL * RLLL + 8.*cc23 * papb * pap2 * p1p2 * LR * LRRR
                                 + 8.*cc23 * m1 * m2 * papb * pap2 * RL * RLRL + 8.*cc23 * m1 * m2 * papb * pap2 * LR *
                                 LRLR - 32.*cc24 * pap2 * pbp1 * RL * RLLL - 32.*cc24 * pap2 * pbp1 * LR * LRRR
                                 + 16.*cc24 * pap1 * pbp2 * RL * RLLL + 16.*cc24 * pap1 * pbp2 * LR * LRRR - 16.*
                                 cc24 * papb * p1p2 * RL * RLLL - 16.*cc24 * papb * p1p2 * LR * LRRR - 24.*cc24 *
                                 m1 * m2 * papb * RL * RLRL - 24.*cc24 * m1 * m2 * papb * LR * LRLR);

    cc->SetIeps(ieps + 2);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    me2 =  + ivs3v2 * ivu2s1 * (8.*cc24 * pap2 * pbp1 * RL * RLLL + 8.*cc24 * pap2 * pbp1 *
                                LR * LRRR - 8.*cc24 * pap1 * pbp2 * RL * RLLL - 8.*cc24 * pap1 * pbp2 * LR * LRRR
                                + 8.*cc24 * papb * p1p2 * RL * RLLL + 8.*cc24 * papb * p1p2 * LR * LRRR + 8.*
                                cc24 * m1 * m2 * papb * RL * RLRL + 8.*cc24 * m1 * m2 * papb * LR * LRLR);

    return real(me0 + me1 + me2);
}

double FI::MVutr3t(int ieps)
{
    complex<double> me0 = complex<double>(0., 0.);
    complex<double> me1 = complex<double>(0., 0.);

    complex<double> la[7];
    Npf *cc = new Npf(0., m2s - 2.*pap2, m2s, m1l, m2l, m3l, mul, ieps);

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
//   complex<double> cc11=la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    me0 =  + ivt1s2 * ivu2s1 * (2.*cc0 * pap2 * pbp1 * ml1 * ml3 * RR * RRRR + 2.*cc0 * pap2
                                * pbp1 * ml1 * ml3 * LL * LLLL + 2.*cc0 * pap1 * pbp2 * ml1 * ml3 * RR * RRRR + 2.*
                                cc0 * pap1 * pbp2 * ml1 * ml3 * LL * LLLL - 2.*cc0 * papb * p1p2 * ml1 * ml3 * RR *
                                RRRR - 2.*cc0 * papb * p1p2 * ml1 * ml3 * LL * LLLL + 2.*cc0 * m2 * pap2 * pbp1 *
                                ml1 * RR * RLRR + 2.*cc0 * m2 * pap2 * pbp1 * ml1 * LL * LRLL + 2.*cc0 * m2 * pap1 *
                                pbp2 * ml1 * RR * RLRR + 2.*cc0 * m2 * pap1 * pbp2 * ml1 * LL * LRLL - 2.*cc0 * m2 *
                                papb * p1p2 * ml1 * RR * RLRR - 2.*cc0 * m2 * papb * p1p2 * ml1 * LL * LRLL + 2.*cc0
                                * m1 * m2 * papb * ml1 * ml3 * RR * LRLR + 2.*cc0 * m1 * m2 * papb * ml1 * ml3 * LL * RLRL
                                + 2.*cc0 * m1 * m2s * papb * ml1 * RR * LLLR + 2.*cc0 * m1 * m2s *
                                papb * ml1 * LL * RRRL - 4.*cc12 * pow(pap2, 2) * pbp1 * RL * LLLL - 4.*cc12 *
                                pow(pap2, 2) * pbp1 * LR * RRRR - 4.*cc12 * pap1 * pap2 * pbp2 * RL * LLLL - 4.*
                                cc12 * pap1 * pap2 * pbp2 * LR * RRRR + 4.*cc12 * papb * pap2 * p1p2 * RL * LLLL +
                                4.*cc12 * papb * pap2 * p1p2 * LR * RRRR + 2.*cc12 * m2 * pap2 * pbp1 * ml3 * RL *
                                LRLL + 2.*cc12 * m2 * pap2 * pbp1 * ml3 * LR * RLRR + 2.*cc12 * m2 * pap2 * pbp1 *
                                ml1 * RR * RLRR + 2.*cc12 * m2 * pap2 * pbp1 * ml1 * LL * LRLL + 2.*cc12 * m2 * pap1
                                * pbp2 * ml3 * RL * LRLL);
    me0 +=  + ivt1s2 * ivu2s1 * (2.*cc12 * m2 * pap1 * pbp2 * ml3 * LR * RLRR + 2.*
                                 cc12 * m2 * pap1 * pbp2 * ml1 * RR * RLRR + 2.*cc12 * m2 * pap1 * pbp2 * ml1 * LL *
                                 LRLL - 2.*cc12 * m2 * papb * p1p2 * ml3 * RL * LRLL - 2.*cc12 * m2 * papb * p1p2 *
                                 ml3 * LR * RLRR - 2.*cc12 * m2 * papb * p1p2 * ml1 * RR * RLRR - 2.*cc12 * m2 * papb
                                 * p1p2 * ml1 * LL * LRLL + 2.*cc12 * m2s * pap2 * pbp1 * RL * LLLL + 2.*
                                 cc12 * m2s * pap2 * pbp1 * LR * RRRR + 2.*cc12 * m2s * pap1 * pbp2 *
                                 RL * LLLL + 2.*cc12 * m2s * pap1 * pbp2 * LR * RRRR - 2.*cc12 * pow(
                                     m2, 2) * papb * p1p2 * RL * LLLL - 2.*cc12 * m2s * papb * p1p2 * LR * RRRR
                                 - 4.*cc12 * m1 * m2 * papb * pap2 * RL * RLRL - 4.*cc12 * m1 * m2 * papb * pap2 * LR *
                                 LRLR + 2.*cc12 * m1 * m2s * papb * ml3 * RL * RRRL + 2.*cc12 * m1 * pow(
                                     m2, 2) * papb * ml3 * LR * LLLR + 2.*cc12 * m1 * m2s * papb * ml1 * RR * LLLR
                                 + 2.*cc12 * m1 * m2s * papb * ml1 * LL * RRRL + 2.*cc12 * m1 * pow(m2, 3) *
                                 papb * RL * RLRL + 2.*cc12 * m1 * pow(m2, 3) * papb * LR * LRLR + 2.*cc22 * pow(
                                     m2, 2) * pap2 * pbp1 * RL * LLLL + 2.*cc22 * m2s * pap2 * pbp1 * LR * RRRR
                                 + 2.*cc22 * m2s * pap1 * pbp2 * RL * LLLL + 2.*cc22 * m2s * pap1 *
                                 pbp2 * LR * RRRR);
    me0 +=  + ivt1s2 * ivu2s1 * (- 2.*cc22 * m2s * papb * p1p2 * RL * LLLL
                                 - 2.*cc22 * m2s * papb * p1p2 * LR * RRRR + 2.*cc22 * m1 * pow(m2, 3) *
                                 papb * RL * RLRL + 2.*cc22 * m1 * pow(m2, 3) * papb * LR * LRLR - 4.*cc23 * pow(
                                     pap2, 2) * pbp1 * RL * LLLL - 4.*cc23 * pow(pap2, 2) * pbp1 * LR * RRRR - 4.*
                                 cc23 * pap1 * pap2 * pbp2 * RL * LLLL - 4.*cc23 * pap1 * pap2 * pbp2 * LR * RRRR +
                                 4.*cc23 * papb * pap2 * p1p2 * RL * LLLL + 4.*cc23 * papb * pap2 * p1p2 * LR * RRRR
                                 - 4.*cc23 * m1 * m2 * papb * pap2 * RL * RLRL - 4.*cc23 * m1 * m2 * papb * pap2 * LR *
                                 LRLR + 8.*cc24 * pap2 * pbp1 * RL * LLLL + 8.*cc24 * pap2 * pbp1 * LR * RRRR + 8.
                                 *cc24 * pap1 * pbp2 * RL * LLLL + 8.*cc24 * pap1 * pbp2 * LR * RRRR - 8.*cc24 *
                                 papb * p1p2 * RL * LLLL - 8.*cc24 * papb * p1p2 * LR * RRRR + 8.*cc24 * m1 * m2 *
                                 papb * RL * RLRL + 8.*cc24 * m1 * m2 * papb * LR * LRLR);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    me1 =  + ivt1s2 * ivu2s1 * (- 4.*cc24 * pap2 * pbp1 * RL * LLLL - 4.*cc24 * pap2 *
                                pbp1 * LR * RRRR - 4.*cc24 * pap1 * pbp2 * RL * LLLL - 4.*cc24 * pap1 * pbp2 * LR *
                                RRRR + 4.*cc24 * papb * p1p2 * RL * LLLL + 4.*cc24 * papb * p1p2 * LR * RRRR - 4.
                                *cc24 * m1 * m2 * papb * RL * RLRL - 4.*cc24 * m1 * m2 * papb * LR * LRLR);

    return real(me0 + me1);
}

double FI::MVutr3u(int ieps)
{
    complex<double> me0 = complex<double>(0., 0.);
    complex<double> me1 = complex<double>(0., 0.);

    complex<double> la[7];
    Npf *cc = new Npf(0., m2s - 2.*pap2, m2s, m1l, m2l, m3l, mul, ieps);

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
//   complex<double> cc11=la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    me0 =  + ivu2s1 * ivu2s2 * (- 4.*cc0 * pap2 * pbp1 * ml1 * ml3 * RR * RRRR - 4.*cc0 *
                                pap2 * pbp1 * ml1 * ml3 * RR * LRLR - 4.*cc0 * pap2 * pbp1 * ml1 * ml3 * LL * RLRL -
                                4.*cc0 * pap2 * pbp1 * ml1 * ml3 * LL * LLLL - 4.*cc0 * m2 * pap2 * pbp1 * ml1 * RR *
                                RLRR - 4.*cc0 * m2 * pap2 * pbp1 * ml1 * RR * LLLR - 4.*cc0 * m2 * pap2 * pbp1 * ml1
                                * LL * RRRL - 4.*cc0 * m2 * pap2 * pbp1 * ml1 * LL * LRLL + 8.*cc12 * pow(pap2, 2)
                                * pbp1 * RL * RLRL + 8.*cc12 * pow(pap2, 2) * pbp1 * RL * LLLL + 8.*cc12 * pow(
                                    pap2, 2) * pbp1 * LR * RRRR + 8.*cc12 * pow(pap2, 2) * pbp1 * LR * LRLR - 4.*
                                cc12 * m2 * pap2 * pbp1 * ml3 * RL * RRRL - 4.*cc12 * m2 * pap2 * pbp1 * ml3 * RL *
                                LRLL - 4.*cc12 * m2 * pap2 * pbp1 * ml3 * LR * RLRR - 4.*cc12 * m2 * pap2 * pbp1 *
                                ml3 * LR * LLLR - 4.*cc12 * m2 * pap2 * pbp1 * ml1 * RR * RLRR - 4.*cc12 * m2 * pap2
                                * pbp1 * ml1 * RR * LLLR - 4.*cc12 * m2 * pap2 * pbp1 * ml1 * LL * RRRL - 4.*cc12 *
                                m2 * pap2 * pbp1 * ml1 * LL * LRLL - 4.*cc12 * m2s * pap2 * pbp1 * RL * RLRL
                                - 4.*cc12 * m2s * pap2 * pbp1 * RL * LLLL - 4.*cc12 * m2s * pap2 *
                                pbp1 * LR * RRRR - 4.*cc12 * m2s * pap2 * pbp1 * LR * LRLR - 4.*cc22 *
                                pow(m2, 2) * pap2 * pbp1 * RL * RLRL - 4.*cc22 * m2s * pap2 * pbp1 * RL *
                                LLLL);
    me0 +=  + ivu2s1 * ivu2s2 * (- 4.*cc22 * m2s * pap2 * pbp1 * LR * RRRR
                                 - 4.*cc22 * m2s * pap2 * pbp1 * LR * LRLR + 8.*cc23 * pow(pap2, 2) *
                                 pbp1 * RL * RLRL + 8.*cc23 * pow(pap2, 2) * pbp1 * RL * LLLL +
                                 8.*cc23 * pow(pap2, 2) * pbp1 * LR * RRRR +
                                 8.*cc23 * pow(pap2, 2) * pbp1 * LR * LRLR - 16.*
                                 cc24 * pap2 * pbp1 * RL * RLRL - 16.*cc24 * pap2 * pbp1 * RL * LLLL - 16.*cc24 *
                                 pap2 * pbp1 * LR * RRRR - 16.*cc24 * pap2 * pbp1 * LR * LRLR);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    me1 = + ivu2s1 * ivu2s2 * (8.*cc24 * pap2 * pbp1 * RL * RLRL + 8.*cc24 * pap2 * pbp1 *
                               RL * LLLL + 8.*cc24 * pap2 * pbp1 * LR * RRRR + 8.*cc24 * pap2 * pbp1 * LR * LRLR);

    return real(me0 + me1);
}

double FI::MVutr4s(int ieps)
{
    complex<double> me0 = complex<double>(0., 0.);
    complex<double> me1 = complex<double>(0., 0.);
    complex<double> me2 = complex<double>(0., 0.);

    complex<double> la[7];
    Npf *cc = new Npf(0., m2s - 2.*pap2, m1s, m1l, m2l, m3l, mul, ieps);

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
//   complex<double> cc11=la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    me0 =  + ivs3v2 * ivu2s1 * (8.*cc0 * pap2 * pbp1 * ml1 * ml3 * RR * RLLL + 8.*cc0 * pap2
                                * pbp1 * ml1 * ml3 * LL * LRRR + 8.*cc0 * m1 * pap2 * pbp1 * ml1 * RR * LLLL + 8.*cc0
                                * m1 * pap2 * pbp1 * ml1 * LL * RRRR + 4.*cc0 * m1 * m2 * papb * ml1 * ml3 * RR * RLRL
                                + 4.*cc0 * m1 * m2 * papb * ml1 * ml3 * LL * LRLR + 4.*cc0 * m1s * m2 * papb *
                                ml1 * RR * LLRL + 4.*cc0 * m1s * m2 * papb * ml1 * LL * RRLR - 16.*cc12 *
                                pap2 * pow(pbp1, 2) * RL * LRRR - 16.*cc12 * pap2 * pow(pbp1, 2) * LR * RLLL +
                                8.*cc12 * m1 * pap2 * pbp1 * ml3 * RL * RRRR + 8.*cc12 * m1 * pap2 * pbp1 * ml3 * LR *
                                LLLL + 8.*cc12 * m1 * pap2 * pbp1 * ml1 * RR * LLLL + 8.*cc12 * m1 * pap2 * pbp1 *
                                ml1 * LL * RRRR - 8.*cc12 * m1 * m2 * papb * pbp1 * RL * LRLR - 8.*cc12 * m1 * m2 *
                                papb * pbp1 * LR * RLRL + 8.*cc12 * m1s * pap2 * pbp1 * RL * LRRR + 8.*
                                cc12 * m1s * pap2 * pbp1 * LR * RLLL + 4.*cc12 * m1s * m2 * papb *
                                ml3 * RL * RRLR + 4.*cc12 * m1s * m2 * papb * ml3 * LR * LLRL + 4.*cc12 *
                                pow(m1, 2) * m2 * papb * ml1 * RR * LLRL + 4.*cc12 * m1s * m2 * papb * ml1 *
                                LL * RRLR + 4.*cc12 * pow(m1, 3) * m2 * papb * RL * LRLR + 4.*cc12 * pow(m1, 3) *
                                m2 * papb * LR * RLRL + 8.*cc22 * m1s * pap2 * pbp1 * RL * LRRR + 8.*cc22 *
                                pow(m1, 2) * pap2 * pbp1 * LR * RLLL);
    me0 +=  + ivs3v2 * ivu2s1 * (4.*cc22 * pow(m1, 3) * m2 * papb * RL * LRLR + 4.*
                                 cc22 * pow(m1, 3) * m2 * papb * LR * RLRL - 16.*cc23 * pap2 * pow(pbp1, 2) * RL *
                                 LRRR - 16.*cc23 * pap2 * pow(pbp1, 2) * LR * RLLL - 8.*cc23 * m1 * m2 * papb *
                                 pbp1 * RL * LRLR - 8.*cc23 * m1 * m2 * papb * pbp1 * LR * RLRL + 32.*cc24 * pap2 *
                                 pbp1 * RL * LRRR + 32.*cc24 * pap2 * pbp1 * LR * RLLL + 16.*cc24 * m1 * m2 * papb *
                                 RL * LRLR + 16.*cc24 * m1 * m2 * papb * LR * RLRL);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc22 = la[1];
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
//   cc11=la[0];
    cc12 = la[1];
    cc->GetNpf(la, 0);
    cc0 = la[0];

    me1 =  + ivs3v2 * ivu2s1 * (- 4.*cc0 * pap2 * pbp1 * ml1 * ml3 * RR * RLLL - 4.*cc0 *
                                pap2 * pbp1 * ml1 * ml3 * LL * LRRR + 4.*cc0 * pap1 * pbp2 * ml1 * ml3 * RR * RLLL +
                                4.*cc0 * pap1 * pbp2 * ml1 * ml3 * LL * LRRR - 4.*cc0 * papb * p1p2 * ml1 * ml3 * RR *
                                RLLL - 4.*cc0 * papb * p1p2 * ml1 * ml3 * LL * LRRR - 4.*cc0 * m1 * pap2 * pbp1 *
                                ml1 * RR * LLLL - 4.*cc0 * m1 * pap2 * pbp1 * ml1 * LL * RRRR + 4.*cc0 * m1 * pap1 *
                                pbp2 * ml1 * RR * LLLL + 4.*cc0 * m1 * pap1 * pbp2 * ml1 * LL * RRRR - 4.*cc0 * m1 *
                                papb * p1p2 * ml1 * RR * LLLL - 4.*cc0 * m1 * papb * p1p2 * ml1 * LL * RRRR - 4.*cc0
                                * m1 * m2 * papb * ml1 * ml3 * RR * RLRL - 4.*cc0 * m1 * m2 * papb * ml1 * ml3 * LL * LRLR
                                - 4.*cc0 * m1s * m2 * papb * ml1 * RR * LLRL - 4.*cc0 * m1s * m2 *
                                papb * ml1 * LL * RRLR + 8.*cc12 * pap2 * pow(pbp1, 2) * RL * LRRR + 8.*cc12 *
                                pap2 * pow(pbp1, 2) * LR * RLLL - 8.*cc12 * pap1 * pbp1 * pbp2 * RL * LRRR - 8.*
                                cc12 * pap1 * pbp1 * pbp2 * LR * RLLL + 8.*cc12 * papb * pbp1 * p1p2 * RL * LRRR +
                                8.*cc12 * papb * pbp1 * p1p2 * LR * RLLL - 4.*cc12 * m1 * pap2 * pbp1 * ml3 * RL *
                                RRRR - 4.*cc12 * m1 * pap2 * pbp1 * ml3 * LR * LLLL - 4.*cc12 * m1 * pap2 * pbp1 *
                                ml1 * RR * LLLL - 4.*cc12 * m1 * pap2 * pbp1 * ml1 * LL * RRRR + 4.*cc12 * m1 * pap1
                                * pbp2 * ml3 * RL * RRRR);
    me1 +=  + ivs3v2 * ivu2s1 * (4.*cc12 * m1 * pap1 * pbp2 * ml3 * LR * LLLL + 4.*
                                 cc12 * m1 * pap1 * pbp2 * ml1 * RR * LLLL + 4.*cc12 * m1 * pap1 * pbp2 * ml1 * LL *
                                 RRRR - 4.*cc12 * m1 * papb * p1p2 * ml3 * RL * RRRR - 4.*cc12 * m1 * papb * p1p2 *
                                 ml3 * LR * LLLL - 4.*cc12 * m1 * papb * p1p2 * ml1 * RR * LLLL - 4.*cc12 * m1 * papb
                                 * p1p2 * ml1 * LL * RRRR + 8.*cc12 * m1 * m2 * papb * pbp1 * RL * LRLR + 8.*cc12 * m1
                                 * m2 * papb * pbp1 * LR * RLRL - 4.*cc12 * m1s * pap2 * pbp1 * RL * LRRR - 4.
                                 *cc12 * m1s * pap2 * pbp1 * LR * RLLL + 4.*cc12 * m1s * pap1 * pbp2
                                 * RL * LRRR + 4.*cc12 * m1s * pap1 * pbp2 * LR * RLLL - 4.*cc12 * pow(
                                     m1, 2) * papb * p1p2 * RL * LRRR - 4.*cc12 * m1s * papb * p1p2 * LR * RLLL
                                 - 4.*cc12 * m1s * m2 * papb * ml3 * RL * RRLR - 4.*cc12 * m1s * m2 *
                                 papb * ml3 * LR * LLRL - 4.*cc12 * m1s * m2 * papb * ml1 * RR * LLRL - 4.*
                                 cc12 * m1s * m2 * papb * ml1 * LL * RRLR - 4.*cc12 * pow(m1, 3) * m2 * papb *
                                 RL * LRLR - 4.*cc12 * pow(m1, 3) * m2 * papb * LR * RLRL - 4.*cc22 * m1s *
                                 pap2 * pbp1 * RL * LRRR - 4.*cc22 * m1s * pap2 * pbp1 * LR * RLLL + 4.*
                                 cc22 * m1s * pap1 * pbp2 * RL * LRRR + 4.*cc22 * m1s * pap1 * pbp2 *
                                 LR * RLLL);
    me1 +=  + ivs3v2 * ivu2s1 * (- 4.*cc22 * m1s * papb * p1p2 * RL * LRRR
                                 - 4.*cc22 * m1s * papb * p1p2 * LR * RLLL - 4.*cc22 * pow(m1, 3) * m2 *
                                 papb * RL * LRLR - 4.*cc22 * pow(m1, 3) * m2 * papb * LR * RLRL + 8.*cc23 * pap2 *
                                 pow(pbp1, 2) * RL * LRRR + 8.*cc23 * pap2 * pow(pbp1, 2) * LR * RLLL - 8.*cc23
                                 * pap1 * pbp1 * pbp2 * RL * LRRR - 8.*cc23 * pap1 * pbp1 * pbp2 * LR * RLLL + 8.*
                                 cc23 * papb * pbp1 * p1p2 * RL * LRRR + 8.*cc23 * papb * pbp1 * p1p2 * LR * RLLL +
                                 8.*cc23 * m1 * m2 * papb * pbp1 * RL * LRLR + 8.*cc23 * m1 * m2 * papb * pbp1 * LR *
                                 RLRL - 32.*cc24 * pap2 * pbp1 * RL * LRRR - 32.*cc24 * pap2 * pbp1 * LR * RLLL
                                 + 16.*cc24 * pap1 * pbp2 * RL * LRRR + 16.*cc24 * pap1 * pbp2 * LR * RLLL - 16.*
                                 cc24 * papb * p1p2 * RL * LRRR - 16.*cc24 * papb * p1p2 * LR * RLLL - 24.*cc24 *
                                 m1 * m2 * papb * RL * LRLR - 24.*cc24 * m1 * m2 * papb * LR * RLRL);

    cc->SetIeps(ieps + 2);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    me2 = + ivs3v2 * ivu2s1 * (8.*cc24 * pap2 * pbp1 * RL * LRRR + 8.*cc24 * pap2 * pbp1 *
                               LR * RLLL - 8.*cc24 * pap1 * pbp2 * RL * LRRR - 8.*cc24 * pap1 * pbp2 * LR * RLLL
                               + 8.*cc24 * papb * p1p2 * RL * LRRR + 8.*cc24 * papb * p1p2 * LR * RLLL + 8.*
                               cc24 * m1 * m2 * papb * RL * LRLR + 8.*cc24 * m1 * m2 * papb * LR * RLRL);

    return real(me0 + me1 + me2);
}

double FI::MVutr4t(int ieps)
{
    complex<double> me0 = complex<double>(0., 0.);
    complex<double> me1 = complex<double>(0., 0.);

    complex<double> la[7];
    Npf *cc = new Npf(0., m2s - 2.*pap2, m1s, m1l, m2l, m3l, mul, ieps);

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
//   complex<double> cc11=la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    me0 =  + ivt1s2 * ivu2s1 * (2.*cc0 * pap2 * pbp1 * ml1 * ml3 * RR * RRRR + 2.*cc0 * pap2
                                * pbp1 * ml1 * ml3 * LL * LLLL + 2.*cc0 * pap1 * pbp2 * ml1 * ml3 * RR * RRRR + 2.*
                                cc0 * pap1 * pbp2 * ml1 * ml3 * LL * LLLL - 2.*cc0 * papb * p1p2 * ml1 * ml3 * RR *
                                RRRR - 2.*cc0 * papb * p1p2 * ml1 * ml3 * LL * LLLL + 2.*cc0 * m1 * pap2 * pbp1 *
                                ml1 * RR * LRRR + 2.*cc0 * m1 * pap2 * pbp1 * ml1 * LL * RLLL + 2.*cc0 * m1 * pap1 *
                                pbp2 * ml1 * RR * LRRR + 2.*cc0 * m1 * pap1 * pbp2 * ml1 * LL * RLLL - 2.*cc0 * m1 *
                                papb * p1p2 * ml1 * RR * LRRR - 2.*cc0 * m1 * papb * p1p2 * ml1 * LL * RLLL + 2.*cc0
                                * m1 * m2 * papb * ml1 * ml3 * RR * RLRL + 2.*cc0 * m1 * m2 * papb * ml1 * ml3 * LL * LRLR
                                + 2.*cc0 * m1s * m2 * papb * ml1 * RR * LLRL + 2.*cc0 * m1s * m2 *
                                papb * ml1 * LL * RRLR - 4.*cc12 * pap2 * pow(pbp1, 2) * RL * LLLL - 4.*cc12 *
                                pap2 * pow(pbp1, 2) * LR * RRRR - 4.*cc12 * pap1 * pbp1 * pbp2 * RL * LLLL - 4.*
                                cc12 * pap1 * pbp1 * pbp2 * LR * RRRR + 4.*cc12 * papb * pbp1 * p1p2 * RL * LLLL +
                                4.*cc12 * papb * pbp1 * p1p2 * LR * RRRR + 2.*cc12 * m1 * pap2 * pbp1 * ml3 * RL *
                                RLLL + 2.*cc12 * m1 * pap2 * pbp1 * ml3 * LR * LRRR + 2.*cc12 * m1 * pap2 * pbp1 *
                                ml1 * RR * LRRR + 2.*cc12 * m1 * pap2 * pbp1 * ml1 * LL * RLLL + 2.*cc12 * m1 * pap1
                                * pbp2 * ml3 * RL * RLLL);
    me0 += + ivt1s2 * ivu2s1 * (2.*cc12 * m1 * pap1 * pbp2 * ml3 * LR * LRRR + 2.*
                                cc12 * m1 * pap1 * pbp2 * ml1 * RR * LRRR + 2.*cc12 * m1 * pap1 * pbp2 * ml1 * LL *
                                RLLL - 2.*cc12 * m1 * papb * p1p2 * ml3 * RL * RLLL - 2.*cc12 * m1 * papb * p1p2 *
                                ml3 * LR * LRRR - 2.*cc12 * m1 * papb * p1p2 * ml1 * RR * LRRR - 2.*cc12 * m1 * papb
                                * p1p2 * ml1 * LL * RLLL - 4.*cc12 * m1 * m2 * papb * pbp1 * RL * LRLR - 4.*cc12 * m1
                                * m2 * papb * pbp1 * LR * RLRL + 2.*cc12 * m1s * pap2 * pbp1 * RL * LLLL + 2.
                                *cc12 * m1s * pap2 * pbp1 * LR * RRRR + 2.*cc12 * m1s * pap1 * pbp2
                                * RL * LLLL + 2.*cc12 * m1s * pap1 * pbp2 * LR * RRRR - 2.*cc12 * pow(
                                    m1, 2) * papb * p1p2 * RL * LLLL - 2.*cc12 * m1s * papb * p1p2 * LR * RRRR
                                + 2.*cc12 * m1s * m2 * papb * ml3 * RL * RRLR + 2.*cc12 * m1s * m2 *
                                papb * ml3 * LR * LLRL + 2.*cc12 * m1s * m2 * papb * ml1 * RR * LLRL + 2.*
                                cc12 * m1s * m2 * papb * ml1 * LL * RRLR + 2.*cc12 * pow(m1, 3) * m2 * papb *
                                RL * LRLR + 2.*cc12 * pow(m1, 3) * m2 * papb * LR * RLRL + 2.*cc22 * m1s *
                                pap2 * pbp1 * RL * LLLL + 2.*cc22 * m1s * pap2 * pbp1 * LR * RRRR + 2.*
                                cc22 * m1s * pap1 * pbp2 * RL * LLLL + 2.*cc22 * m1s * pap1 * pbp2 *
                                LR * RRRR);
    me0 += + ivt1s2 * ivu2s1 * (- 2.*cc22 * m1s * papb * p1p2 * RL * LLLL
                                - 2.*cc22 * m1s * papb * p1p2 * LR * RRRR + 2.*cc22 * pow(m1, 3) * m2 *
                                papb * RL * LRLR + 2.*cc22 * pow(m1, 3) * m2 * papb * LR * RLRL - 4.*cc23 * pap2 *
                                pow(pbp1, 2) * RL * LLLL - 4.*cc23 * pap2 * pow(pbp1, 2) * LR * RRRR - 4.*cc23
                                * pap1 * pbp1 * pbp2 * RL * LLLL - 4.*cc23 * pap1 * pbp1 * pbp2 * LR * RRRR + 4.*
                                cc23 * papb * pbp1 * p1p2 * RL * LLLL + 4.*cc23 * papb * pbp1 * p1p2 * LR * RRRR -
                                4.*cc23 * m1 * m2 * papb * pbp1 * RL * LRLR - 4.*cc23 * m1 * m2 * papb * pbp1 * LR *
                                RLRL + 8.*cc24 * pap2 * pbp1 * RL * LLLL + 8.*cc24 * pap2 * pbp1 * LR * RRRR + 8.
                                *cc24 * pap1 * pbp2 * RL * LLLL + 8.*cc24 * pap1 * pbp2 * LR * RRRR - 8.*cc24 *
                                papb * p1p2 * RL * LLLL - 8.*cc24 * papb * p1p2 * LR * RRRR + 8.*cc24 * m1 * m2 *
                                papb * RL * LRLR + 8.*cc24 * m1 * m2 * papb * LR * RLRL);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    me1 =  + ivt1s2 * ivu2s1 * (- 4.*cc24 * pap2 * pbp1 * RL * LLLL - 4.*cc24 * pap2 *
                                pbp1 * LR * RRRR - 4.*cc24 * pap1 * pbp2 * RL * LLLL - 4.*cc24 * pap1 * pbp2 * LR *
                                RRRR + 4.*cc24 * papb * p1p2 * RL * LLLL + 4.*cc24 * papb * p1p2 * LR * RRRR - 4.
                                *cc24 * m1 * m2 * papb * RL * LRLR - 4.*cc24 * m1 * m2 * papb * LR * RLRL);

    return real(me0 + me1);
}

double FI::MVutr4u(int ieps)
{
    complex<double> me0 = complex<double>(0., 0.);
    complex<double> me1 = complex<double>(0., 0.);

    complex<double> la[7];
    Npf *cc = new Npf(0., m2s - 2.*pap2, m1s, m1l, m2l, m3l, mul, ieps);

    cc->GetNpf(la, 2);
    complex<double> cc22 = la[1];
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
//   complex<double> cc11=la[0];
    complex<double> cc12 = la[1];
    cc->GetNpf(la, 0);
    complex<double> cc0 = la[0];

    me0 =  + ivu2s1 * ivu2s2 * (- 4.*cc0 * pap2 * pbp1 * ml1 * ml3 * RR * RLRL - 4.*cc0 *
                                pap2 * pbp1 * ml1 * ml3 * RR * RRRR - 4.*cc0 * pap2 * pbp1 * ml1 * ml3 * LL * LRLR -
                                4.*cc0 * pap2 * pbp1 * ml1 * ml3 * LL * LLLL - 4.*cc0 * m1 * pap2 * pbp1 * ml1 * RR *
                                LRRR - 4.*cc0 * m1 * pap2 * pbp1 * ml1 * RR * LLRL - 4.*cc0 * m1 * pap2 * pbp1 * ml1
                                * LL * RRLR - 4.*cc0 * m1 * pap2 * pbp1 * ml1 * LL * RLLL + 8.*cc12 * pap2 * pow(
                                    pbp1, 2) * RL * LRLR + 8.*cc12 * pap2 * pow(pbp1, 2) * RL * LLLL + 8.*cc12 *
                                pap2 * pow(pbp1, 2) * LR * RLRL + 8.*cc12 * pap2 * pow(pbp1, 2) * LR * RRRR - 4.
                                *cc12 * m1 * pap2 * pbp1 * ml3 * RL * RRLR - 4.*cc12 * m1 * pap2 * pbp1 * ml3 * RL *
                                RLLL - 4.*cc12 * m1 * pap2 * pbp1 * ml3 * LR * LRRR - 4.*cc12 * m1 * pap2 * pbp1 *
                                ml3 * LR * LLRL - 4.*cc12 * m1 * pap2 * pbp1 * ml1 * RR * LRRR - 4.*cc12 * m1 * pap2
                                * pbp1 * ml1 * RR * LLRL - 4.*cc12 * m1 * pap2 * pbp1 * ml1 * LL * RRLR - 4.*cc12 *
                                m1 * pap2 * pbp1 * ml1 * LL * RLLL - 4.*cc12 * m1s * pap2 * pbp1 * RL * LRLR
                                - 4.*cc12 * m1s * pap2 * pbp1 * RL * LLLL - 4.*cc12 * m1s * pap2 *
                                pbp1 * LR * RLRL - 4.*cc12 * m1s * pap2 * pbp1 * LR * RRRR - 4.*cc22 *
                                pow(m1, 2) * pap2 * pbp1 * RL * LRLR - 4.*cc22 * m1s * pap2 * pbp1 * RL *
                                LLLL);
    me0 += + ivu2s1 * ivu2s2 * (- 4.*cc22 * m1s * pap2 * pbp1 * LR * RLRL
                                - 4.*cc22 * m1s * pap2 * pbp1 * LR * RRRR + 8.*cc23 * pap2 * pow(
                                    pbp1, 2) * RL * LRLR + 8.*cc23 * pap2 * pow(pbp1, 2) * RL * LLLL + 8.*cc23 *
                                pap2 * pow(pbp1, 2) * LR * RLRL + 8.*cc23 * pap2 * pow(pbp1, 2) * LR * RRRR -
                                16.*cc24 * pap2 * pbp1 * RL * LRLR - 16.*cc24 * pap2 * pbp1 * RL * LLLL - 16.*
                                cc24 * pap2 * pbp1 * LR * RLRL - 16.*cc24 * pap2 * pbp1 * LR * RRRR);

    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc24 = la[3];
    delete cc;

    me1 =  + ivu2s1 * ivu2s2 * (8.*cc24 * pap2 * pbp1 * RL * LRLR + 8.*cc24 * pap2 * pbp1 *
                                RL * LLLL + 8.*cc24 * pap2 * pbp1 * LR * RLRL + 8.*cc24 * pap2 * pbp1 * LR * RRRR);

    return real(me0 + me1);
}

// Process: q + \bar{q} -> sl + \bar{sl}.

double FI::MVstr1sSL(int ieps)
{
    Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc23 = la[2];
    complex<double> cc24 = la[3];
    cc->GetNpf(la, 1);
    complex<double> cc11 = la[0];
    complex<double> cc12 = la[1];

    complex<double> me0 = (4.*cc24 + 4.*(cc11 + cc23) * papb)
                          * (8.*ivs3v1 * ivs3v2 * (LRLR * LR + LLLL * RL) * (2.*pap1 * pap2 - pap1 * m2s - pap2 * m1s));


    cc->SetIeps(ieps + 1);
    cc->GetNpf(la, 2);
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
    cc11 = la[0];
    cc12 = la[1];

    complex<double> me1 = (-8.*cc24 - 4.*(cc12 + cc23) * papb)
                          * (8.*ivs3v1 * ivs3v2 * (LRLR * LR + LLLL * RL) * (2.*pap1 * pap2 - pap1 * m2s - pap2 * m1s));


    cc->SetIeps(ieps + 2);
    cc->GetNpf(la, 2);
    cc23 = la[2];
    cc24 = la[3];
    cc->GetNpf(la, 1);
    cc11 = la[0];
    cc12 = la[1];
    delete cc;

    complex<double> me2 = (4.*cc24)
                          * (8.*ivs3v1 * ivs3v2 * (LRLR * LR + LLLL * RL) * (2.*pap1 * pap2 - pap1 * m2s - pap2 * m1s));


    return real(me0 + me1 + me2);
}


double FI::MVstr2sSL(int ieps)
{
    Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);
    complex<double> la[7];

    cc->GetNpf(la, 2);
    complex<double> cc24 = la[3];
    delete cc;


    complex<double> me0 = 2.* cc24
                          * (8.*ivs3v1 * ivs3v2 * (LRLR * LR + LLLL * RL) * (2.*pap1 * pap2 - pap1 * m2s - pap2 * m1s));

    return real(me0);
}
