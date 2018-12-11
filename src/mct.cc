// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Counter-term matrix elements.

#include <complex>
#include "fi.h"
#include "npf.h"

double FI::MVsct1t(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);
    std::complex<double> me1 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = - ivs3v1 * ivt1s2 * (4 * m2 * pbp1 * ED * RRLL + 4 * m2 * pbp1 * ED * LLRR + 8 * m1
                               * pbp2 * ED * RLRR + 8 * m1 * pbp2 * ED * LRLL);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;
    ED = -.5 * bb0 * ml1 * (LL + RR);

    me1 = - ivs3v1 * ivt1s2 * (- 4 * m2 * pbp1 * ED * RRLL - 4 * m2 * pbp1 * ED * LLRR - 4
                               * m1 * pbp2 * ED * RLRR - 4 * m1 * pbp2 * ED * LRLL);

    return std::real(me0 + me1);
}

double FI::MVsct1u(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);
    std::complex<double> me1 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivs3v1 * ivu2s2 * (8 * m2 * pbp1 * ED * RRLL + 8 * m2 * pbp1 * ED * LLRR + 4 * m1
                               * pbp2 * ED * RLRR + 4 * m1 * pbp2 * ED * LRLL);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;
    ED = -.5 * bb0 * ml1 * (LL + RR);

    me1 = + ivs3v1 * ivu2s2 * (- 4 * m2 * pbp1 * ED * RRLL - 4 * m2 * pbp1 * ED * LLRR - 4
                               * m1 * pbp2 * ED * RLRR - 4 * m1 * pbp2 * ED * LRLL);

    return std::real(me0 + me1);
}

double FI::MVtct1s(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);
    std::complex<double> me1 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivs3v2 * ivt1s1 * (- 4 * m2 * pbp1 * ED * RRLL - 4 * m2 * pbp1 * ED * LLRR - 8
                               * m1 * pbp2 * ED * RRRL - 8 * m1 * pbp2 * ED * LLLR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;
    ED = -.5 * bb0 * ml1 * (LL + RR);

    me1 = + ivs3v2 * ivt1s1 * (4 * m2 * pbp1 * ED * RRLL + 4 * m2 * pbp1 * ED * LLRR + 4 * m1
                               * pbp2 * ED * RRRL + 4 * m1 * pbp2 * ED * LLLR);

    return std::real(me0 + me1);
}

double FI::MVtct1t(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];
    delete bb;

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivt1s1 * ivt1s2 * (- 4 * m1 * pbp2 * ED * RLRR - 4 * m1 * pbp2 * ED * RRRL - 4
                               * m1 * pbp2 * ED * LRLL - 4 * m1 * pbp2 * ED * LLLR);

    return std::real(me0);
}

double FI::MVtct1u(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];
    delete bb;

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivt1s1 * ivu2s2 * (2 * m2 * pbp1 * ED * RRRL + 2 * m2 * pbp1 * ED * LLLR + 2 * m1
                               * pbp2 * ED * RLRR + 2 * m1 * pbp2 * ED * LRLL);

    return std::real(me0);
}

double FI::MVuct1s(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);
    std::complex<double> me1 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivs3v2 * ivu2s1 * (- 8 * m2 * pbp1 * ED * RRLL - 8 * m2 * pbp1 * ED * LLRR - 4
                               * m1 * pbp2 * ED * RRRL - 4 * m1 * pbp2 * ED * LLLR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;
    ED = -.5 * bb0 * ml1 * (LL + RR);

    me1 = + ivs3v2 * ivu2s1 * (4 * m2 * pbp1 * ED * RRLL + 4 * m2 * pbp1 * ED * LLRR + 4 * m1
                               * pbp2 * ED * RRRL + 4 * m1 * pbp2 * ED * LLLR);

    return std::real(me0);
}

double FI::MVuct1t(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];
    delete bb;

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivt1s2 * ivu2s1 * (- 2 * m2 * pbp1 * ED * RLRR - 2 * m2 * pbp1 * ED * LRLL - 2
                               * m1 * pbp2 * ED * RRRL - 2 * m1 * pbp2 * ED * LLLR);

    return std::real(me0);
}

double FI::MVuct1u(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];
    delete bb;

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivu2s1 * ivu2s2 * (4 * m2 * pbp1 * ED * RLRR + 4 * m2 * pbp1 * ED * RRRL + 4 * m2
                               * pbp1 * ED * LRLL + 4 * m2 * pbp1 * ED * LLLR);

    return std::real(me0);
}

double FI::MVsct2t(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);
    std::complex<double> me1 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivs3v1 * ivt1s2 * (- 8 * m2 * pap1 * ED * LRRR - 8 * m2 * pap1 * ED * RLLL - 4
                               * m1 * pap2 * ED * RRRR - 4 * m1 * pap2 * ED * LLLL);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;
    ED = -.5 * bb0 * ml1 * (LL + RR);

    me1 = + ivs3v1 * ivt1s2 * (4 * m2 * pap1 * ED * LRRR + 4 * m2 * pap1 * ED * RLLL + 4 * m1
                               * pap2 * ED * RRRR + 4 * m1 * pap2 * ED * LLLL);

    return std::real(me0 + me1);
}

double FI::MVsct2u(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);
    std::complex<double> me1 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivs3v1 * ivu2s2 * (4 * m2 * pap1 * ED * LRRR + 4 * m2 * pap1 * ED * RLLL + 8 * m1
                               * pap2 * ED * RRRR + 8 * m1 * pap2 * ED * LLLL);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;
    ED = -.5 * bb0 * ml1 * (LL + RR);

    me1 = + ivs3v1 * ivu2s2 * (- 4 * m2 * pap1 * ED * LRRR - 4 * m2 * pap1 * ED * RLLL - 4
                               * m1 * pap2 * ED * RRRR - 4 * m1 * pap2 * ED * LLLL);

    return std::real(me0 + me1);
}

double FI::MVtct2s(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);
    std::complex<double> me1 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivs3v2 * ivt1s1 * (- 8 * m2 * pap1 * ED * RRLR - 8 * m2 * pap1 * ED * LLRL - 4
                               * m1 * pap2 * ED * RRRR - 4 * m1 * pap2 * ED * LLLL);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;
    ED = -.5 * bb0 * ml1 * (LL + RR);

    me1 = + ivs3v2 * ivt1s1 * (4 * m2 * pap1 * ED * RRLR + 4 * m2 * pap1 * ED * LLRL + 4 * m1
                               * pap2 * ED * RRRR + 4 * m1 * pap2 * ED * LLLL);

    return std::real(me0 + me1);
}

double FI::MVtct2t(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];
    delete bb;

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivt1s1 * ivt1s2 * (- 4 * m2 * pap1 * ED * LRRR - 4 * m2 * pap1 * ED * RRLR - 4
                               * m2 * pap1 * ED * RLLL - 4 * m2 * pap1 * ED * LLRL);

    return std::real(me0);
}

double FI::MVtct2u(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];
    delete bb;

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivt1s1 * ivu2s2 * (2 * m2 * pap1 * ED * LRRR + 2 * m2 * pap1 * ED * RLLL + 2 * m1
                               * pap2 * ED * RRLR + 2 * m1 * pap2 * ED * LLRL);

    return std::real(me0);
}

double FI::MVuct2s(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);
    std::complex<double> me1 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivs3v2 * ivu2s1 * (- 4 * m2 * pap1 * ED * RRLR - 4 * m2 * pap1 * ED * LLRL - 8
                               * m1 * pap2 * ED * RRRR - 8 * m1 * pap2 * ED * LLLL);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 0);
    bb0 = la[0];
    delete bb;
    ED = -.5 * bb0 * ml1 * (LL + RR);

    me1 = + ivs3v2 * ivu2s1 * (4 * m2 * pap1 * ED * RRLR + 4 * m2 * pap1 * ED * LLRL + 4 * m1
                               * pap2 * ED * RRRR + 4 * m1 * pap2 * ED * LLLL);

    return std::real(me0 + me1);
}

double FI::MVuct2t(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];
    delete bb;

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivt1s2 * ivu2s1 * (- 2 * m2 * pap1 * ED * RRLR - 2 * m2 * pap1 * ED * LLRL - 2
                               * m1 * pap2 * ED * LRRR - 2 * m1 * pap2 * ED * RLLL);

    return std::real(me0);
}

double FI::MVuct2u(int ieps)
{
    std::complex<double> me0 = std::complex<double>(0., 0.);

    Npf * bb = new Npf(0., m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];
    delete bb;

    std::complex<double> ED = -.5 * bb0 * ml1 * (LL + RR);

    me0 = + ivu2s1 * ivu2s2 * (4 * m1 * pap2 * ED * LRRR + 4 * m1 * pap2 * ED * RRLR + 4 * m1
                               * pap2 * ED * RLLL + 4 * m1 * pap2 * ED * LLRL);

    return std::real(me0);
}

double FI::MVtct3s(int ieps)
{
    Npf * bb = new Npf(m1l, m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 2);
    std::complex<double> bb21 = la[0];
    std::complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    std::complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> me0 = 4.*ivs3v2 * std::pow(ivt1s1, 2) *
                               (4.*bb22 + (bb0 - 2.*bb1 + bb21) * m1l) *
                               (m1 * m2 * papb * (LRRR + RLLL) + 2.*pap1 * pbp2 * (LRLR + RLRL)) * RR;

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb21 = la[0];
    bb22 = la[1];
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];

    std::complex<double> me1 = -4.*ivs3v2 * std::pow(ivt1s1, 2) *
                               ((bb0 - 2.*bb1 + bb21) * m1l *
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

    std::complex<double> me2 = 8.*bb22 * ivs3v2 * std::pow(ivt1s1, 2) *
                               (LRRR * m1 * m2 * papb +
                                LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLLL +
                                p1p2 * papb * RLRL - pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL) * RR;

    return -real(me0 + me1 + me2);
}

double FI::MVtct3t(int ieps)
{
    Npf * bb = new Npf(m1l, m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 2);
    std::complex<double> bb21 = la[0];
    std::complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    std::complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> me0 = 4.*std::pow(ivt1s1, 2) * ivt1s2 *
                               (4.*bb22 + (bb0 - 2.*bb1 + bb21) * m1l) * pap1 * pbp2 * RR *
                               (LLLL + LRLR + RLRL + RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    std::complex<double> me1 = -8.*bb22 * std::pow(ivt1s1, 2) * ivt1s2 * pap1 * pbp2
                               * RR * (LLLL + LRLR + RLRL + RRRR);

    return -real(me0 + me1);
}

double FI::MVtct3u(int ieps)
{
    Npf * bb = new Npf(m1l, m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 2);
    std::complex<double> bb21 = la[0];
    std::complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    std::complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> me0 = -2.*std::pow(ivt1s1, 2) * ivu2s2 *
                               (4.*bb22 + (bb0 - 2.*bb1 + bb21) * m1l) * RR *
                               (LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                                m1 * m2 * papb * RLRL - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    std::complex<double> me1 = 4.*bb22 * std::pow(ivt1s1, 2) * ivu2s2 *
                               RR * (LRLR * m1 * m2 * papb +
                                     LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL -
                                     p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    return -real(me0 + me1);
}

double FI::MVuct3s(int ieps)
{
    Npf * bb = new Npf(m1l, m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 2);
    std::complex<double> bb21 = la[0];
    std::complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    std::complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> me0 = -4.*ivs3v2 * std::pow(ivu2s1, 2) *
                               (4.*bb22 + (bb0 - 2.*bb1 + bb21) * m1l) *
                               (2.*pap2 * pbp1 * (LRRR + RLLL) + m1 * m2 * papb * (LRLR + RLRL)) * RR;

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb21 = la[0];
    bb22 = la[1];
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];

    std::complex<double> me1 = 4.*ivs3v2 * std::pow(ivu2s1, 2) *
                               ((bb0 - 2.*bb1 + bb21) * m1l *
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

    std::complex<double> me2 = -8.*bb22 * ivs3v2 * std::pow(ivu2s1, 2) *
                               (LRLR * m1 * m2 * papb + (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) *
                                (LRRR + RLLL) + m1 * m2 * papb * RLRL) * RR;

    return -real(me0 + me1 + me2);
}

double FI::MVuct3t(int ieps)
{
    Npf * bb = new Npf(m1l, m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 2);
    std::complex<double> bb21 = la[0];
    std::complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    std::complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> me0 = -2.*ivt1s2 * std::pow(ivu2s1, 2) *
                               (4.*bb22 + (bb0 - 2.*bb1 + bb21) * m1l) * RR *
                               (LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                                m1 * m2 * papb * RLRL - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    std::complex<double> me1 = 4.*bb22 * ivt1s2 * std::pow(ivu2s1, 2) *
                               RR * (LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                                     m1 * m2 * papb * RLRL - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    return -real(me0 + me1);
}

double FI::MVuct3u(int ieps)
{
    Npf * bb = new Npf(m1l, m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 2);
    std::complex<double> bb21 = la[0];
    std::complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    std::complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> me0 = 4.*std::pow(ivu2s1, 2) * ivu2s2 *
                               (4.*bb22 + (bb0 - 2.*bb1 + bb21) * m1l) * pap2 * pbp1 * RR *
                               (LLLL + LRLR + RLRL + RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    std::complex<double> me1 = -8.*bb22 * std::pow(ivu2s1, 2) * ivu2s2 * pap2 * pbp1 *
                               RR * (LLLL + LRLR + RLRL + RRRR);

    return -real(me0 + me1);
}

double FI::MVtct4s(int ieps)
{
    ivt1s1 = 1. / (m1s - 2.*pap1 - m3l);
    ivt2s1 = 1. / (m1s - 2.*pap1 - m4l);

    Npf * bb = new Npf(m4l, m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 2);
    std::complex<double> bb21 = la[0];
    std::complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    std::complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> me0 = 8.*ivs3v2 * ivt1s1 * ivt2s1 *
                               (m1 * m2 * papb * (LRRR + RLLL) + 2 * pap1 * pbp2 * (LRLR + RLRL)) *
                               ((4.*bb22 + (bb1 + bb21) * m4l) * (LR + RL) +
                                bb0 * ml1 * ml2 * (LL + RR));

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb21 = la[0];
    bb22 = la[1];
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];

    std::complex<double> me1 = -8.*ivs3v2 * ivt1s1 * ivt2s1 *
                               (2.*bb22 * (LR + RL) * (3.*LRRR * m1 * m2 * papb + 2.*LRLR * p1p2 * papb -
                                       2.*LRLR * pap2 * pbp1 + 4.*LRLR * pap1 * pbp2 + 3.*m1 * m2 * papb * RLLL +
                                       2.*p1p2 * papb * RLRL - 2.*pap2 * pbp1 * RLRL + 4.*pap1 * pbp2 * RLRL) +
                                (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - LRLR * pap2 * pbp1 + LRLR * pap1 * pbp2 +
                                 m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL) *
                                ((bb1 + bb21) * m4l * (LR + RL) + bb0 * ml1 * ml2 * (LL + RR)));

    bb->SetIeps(ieps + 2);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    std::complex<double> me2 = 16.*bb22 * ivs3v2 * ivt1s1 * ivt2s1 *
                               (LR + RL) * (LRRR * m1 * m2 * papb + LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)
                                            + m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL);

    return -real(me0 + me1 + me2);
}

double FI::MVtct4t(int ieps)
{
    ivt1s1 = 1. / (m1s - 2.*pap1 - m3l);
    ivt2s1 = 1. / (m1s - 2.*pap1 - m4l);

    Npf * bb = new Npf(m4l, m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 2);
    std::complex<double> bb21 = la[0];
    std::complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    std::complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> me0 = 8.*ivt1s1 * ivt2s1 * ivt1s2 * pap1 * pbp2 *
                               ((4.*bb22 + (bb1 + bb21) * m4l) *
                                (LR + RL) + bb0 * ml1 * ml2 * (LL + RR)) * (LLLL + LRLR + RLRL + RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    std::complex<double> me1 = -16.*bb22 * ivt1s1 * ivt2s1 * ivt1s2 *
                               pap1 * pbp2 * (LR + RL) * (LLLL + LRLR + RLRL + RRRR);

    return -real(me0 + me1);
}

double FI::MVtct4u(int ieps)
{
    ivt1s1 = 1. / (m1s - 2.*pap1 - m3l);
    ivt2s1 = 1. / (m1s - 2.*pap1 - m4l);

    Npf * bb = new Npf(m4l, m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 2);
    std::complex<double> bb21 = la[0];
    std::complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    std::complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> me0 = -4.*ivt1s1 * ivt2s1 * ivu2s2 *
                               ((4.*bb22 + (bb1 + bb21) * m4l) * (LR + RL) +
                                bb0 * ml1 * ml2 * (LL + RR)) * (LRLR * m1 * m2 * papb +
                                        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL -
                                        p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    std::complex<double> me1 = 8.*bb22 * ivt1s1 * ivt2s1 * ivu2s2 *
                               (LR + RL) * (LRLR * m1 * m2 * papb +
                                            LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL -
                                            p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    return -real(me0 + me1);
}

double FI::MVuct4s(int ieps)
{
    ivu1s1 = 1. / (m2s - 2.*pap2 - m3l);
    ivu2s1 = 1. / (m2s - 2.*pap2 - m4l);

    Npf * bb = new Npf(m4l, m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 2);
    std::complex<double> bb21 = la[0];
    std::complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    std::complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> me0 = -8.*ivs3v2 * ivu1s1 * ivu2s1 *
                               (2.*pap2 * pbp1 * (LRRR + RLLL) + m1 * m2 * papb * (LRLR + RLRL)) *
                               ((4.*bb22 + (bb1 + bb21) * m4l) * (LR + RL) + bb0 * ml1 * ml2 * (LL + RR));

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb21 = la[0];
    bb22 = la[1];
    bb->GetNpf(la, 1);
    bb1 = la[0];
    bb->GetNpf(la, 0);
    bb0 = la[0];

    std::complex<double> me1 = 8.*ivs3v2 * ivu1s1 * ivu2s1 *
                               (2.*bb22 * (LR + RL) * (3.*LRLR * m1 * m2 * papb +
                                       2.*(p1p2 * papb + 2.*pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) +
                                       3.*m1 * m2 * papb * RLRL) +
                                (LRLR * m1 * m2 * papb + (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) +
                                 m1 * m2 * papb * RLRL) * ((bb1 + bb21) * m4l * (LR + RL) +
                                         bb0 * ml1 * ml2 * (LL + RR)));

    bb->SetIeps(ieps + 2);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    std::complex<double> me2 = -16.*bb22 * ivs3v2 * ivu1s1 * ivu2s1 *
                               (LR + RL) * (LRLR * m1 * m2 * papb +
                                            (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) + m1 * m2 * papb * RLRL);

    return -real(me0 + me1 + me2);
}

double FI::MVuct4t(int ieps)
{
    ivu1s1 = 1. / (m2s - 2.*pap2 - m3l);
    ivu2s1 = 1. / (m2s - 2.*pap2 - m4l);

    Npf * bb = new Npf(m4l, m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 2);
    std::complex<double> bb21 = la[0];
    std::complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    std::complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> me0 = -4.*ivt1s2 * ivu1s1 * ivu2s1 *
                               ((4.*bb22 + (bb1 + bb21) * m4l) * (LR + RL) +
                                bb0 * ml1 * ml2 * (LL + RR)) * (LRLR * m1 * m2 * papb +
                                        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL -
                                        p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    std::complex<double> me1 = 8.*bb22 * ivt1s2 * ivu1s1 * ivu2s1 *
                               (LR + RL) * (LRLR * m1 * m2 * papb +
                                            LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL -
                                            p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

    return -real(me0 + me1);
}

double FI::MVuct4u(int ieps)
{
    ivu1s1 = 1. / (m2s - 2.*pap2 - m3l);
    ivu2s1 = 1. / (m2s - 2.*pap2 - m4l);

    Npf * bb = new Npf(m4l, m1l, m2l, mul, ieps);
    std::complex<double> la[7];
    bb->GetNpf(la, 2);
    std::complex<double> bb21 = la[0];
    std::complex<double> bb22 = la[1];
    bb->GetNpf(la, 1);
    std::complex<double> bb1 = la[0];
    bb->GetNpf(la, 0);
    std::complex<double> bb0 = la[0];

    std::complex<double> me0 = 8.*ivu1s1 * ivu2s1 * ivu2s2 *
                               pap2 * pbp1 * ((4.*bb22 + (bb1 + bb21) * m4l) *
                                       (LR + RL) + bb0 * ml1 * ml2 * (LL + RR)) * (LLLL + LRLR + RLRL + RRRR);

    bb->SetIeps(ieps + 1);
    bb->GetNpf(la, 2);
    bb22 = la[1];
    delete bb;

    std::complex<double> me1 = -16.*bb22 * ivu1s1 * ivu2s1 * ivu2s2 *
                               pap2 * pbp1 * (LR + RL) * (LLLL + LRLR + RLRL + RRRR);

    return -real(me0 + me1);
}
