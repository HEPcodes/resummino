// This file is part of Resummino.
//
// Copyright 2008-2011 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Box matrix elements.
//
// See mbu.cc and mbn.cc for notation details.

#include <complex>
#include "fi.h"
#include "npf.h"

// Process: q + \bar{q} -> f + \bar{f}.

double FI::MVtbo1s(int ieps)
{

    // Constructs general four point function.
    Npf *dd = new Npf(0., 0., m1s, m2s, 2.*papb, m1s - 2.*pap1,
                      m1l, m2l, m3l, m4l, mul, ieps);

    std::complex<double> la[7];
    dd->GetNpf(la, 2);
    std::complex<double> dd23 = la[2];
    std::complex<double> dd24 = la[3];
    std::complex<double> dd25 = la[4];
    std::complex<double> dd26 = la[5];
    std::complex<double> dd27 = la[6];

    dd->GetNpf(la, 1);
    std::complex<double> dd11 = la[0];
    std::complex<double> dd12 = la[1];
    std::complex<double> dd13 = la[2];

    // Part of scalar integrals of O(1/eps^0) times squared matrix element of
    // order O(eps^0).
    std::complex<double> me0
        = 8.*ivs3v2 *
          (-2.*LL * m1 * m2 * (dd26 * pap1 * papb - (dd11 + dd24) * std::pow(papb, 2) -
                               dd23 * pap1 * pbp1 + dd25 * papb * pbp1) * RLLL +
           LL * ((dd23 * m1s - 2.*dd26 * pap1) * (p1p2 * papb - pap2 * pbp1) +
                 (dd23 * m1s * pap1 - 2.*dd26 * std::pow(pap1, 2) - 2.*dd25 * m1s * papb +
                  4.*(dd11 + dd24) * pap1 * papb) * pbp2) * RLRL +
           (dd23 * (2.*LRRR * m1 * m2 * pap1 * pbp1 + LRLR * m1s * (p1p2 * papb - pap2 * pbp1 +
                    pap1 * pbp2)) + 2.*(LRRR * m1 * m2 * papb * ((dd11 + dd24) * papb - dd25 * pbp1) -
                                        LRLR * (dd25 * m1s - 2.*(dd11 + dd24) * pap1) * papb * pbp2 -
                                        dd26 * pap1 * (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - LRLR * pap2 * pbp1 +
                                                LRLR * pap1 * pbp2))) * RR + 2.*dd27 * (LL * m1 * m2 * papb * RLLL +
                                                        2.*LL * pap1 * pbp2 * RLRL + LRRR * m1 * m2 * papb * RR + 2.*LRLR * pap1 * pbp2 * RR) -
           2.*dd13 * pap1 * (LL * (m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL +
                                   pap1 * pbp2 * RLRL) + (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - LRLR * pap2 * pbp1 +
                                           LRLR * pap1 * pbp2) * RR));

    dd->SetIeps(ieps + 1); // Gets four-point function prop. to O(1/epsilon^1).
    dd->GetNpf(la, 2);
    dd23 = la[2];
    dd24 = la[3];
    dd25 = la[4];
    dd26 = la[5];
    dd27 = la[6];
    dd->GetNpf(la, 1);
    dd11 = la[0];
    dd12 = la[1];
    dd13 = la[2];

    // Part of scalar integrals of O(1/eps^1) times squared matrix element of O(eps^1).
    std::complex<double> me1 = -8.*ivs3v2 *
                               (dd23 * (LL * m1 * m2 * (m1s * papb + 2.*pap1 * pbp1) * RLLL +
                                        2.*LL * m1s * p1p2 * papb * RLRL + (2.*LRLR * m1s * p1p2 * papb +
                                                LRRR * m1 * m2 * (m1s * papb + 2.*pap1 * pbp1)) * RR) +
                                2.*dd27 * (LL * (3.*m1 * m2 * papb * RLLL + (p1p2 * papb + 3.*pap2 * pbp1 + pap1 * pbp2) *
                                           RLRL) + (3.*LRRR * m1 * m2 * papb + LRLR * (p1p2 * papb + 3.*pap2 * pbp1 +
                                                   pap1 * pbp2)) * RR) +
                                2.*(LL * m1 * m2 * papb * ((dd11 + dd12 + 2.*dd24) * papb - 2.*dd25 * pbp1) * RLLL +
                                    LL * (-(dd25 * p1p2 * papb * pbp1) + 2.*dd12 * pap2 * papb * pbp1 -
                                          dd25 * pap2 * std::pow(pbp1, 2) -
                                          dd25 * m1s * papb * pbp2 + dd25 * pap1 * pbp1 * pbp2 +
                                          dd11 * papb * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                                          dd24 * papb * (p1p2 * papb + pap2 * pbp1 + pap1 * pbp2)) * RLRL +
                                    dd13 * LL * (-(m1 * m2 * papb * (pap1 + pbp1) * RLLL) +
                                            (m1s * pap2 * papb - pap2 * std::pow(pbp1, 2) - p1p2 * papb * (2.*pap1 + pbp1) +
                                                    pap1 * pbp1 * pbp2) * RLRL) - dd13 * (LRRR * m1 * m2 * papb * (pap1 + pbp1) +
                                                            LRLR * (-(m1s * pap2 * papb) + pap2 * std::pow(pbp1, 2) +
                                                                    p1p2 * papb * (2.*pap1 + pbp1) - pap1 * pbp1 * pbp2)) * RR +
                                    ((dd24 * papb - dd25 * pbp1) * (2.*LRRR * m1 * m2 * papb + LRLR * p1p2 * papb +
                                            LRLR * pap2 * pbp1) + dd12 * papb * (LRRR * m1 * m2 * papb + 2.*LRLR * pap2 * pbp1) +
                                     LRLR * (-(dd25 * m1s * papb) + dd24 * pap1 * papb + dd25 * pap1 * pbp1) * pbp2 +
                                     dd11 * papb * (LRRR * m1 * m2 * papb + LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2))) *
                                    RR - 2.*dd26 * pap1 * papb * (LL * m1 * m2 * RLLL + LL * p1p2 * RLRL + LRRR * m1 * m2 * RR +
                                            LRLR * p1p2 * RR)));

    dd->SetIeps(ieps + 2); // Gets three-point function prop. to O(1/epsilon^2).
    dd->GetNpf(la, 2);
    dd23 = la[2];
    dd24 = la[3];
    dd25 = la[4];
    dd26 = la[5];
    dd27 = la[6];
    dd->GetNpf(la, 1);
    dd11 = la[0];
    dd12 = la[1];
    dd13 = la[2];
    delete dd;

    // Part of scalar integrals of O(1/eps^2) times squared matrix element of O(eps^2).
    std::complex<double> me2 = 8.*ivs3v2 *
                               (dd27 * (6.*LL * m1 * m2 * papb * RLLL + 4.*LL * p1p2 * papb * RLRL +
                                        4.*LL * pap2 * pbp1 * RLRL + 6.*LRRR * m1 * m2 * papb * RR + 4.*LRLR * p1p2 * papb * RR +
                                        4.*LRLR * pap2 * pbp1 * RR) +
                                2.*papb * (-(LL * m1 * m2 * (dd26 * pap1 - (dd12 + dd24) * papb + (dd13 + dd25) * pbp1) *
                                           RLLL) + LL * (-(dd26 * m1s * pap2) + (dd12 + dd24) * p1p2 * papb -
                                                   2.*(dd13 + dd25) * p1p2 * pbp1 + (dd13 + dd25) * m1s * pbp2 +
                                                   (dd12 + dd24) * (pap2 * pbp1 - pap1 * pbp2)) * RLRL +
                                           (-(dd26 * (LRRR * m1 * m2 * pap1 + LRLR * m1s * pap2)) + dd12 * LRRR * m1 * m2 * papb +
                                            dd24 * LRRR * m1 * m2 * papb + dd12 * LRLR * p1p2 * papb + dd24 * LRLR * p1p2 * papb -
                                            dd13 * LRRR * m1 * m2 * pbp1 - dd25 * LRRR * m1 * m2 * pbp1 - 2.*dd13 * LRLR * p1p2 * pbp1 -
                                            2.*dd25 * LRLR * p1p2 * pbp1 + dd12 * LRLR * pap2 * pbp1 + dd24 * LRLR * pap2 * pbp1 +
                                            LRLR * ((dd13 + dd25) * m1s - (dd12 + dd24) * pap1) * pbp2) * RR) +
                                dd23 * m1s * (LL * (m1 * m2 * papb * RLLL + (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RLRL) +
                                        (LRRR * m1 * m2 * papb + LRLR * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) * RR));

    return std::real(me0 + me1 + me2);
}

double FI::MVtbo1t(int ieps)
{
    Npf * dd = new Npf(0., 0., m1s, m2s, 2.*papb, m1s - 2.*pap1,
                       m1l, m2l, m3l, m4l, mul, ieps);
    std::complex<double> la[7];
    dd->GetNpf(la, 2);
    std::complex<double> dd23 = la[2];
    std::complex<double> dd24 = la[3];
    std::complex<double> dd25 = la[4];
    std::complex<double> dd26 = la[5];
    std::complex<double> dd27 = la[6];

    dd->GetNpf(la, 1);
    std::complex<double> dd11 = la[0];
    std::complex<double> dd12 = la[1];
    std::complex<double> dd13 = la[2];

    std::complex<double> me0 = 4.*ivt1s2 *
                               (dd23 * m1s * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) *
                                (LLLL * LR + LL * RLRL + LRLR * RR + RL * RRRR) +
                                4.*dd27 * (LLLL * LR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) + LL * pap1 * pbp2 * RLRL +
                                           LRLR * pap1 * pbp2 * RR + (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RL * RRRR) -
                                2.*(-(dd12 * LLLL * LR * p1p2 * std::pow(papb, 2)) -
                                    dd24 * LLLL * LR * p1p2 * std::pow(papb, 2) +
                                    dd25 * LLLL * LR * p1p2 * papb * pbp1 + dd12 * LLLL * LR * pap2 * papb * pbp1 +
                                    dd24 * LLLL * LR * pap2 * papb * pbp1 - dd25 * LLLL * LR * pap2 * std::pow(pbp1, 2) -
                                    2.*dd11 * LLLL * LR * pap1 * papb * pbp2 + dd12 * LLLL * LR * pap1 * papb * pbp2 -
                                    dd24 * LLLL * LR * pap1 * papb * pbp2 + dd25 * LLLL * LR * pap1 * pbp1 * pbp2 +
                                    dd25 * LL * m1s * papb * pbp2 * RLRL - 2.*dd11 * LL * pap1 * papb * pbp2 * RLRL -
                                    2.*dd24 * LL * pap1 * papb * pbp2 * RLRL + dd25 * LRLR * m1s * papb * pbp2 * RR -
                                    2.*dd11 * LRLR * pap1 * papb * pbp2 * RR - 2.*dd24 * LRLR * pap1 * papb * pbp2 * RR -
                                    (((dd12 + dd24) * papb - dd25 * pbp1) * (p1p2 * papb - pap2 * pbp1) +
                                     pap1 * ((2.*dd11 - dd12 + dd24) * papb - dd25 * pbp1) * pbp2) * RL * RRRR +
                                    dd26 * pap1 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * (LLLL * LR + LL * RLRL +
                                            LRLR * RR + RL * RRRR) +
                                    dd13 * (LLLL * LR * ((pap1 + pbp1) * (p1p2 * papb - pap2 * pbp1) +
                                            (-(m1s * papb) + pap1 * (pap1 + pbp1)) * pbp2) +
                                            pap1 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * (LL * RLRL + LRLR * RR) +
                                            ((pap1 + pbp1) * (p1p2 * papb - pap2 * pbp1) +
                                                    (-(m1s * papb) + pap1 * (pap1 + pbp1)) * pbp2) * RL * RRRR)));

    dd->SetIeps(ieps + 1);
    dd->GetNpf(la, 2);
    dd23 = la[2];
    dd24 = la[3];
    dd25 = la[4];
    dd26 = la[5];
    dd27 = la[6];
    dd->GetNpf(la, 1);
    dd11 = la[0];
    dd12 = la[1];
    dd13 = la[2];

    std::complex<double> me1 = -8.*ivt1s2 *
                               ((((dd12 + dd24) * papb - (dd13 + dd25) * pbp1) *
                                 (p1p2 * papb - pap2 * pbp1) + (((dd13 + dd25) * m1s - (dd12 + dd24) * pap1) *
                                         papb - (dd13 + dd25) * pap1 * pbp1) * pbp2) * (LLLL * LR - LL * RLRL - LRLR * RR +
                                                 RL * RRRR) + dd27 * (LLLL * LR * (3.*p1p2 * papb - 3.*pap2 * pbp1 + pap1 * pbp2) -
                                                         (3.*p1p2 * papb - 3.*pap2 * pbp1 - pap1 * pbp2) * (LL * RLRL + LRLR * RR) +
                                                         (3.*p1p2 * papb - 3.*pap2 * pbp1 + pap1 * pbp2) * RL * RRRR));

    dd->SetIeps(ieps + 2);
    dd->GetNpf(la, 2);
    dd27 = la[6];
    delete dd;

    std::complex<double> me2 = 16.*dd27 * ivt1s2 * (p1p2 * papb - pap2 * pbp1)
                               * (LLLL * LR - LL * RLRL - LRLR * RR + RL * RRRR);

    return std::real(me0 + me1 + me2);
}

double FI::MVtbo1u(int ieps)
{
    // LL<->RL; RR<->LR;
    Npf * dd = new Npf(0., 0., m1s, m2s, 2.*papb, m1s - 2.*pap1,
                       m1l, m2l, m3l, m4l, mul, ieps);
    std::complex<double> la[7];
    dd->GetNpf(la, 2);
    std::complex<double> dd23 = la[2];
    std::complex<double> dd24 = la[3];
    std::complex<double> dd25 = la[4];
    std::complex<double> dd26 = la[5];
    std::complex<double> dd27 = la[6];
    dd->GetNpf(la, 1);
    std::complex<double> dd11 = la[0];
    std::complex<double> dd12 = la[1];
    std::complex<double> dd13 = la[2];

    std::complex<double> me0 = 4.*ivu2s2 *
                               (2.*dd26 * LR * LRLR * m1 * m2 * pap1 * papb - 2.*dd11 * LR * LRLR * m1 * m2 * std::pow(papb, 2) -
                                2.*dd24 * LR * LRLR * m1 * m2 * std::pow(papb, 2) - 2.*dd23 * LR * LRLR * m1 * m2 * pap1 * pbp1 +
                                2.*dd25 * LR * LRLR * m1 * m2 * papb * pbp1 + 2.*dd26 * m1 * m2 * pap1 * papb * RL * RLRL -
                                2.*dd11 * m1 * m2 * std::pow(papb, 2) * RL * RLRL -
                                2.*dd24 * m1 * m2 * std::pow(papb, 2) * RL * RLRL -
                                2.*dd23 * m1 * m2 * pap1 * pbp1 * RL * RLRL + 2.*dd25 * m1 * m2 * papb * pbp1 * RL * RLRL +
                                dd23 * LLLL * m1s * p1p2 * papb * RR - 2.*dd26 * LLLL * p1p2 * pap1 * papb * RR +
                                2.*dd11 * LLLL * p1p2 * std::pow(papb, 2) * RR +
                                2.*dd24 * LLLL * p1p2 * std::pow(papb, 2) * RR +
                                dd23 * LLLL * m1s * pap2 * pbp1 * RR - 2.*dd26 * LLLL * pap1 * pap2 * pbp1 * RR -
                                2.*dd25 * LLLL * p1p2 * papb * pbp1 * RR - 2.*dd11 * LLLL * pap2 * papb * pbp1 * RR +
                                4.*dd12 * LLLL * pap2 * papb * pbp1 * RR + 2.*dd24 * LLLL * pap2 * papb * pbp1 * RR -
                                2.*dd25 * LLLL * pap2 * std::pow(pbp1, 2) * RR - dd23 * LLLL * m1s * pap1 * pbp2 * RR +
                                2.*dd26 * LLLL * std::pow(pap1, 2) * pbp2 * RR - 2.*dd11 * LLLL * pap1 * papb * pbp2 * RR -
                                2.*dd24 * LLLL * pap1 * papb * pbp2 * RR + 2.*dd25 * LLLL * pap1 * pbp1 * pbp2 * RR +
                                dd23 * LL * m1s * p1p2 * papb * RRRR - 2.*dd26 * LL * p1p2 * pap1 * papb * RRRR +
                                2.*dd11 * LL * p1p2 * std::pow(papb, 2) * RRRR +
                                2.*dd24 * LL * p1p2 * std::pow(papb, 2) * RRRR +
                                dd23 * LL * m1s * pap2 * pbp1 * RRRR - 2.*dd26 * LL * pap1 * pap2 * pbp1 * RRRR -
                                2.*dd25 * LL * p1p2 * papb * pbp1 * RRRR - 2.*dd11 * LL * pap2 * papb * pbp1 * RRRR +
                                4.*dd12 * LL * pap2 * papb * pbp1 * RRRR + 2.*dd24 * LL * pap2 * papb * pbp1 * RRRR -
                                2.*dd25 * LL * pap2 * std::pow(pbp1, 2) * RRRR - dd23 * LL * m1s * pap1 * pbp2 * RRRR +
                                2.*dd26 * LL * std::pow(pap1, 2) * pbp2 * RRRR - 2.*dd11 * LL * pap1 * papb * pbp2 * RRRR -
                                2.*dd24 * LL * pap1 * papb * pbp2 * RRRR + 2.*dd25 * LL * pap1 * pbp1 * pbp2 * RRRR -
                                2.*dd27 * (LR * LRLR * m1 * m2 * papb + m1 * m2 * papb * RL * RLRL -
                                           2.*(p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LLLL * RR + LL * RRRR)) +
                                2.*dd13 * (LR * LRLR * m1 * m2 * pap1 * papb + m1 * m2 * pap1 * papb * RL * RLRL -
                                           (p1p2 * papb * (pap1 + pbp1) + pap2 * (-(m1s * papb) + pbp1 * (pap1 + pbp1)) -
                                            pap1 * (pap1 + pbp1) * pbp2) * (LLLL * RR + LL * RRRR)));

    dd->SetIeps(ieps + 1);
    dd->GetNpf(la, 2);
    dd23 = la[2];
    dd24 = la[3];
    dd25 = la[4];
    dd26 = la[5];
    dd27 = la[6];
    dd->GetNpf(la, 1);
    dd11 = la[0];
    dd12 = la[1];
    dd13 = la[2];

    std::complex<double> me1 = -4.*ivu2s2 *
                               (-2.*(dd26 * m1s * pap2 * papb + pbp1 * (dd13 * p1p2 * papb + dd25 * p1p2 * papb -
                                       2.*dd12 * pap2 * papb - 2.*dd24 * pap2 * papb + dd13 * pap2 * pbp1 + dd25 * pap2 * pbp1 -
                                       (dd13 + dd25) * pap1 * pbp2)) * (LLLL * RR + LL * RRRR) +
                                dd23 * (LR * LRLR * m1 * m2 * (m1s * papb - 2.*pap1 * pbp1) +
                                        m1 * m2 * (m1s * papb - 2.*pap1 * pbp1) * RL * RLRL +
                                        m1s * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LLLL * RR + LL * RRRR)) +
                                2.*dd27 * (LR * LRLR * m1 * m2 * papb + m1 * m2 * papb * RL * RLRL +
                                           (p1p2 * papb + 3.*pap2 * pbp1 - pap1 * pbp2) * (LLLL * RR + LL * RRRR)));

    dd->SetIeps(ieps + 2);
    dd->GetNpf(la, 2);
    dd27 = la[6];
    delete dd;

    std::complex<double> me2 = 8.*dd27 * ivu2s2 *
                               (LR * LRLR * m1 * m2 * papb + m1 * m2 * papb * RL * RLRL -
                                (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) * (LLLL * RR + LL * RRRR));

    return std::real(me0 + me1 + me2);
}

double FI::MVubo1s(int ieps)
{
    Npf * dd = new Npf(0., 0., m2s, m1s, 2.*papb, m2s - 2.*pap2,
                       m1l, m2l, m3l, m4l, mul, ieps);
    std::complex<double> la[7];
    dd->GetNpf(la, 2);
    std::complex<double> dd23 = la[2];
    std::complex<double> dd24 = la[3];
    std::complex<double> dd25 = la[4];
    std::complex<double> dd26 = la[5];
    std::complex<double> dd27 = la[6];
    dd->GetNpf(la, 1);
    std::complex<double> dd11 = la[0];
    std::complex<double> dd12 = la[1];
    std::complex<double> dd13 = la[2];

    std::complex<double> me0 = -8.*ivs3v2 *
                               (dd23 * LL * (m2s * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RLLL +
                                       2.*m1 * m2 * pap2 * pbp2 * RLRL) + dd23 * (2.*LRLR * m1 * m2 * pap2 * pbp2 +
                                               LRRR * m2s * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) * RR +
                                2.*(-(LL * ((dd25 * m2s * papb - 2.*pap2 * (dd27 + (dd11 + dd24) * papb)) * pbp1 +
                                            dd26 * pap2 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) * RLLL) +
                                    LL * m1 * m2 * papb * (dd27 - dd26 * pap2 + (dd11 + dd24) * papb - dd25 * pbp2) * RLRL +
                                    (dd27 * (LRLR * m1 * m2 * papb + 2.*LRRR * pap2 * pbp1) +
                                     papb * (-(dd25 * LRRR * m2s * pbp1) + (dd11 + dd24) * (LRLR * m1 * m2 * papb +
                                             2.*LRRR * pap2 * pbp1) - dd25 * LRLR * m1 * m2 * pbp2) -
                                     dd26 * pap2 * (LRLR * m1 * m2 * papb + LRRR * p1p2 * papb + LRRR * pap2 * pbp1 -
                                             LRRR * pap1 * pbp2)) * RR - dd13 * pap2 *
                                    (LL * (p1p2 * papb * RLLL + pap2 * pbp1 * RLLL - pap1 * pbp2 * RLLL +
                                           m1 * m2 * papb * RLRL) + (LRLR * m1 * m2 * papb + LRRR * (p1p2 * papb + pap2 * pbp1 -
                                                   pap1 * pbp2)) * RR)));

    dd->SetIeps(ieps + 1);
    dd->GetNpf(la, 2);
    dd23 = la[2];
    dd24 = la[3];
    dd25 = la[4];
    dd26 = la[5];
    dd27 = la[6];
    dd->GetNpf(la, 1);
    dd11 = la[0];
    dd12 = la[1];
    dd13 = la[2];

    std::complex<double> me1 = 8.*ivs3v2 *
                               (2.*dd27 * (LL * (p1p2 * papb * RLLL + pap2 * pbp1 * RLLL + 3.*pap1 * pbp2 * RLLL +
                                           3.*m1 * m2 * papb * RLRL) + (3.*LRLR * m1 * m2 * papb +
                                                   LRRR * (p1p2 * papb + pap2 * pbp1 + 3.*pap1 * pbp2)) * RR) +
                                dd23 * (2.*LL * m2s * p1p2 * papb * RLLL + LL * m1 * m2 * (m2s * papb + 2.*pap2 * pbp2) * RLRL +
                                        (2.*LRRR * m2s * p1p2 * papb + LRLR * m1 * m2 * (m2s * papb + 2.*pap2 * pbp2)) * RR) +
                                2.*(LL * (papb * (-(dd25 * m2s * pbp1) + (dd11 + dd24) * (p1p2 * papb + pap2 * pbp1)) -
                                          ((dd11 - 2.*dd12 - dd24) * pap1 * papb + dd25 * (p1p2 * papb - pap2 * pbp1)) *
                                          pbp2 - dd25 * pap1 * std::pow(pbp2, 2)) * RLLL + LL * m1 * m2 * papb *
                                    ((dd11 + dd12 + 2.*dd24) * papb - 2.*dd25 * pbp2) * RLRL +
                                    dd13 * LL * (m2s * pap1 * papb * RLLL - p1p2 * papb * (2.*pap2 + pbp2) * RLLL +
                                            pbp2 * (pap2 * pbp1 - pap1 * pbp2) * RLLL - m1 * m2 * papb * (pap2 + pbp2) * RLRL) +
                                    (papb * (dd12 * LRLR * m1 * m2 * papb + 2.*dd24 * LRLR * m1 * m2 * papb +
                                            dd24 * LRRR * p1p2 * papb - dd25 * LRRR * m2s * pbp1 + dd24 * LRRR * pap2 * pbp1 +
                                            dd11 * (LRLR * m1 * m2 * papb + LRRR * p1p2 * papb + LRRR * pap2 * pbp1)) -
                                     ((dd11 - 2.*dd12 - dd24) * LRRR * pap1 * papb +
                                      dd25 * (2.*LRLR * m1 * m2 * papb + LRRR * p1p2 * papb - LRRR * pap2 * pbp1)) * pbp2 -
                                     dd25 * LRRR * pap1 * std::pow(pbp2, 2)) * RR +
                                    dd13 * (-(LRLR * m1 * m2 * papb * (pap2 + pbp2)) +
                                            LRRR * (m2s * pap1 * papb - p1p2 * papb * (2.*pap2 + pbp2) +
                                                    pbp2 * (pap2 * pbp1 - pap1 * pbp2))) * RR - 2.*dd26 * pap2 * papb *
                                    (LL * p1p2 * RLLL + LL * m1 * m2 * RLRL + LRLR * m1 * m2 * RR + LRRR * p1p2 * RR)));

    dd->SetIeps(ieps + 2);
    dd->GetNpf(la, 2);
    dd23 = la[2];
    dd24 = la[3];
    dd25 = la[4];
    dd26 = la[5];
    dd27 = la[6];
    dd->GetNpf(la, 1);
    dd11 = la[0];
    dd12 = la[1];
    dd13 = la[2];
    delete dd;

    std::complex<double> me2 = -8.*ivs3v2 *
                               (dd27 * (4.*LL * p1p2 * papb * RLLL + 4.*LL * pap1 * pbp2 * RLLL +
                                        6.*LL * m1 * m2 * papb * RLRL + 6.*LRLR * m1 * m2 * papb * RR + 4.*LRRR * p1p2 * papb * RR +
                                        4.*LRRR * pap1 * pbp2 * RR) +
                                dd23 * m2s * (LL * (p1p2 * papb * RLLL - pap2 * pbp1 * RLLL + pap1 * pbp2 * RLLL +
                                        m1 * m2 * papb * RLRL) + (LRLR * m1 * m2 * papb +
                                                LRRR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) * RR) +
                                2.*papb * (LL * ((dd13 + dd25) * m2s * pbp1 + (dd12 + dd24) *
                                           (p1p2 * papb - pap2 * pbp1) + (-2.*(dd13 + dd25) * p1p2 + (dd12 + dd24) * pap1) *
                                           pbp2) * RLLL + LL * m1 * m2 * ((dd12 + dd24) * papb - (dd13 + dd25) * pbp2) * RLRL +
                                           ((dd13 + dd25) * (LRRR * m2s * pbp1 - LRLR * m1 * m2 * pbp2 - 2.*LRRR * p1p2 * pbp2) +
                                            dd24 * (LRLR * m1 * m2 * papb + LRRR * p1p2 * papb - LRRR * pap2 * pbp1 +
                                                    LRRR * pap1 * pbp2) + dd12 * (LRLR * m1 * m2 * papb +
                                                            LRRR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2))) * RR -
                                           dd26 * (LL * m2s * pap1 * RLLL + LL * m1 * m2 * pap2 * RLRL + LRRR * m2s * pap1 * RR +
                                                   LRLR * m1 * m2 * pap2 * RR)));

    return std::real(me0 + me1 + me2);
}

double FI::MVubo1t(int ieps)
{
    Npf * dd = new Npf(0., 0., m2s, m1s, 2.*papb, m2s - 2.*pap2,
                       m1l, m2l, m3l, m4l, mul, ieps);
    std::complex<double> la[7];
    dd->GetNpf(la, 2);
    std::complex<double> dd23 = la[2];
    std::complex<double> dd24 = la[3];
    std::complex<double> dd25 = la[4];
    std::complex<double> dd26 = la[5];
    std::complex<double> dd27 = la[6];
    dd->GetNpf(la, 1);
    std::complex<double> dd11 = la[0];
    std::complex<double> dd12 = la[1];
    std::complex<double> dd13 = la[2];

    std::complex<double> me0 = -4.*ivt1s2 *
                               (2.*dd27 * (-2.*LLLL * LR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                                           m1 * m2 * papb * (LL * RLRL + LRLR * RR) - 2.*(p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RL *
                                           RRRR) - dd23 * (LLLL * LR * m2s * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) -
                                                   2.*m1 * m2 * pap2 * pbp2 * (LL * RLRL + LRLR * RR) +
                                                   m2s * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RL * RRRR) +
                                2.*(-(dd11 * LLLL * LR * p1p2 * std::pow(papb, 2)) -
                                    dd24 * LLLL * LR * p1p2 * std::pow(papb, 2) +
                                    dd11 * LLLL * LR * pap2 * papb * pbp1 + dd24 * LLLL * LR * pap2 * papb * pbp1 +
                                    dd25 * LLLL * LR * p1p2 * papb * pbp2 + dd11 * LLLL * LR * pap1 * papb * pbp2 -
                                    2.*dd12 * LLLL * LR * pap1 * papb * pbp2 - dd24 * LLLL * LR * pap1 * papb * pbp2 -
                                    dd25 * LLLL * LR * pap2 * pbp1 * pbp2 + dd25 * LLLL * LR * pap1 * std::pow(pbp2, 2) +
                                    dd26 * LLLL * LR * pap2 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                                    dd13 * LLLL * LR * (-(m2s * pap1 * papb) + (pap2 + pbp2) * (p1p2 * papb - pap2 * pbp1 +
                                            pap1 * pbp2)) - dd13 * LL * m1 * m2 * pap2 * papb * RLRL -
                                    dd26 * LL * m1 * m2 * pap2 * papb * RLRL + dd11 * LL * m1 * m2 * std::pow(papb, 2) * RLRL +
                                    dd24 * LL * m1 * m2 * std::pow(papb, 2) * RLRL - dd25 * LL * m1 * m2 * papb * pbp2 * RLRL -
                                    dd13 * LRLR * m1 * m2 * pap2 * papb * RR - dd26 * LRLR * m1 * m2 * pap2 * papb * RR +
                                    dd11 * LRLR * m1 * m2 * std::pow(papb, 2) * RR + dd24 * LRLR * m1 * m2 * std::pow(papb, 2) * RR -
                                    dd25 * LRLR * m1 * m2 * papb * pbp2 * RR +
                                    ((dd26 * pap2 - (dd11 + dd24) * papb) * (p1p2 * papb - pap2 * pbp1) +
                                     (dd26 * pap1 * pap2 + dd25 * p1p2 * papb + (dd11 - 2.*dd12 - dd24) * pap1 * papb -
                                      dd25 * pap2 * pbp1) * pbp2 + dd25 * pap1 * std::pow(pbp2, 2) +
                                     dd13 * (-(m2s * pap1 * papb) + (pap2 + pbp2) * (p1p2 * papb - pap2 * pbp1 +
                                             pap1 * pbp2))) * RL * RRRR));

    dd->SetIeps(ieps + 1);
    dd->GetNpf(la, 2);
    dd23 = la[2];
    dd24 = la[3];
    dd25 = la[4];
    dd26 = la[5];
    dd27 = la[6];
    dd->GetNpf(la, 1);
    dd11 = la[0];
    dd12 = la[1];
    dd13 = la[2];

    std::complex<double> me1 = -4.*ivt1s2 * (-2.*(dd26 * m2s * pap1 * papb +
                               (((dd13 + dd25) * p1p2 - 2.*(dd12 + dd24) * pap1) * papb -
                                (dd13 + dd25) * pap2 * pbp1) * pbp2 + (dd13 + dd25) * pap1 * std::pow(pbp2, 2)) *
                               (LLLL * LR + RL * RRRR) + dd23 * (LLLL * LR * m2s * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                                       m1 * m2 * (m2s * papb - 2.*pap2 * pbp2) * (LL * RLRL + LRLR * RR) +
                                       m2s * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RL * RRRR) +
                               2.*dd27 * (LLLL * LR * (p1p2 * papb - pap2 * pbp1 + 3.*pap1 * pbp2) +
                                          m1 * m2 * papb * (LL * RLRL + LRLR * RR) + (p1p2 * papb - pap2 * pbp1 + 3.*pap1 * pbp2) * RL *
                                          RRRR));

    dd->SetIeps(ieps + 2);
    dd->GetNpf(la, 2);
    dd27 = la[6];
    delete dd;

    std::complex<double> me2 = 8.*dd27 * ivt1s2 *
                               (LLLL * LR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                                LL * m1 * m2 * papb * RLRL + LRLR * m1 * m2 * papb * RR - p1p2 * papb * RL * RRRR +
                                pap2 * pbp1 * RL * RRRR + pap1 * pbp2 * RL * RRRR);

    return std::real(me0 + me1 + me2);
}

double FI::MVubo1u(int ieps)
{
    Npf * dd = new Npf(0., 0., m2s, m1s, 2.*papb, m2s - 2.*pap2,
                       m1l, m2l, m3l, m4l, mul, ieps);
    std::complex<double> la[7];
    dd->GetNpf(la, 2);
    std::complex<double> dd23 = la[2];
    std::complex<double> dd24 = la[3];
    std::complex<double> dd25 = la[4];
    std::complex<double> dd26 = la[5];
    std::complex<double> dd27 = la[6];
    dd->GetNpf(la, 1);
    std::complex<double> dd11 = la[0];
    std::complex<double> dd12 = la[1];
    std::complex<double> dd13 = la[2];

    std::complex<double> me0 = 4.*ivu2s2 *
                               (dd23 * m2s * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) *
                                (LLLL * LR + LL * RLRL + LRLR * RR + RL * RRRR) +
                                4.*dd27 * (LLLL * LR * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) + LL * pap2 * pbp1 * RLRL +
                                           LRLR * pap2 * pbp1 * RR + (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RL * RRRR) -
                                2.*(-(dd12 * LLLL * LR * p1p2 * std::pow(papb, 2)) -
                                    dd24 * LLLL * LR * p1p2 * std::pow(papb, 2) -
                                    2.*dd11 * LLLL * LR * pap2 * papb * pbp1 + dd12 * LLLL * LR * pap2 * papb * pbp1 -
                                    dd24 * LLLL * LR * pap2 * papb * pbp1 + dd25 * LLLL * LR * p1p2 * papb * pbp2 +
                                    dd12 * LLLL * LR * pap1 * papb * pbp2 + dd24 * LLLL * LR * pap1 * papb * pbp2 +
                                    dd25 * LLLL * LR * pap2 * pbp1 * pbp2 - dd25 * LLLL * LR * pap1 * std::pow(pbp2, 2) +
                                    dd25 * LL * m2s * papb * pbp1 * RLRL - 2.*dd11 * LL * pap2 * papb * pbp1 * RLRL -
                                    2.*dd24 * LL * pap2 * papb * pbp1 * RLRL + dd25 * LRLR * m2s * papb * pbp1 * RR -
                                    2.*dd11 * LRLR * pap2 * papb * pbp1 * RR - 2.*dd24 * LRLR * pap2 * papb * pbp1 * RR -
                                    (papb * ((dd12 + dd24) * p1p2 * papb + (2.*dd11 - dd12 + dd24) * pap2 * pbp1) -
                                     ((dd12 + dd24) * pap1 * papb + dd25 * (p1p2 * papb + pap2 * pbp1)) * pbp2 +
                                     dd25 * pap1 * std::pow(pbp2, 2)) * RL * RRRR + dd26 * pap2 * (p1p2 * papb + pap2 * pbp1 -
                                             pap1 * pbp2) * (LLLL * LR + LL * RLRL + LRLR * RR + RL * RRRR) +
                                    dd13 * (LLLL * LR * ((std::pow(pap2, 2) - m2s * papb) * pbp1 +
                                            pap2 * (-pap1 + pbp1) * pbp2 - pap1 * std::pow(pbp2, 2) +
                                            p1p2 * papb * (pap2 + pbp2)) +
                                            pap2 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LL * RLRL + LRLR * RR) +
                                            ((std::pow(pap2, 2) - m2s * papb) * pbp1 + pap2 * (-pap1 + pbp1) * pbp2 -
                                                    pap1 * std::pow(pbp2, 2) + p1p2 * papb * (pap2 + pbp2)) * RL * RRRR)));

    dd->SetIeps(ieps + 1);
    dd->GetNpf(la, 2);
    dd23 = la[2];
    dd24 = la[3];
    dd25 = la[4];
    dd26 = la[5];
    dd27 = la[6];
    dd->GetNpf(la, 1);
    dd11 = la[0];
    dd12 = la[1];
    dd13 = la[2];

    std::complex<double> me1 = -8.*ivu2s2 *
                               ((papb * ((dd13 + dd25) * m2s * pbp1 + (dd12 + dd24) * (p1p2 * papb - pap2 * pbp1)) -
                                 (((dd13 + dd25) * p1p2 + (dd12 + dd24) * pap1) * papb +
                                  (dd13 + dd25) * pap2 * pbp1) * pbp2 + (dd13 + dd25) * pap1 * std::pow(pbp2, 2)) *
                                (LLLL * LR - LL * RLRL - LRLR * RR + RL * RRRR) +
                                dd27 * (LLLL * LR * (3.*p1p2 * papb + pap2 * pbp1 - 3.*pap1 * pbp2) -
                                        (3.*p1p2 * papb - pap2 * pbp1 - 3.*pap1 * pbp2) * (LL * RLRL + LRLR * RR) +
                                        (3.*p1p2 * papb + pap2 * pbp1 - 3.*pap1 * pbp2) * RL * RRRR));

    dd->SetIeps(ieps + 2);
    dd->GetNpf(la, 2);
    dd27 = la[6];
    delete dd;

    std::complex<double> me2 = 16.*dd27 * ivu2s2 *
                               (p1p2 * papb - pap1 * pbp2) * (LLLL * LR - LL * RLRL - LRLR * RR + RL * RRRR);

    return std::real(me0 + me1 + me2);
}

double FI::MVtbo2s(int ieps)
{
    Npf * dd = new Npf(0., 0., m1s, m2s, 2.*papb, m1s - 2.*pap1,
                       m1l, m2l, m3l, m4l, mul, ieps);
    std::complex<double> la[7];
    dd->GetNpf(la, 2);
    std::complex<double> dd23 = la[2];
    std::complex<double> dd25 = la[4];
    std::complex<double> dd26 = la[5];
    std::complex<double> dd27 = la[6];
    dd->GetNpf(la, 1);
    std::complex<double> dd13 = la[2];
    delete dd;

    std::complex<double> me0 = 4.*ivs3v2 *
                               (dd23 * LR * LRLR * m1 * m1s * m2 * papb - dd23 * LR * LRRR * m1s * p1p2 * papb +
                                2.*dd26 * LR * LRRR * m1s * pap2 * papb - 2.*dd23 * LR * LRLR * m1 * m2 * pap1 * pbp1 +
                                dd23 * LR * LRRR * m1s * pap2 * pbp1 - 4.*dd26 * LR * LRRR * pap1 * pap2 * pbp1 +
                                2.*dd25 * LR * LRRR * p1p2 * papb * pbp1 - 2.*dd25 * LR * LRRR * pap2 * std::pow(pbp1, 2) +
                                dd23 * LR * LRRR * m1s * pap1 * pbp2 - 2.*dd25 * LR * LRRR * pap1 * pbp1 * pbp2 +
                                dd23 * LRLL * m1 * m1s * m2 * papb * RL - dd23 * LRRL * m1s * p1p2 * papb * RL +
                                2.*dd26 * LRRL * p1p2 * pap1 * papb * RL - 2.*dd23 * LRLL * m1 * m2 * pap1 * pbp1 * RL +
                                dd23 * LRRL * m1s * pap2 * pbp1 * RL - 2.*dd26 * LRRL * pap1 * pap2 * pbp1 * RL +
                                dd23 * LRRL * m1s * pap1 * pbp2 * RL - 2.*dd26 * LRRL * std::pow(pap1, 2) * pbp2 * RL +
                                2.*dd25 * LRRL * m1s * papb * pbp2 * RL - 4.*dd25 * LRRL * pap1 * pbp1 * pbp2 * RL -
                                dd23 * m1s * p1p2 * papb * RL * RLLL + 2.*dd26 * m1s * pap2 * papb * RL * RLLL +
                                dd23 * m1s * pap2 * pbp1 * RL * RLLL - 4.*dd26 * pap1 * pap2 * pbp1 * RL * RLLL +
                                2.*dd25 * p1p2 * papb * pbp1 * RL * RLLL - 2.*dd25 * pap2 * std::pow(pbp1, 2) * RL * RLLL +
                                dd23 * m1s * pap1 * pbp2 * RL * RLLL - 2.*dd25 * pap1 * pbp1 * pbp2 * RL * RLLL -
                                dd23 * LR * m1s * p1p2 * papb * RLLR + 2.*dd26 * LR * p1p2 * pap1 * papb * RLLR +
                                dd23 * LR * m1s * pap2 * pbp1 * RLLR - 2.*dd26 * LR * pap1 * pap2 * pbp1 * RLLR +
                                dd23 * LR * m1s * pap1 * pbp2 * RLLR - 2.*dd26 * LR * std::pow(pap1, 2) * pbp2 * RLLR +
                                2.*dd25 * LR * m1s * papb * pbp2 * RLLR - 4.*dd25 * LR * pap1 * pbp1 * pbp2 * RLLR +
                                dd23 * m1 * m1s * m2 * papb * RL * RLRL - 2.*dd23 * m1 * m2 * pap1 * pbp1 * RL * RLRL +
                                dd23 * LR * m1 * m2 * (m1s * papb - 2.*pap1 * pbp1) * RLRR +
                                2.*dd27 * (RL * (2.*(LRRL * pap1 * pbp2 + pap2 * pbp1 * RLLL) +
                                           m1 * m2 * papb * (LRLL + RLRL)) + LR * (2.*(LRRR * pap2 * pbp1 + pap1 * pbp2 * RLLR) +
                                                   m1 * m2 * papb * (LRLR + RLRR))) +
                                dd13 * (RL * (-(LRRL * m1s * p1p2 * papb) - LLLL * m1 * ml4 * p1p2 * papb +
                                        2.*LRRL * p1p2 * pap1 * papb + LRRL * m1s * pap2 * pbp1 + LLLL * m1 * ml4 * pap2 * pbp1 -
                                        2.*LRRL * pap1 * pap2 * pbp1 + LRLL * m1 * m2 * (m1s * papb - 2.*pap1 * pbp1) +
                                        LLRL * m2 * ml4 * (m1s * papb - 2.*pap1 * pbp1) + LRRL * m1s * pap1 * pbp2 +
                                        LLLL * m1 * ml4 * pap1 * pbp2 - 2.*LRRL * std::pow(pap1, 2) * pbp2 +
                                        2.*LRRL * m1s * papb * pbp2 -
                                        4.*LRRL * pap1 * pbp1 * pbp2 - m1s * p1p2 * papb * RLLL + 2.*m1s * pap2 * papb * RLLL +
                                        m1s * pap2 * pbp1 * RLLL - 4.*pap1 * pap2 * pbp1 * RLLL + 2.*p1p2 * papb * pbp1 * RLLL -
                                        2.*pap2 * std::pow(pbp1, 2) * RLLL + m1s * pap1 * pbp2 * RLLL -
                                        2.*pap1 * pbp1 * pbp2 * RLLL +
                                        m1 * m1s * m2 * papb * RLRL - 2.*m1 * m2 * pap1 * pbp1 * RLRL + m1s * m2 * ml4 * papb * RRLL -
                                        2.*m2 * ml4 * pap1 * pbp1 * RRLL - m1 * ml4 * p1p2 * papb * RRRL +
                                        m1 * ml4 * pap2 * pbp1 * RRRL + m1 * ml4 * pap1 * pbp2 * RRRL) +
                                        LR * (-(LRRR * m1s * p1p2 * papb) - LLLR * m1 * ml4 * p1p2 * papb +
                                                2.*LRRR * m1s * pap2 * papb + LRRR * m1s * pap2 * pbp1 + LLLR * m1 * ml4 * pap2 * pbp1 -
                                                4.*LRRR * pap1 * pap2 * pbp1 + 2.*LRRR * p1p2 * papb * pbp1 -
                                                2.*LRRR * pap2 * std::pow(pbp1, 2) +
                                                LRLR * m1 * m2 * (m1s * papb - 2.*pap1 * pbp1) + LLRR * m2 * ml4 *
                                                (m1s * papb - 2.*pap1 * pbp1) + LRRR * m1s * pap1 * pbp2 +
                                                LLLR * m1 * ml4 * pap1 * pbp2 - 2.*LRRR * pap1 * pbp1 * pbp2 - m1s * p1p2 * papb * RLLR +
                                                2.*p1p2 * pap1 * papb * RLLR + m1s * pap2 * pbp1 * RLLR - 2.*pap1 * pap2 * pbp1 * RLLR +
                                                m1s * pap1 * pbp2 * RLLR - 2.*std::pow(pap1, 2) * pbp2 * RLLR +
                                                2.*m1s * papb * pbp2 * RLLR -
                                                4.*pap1 * pbp1 * pbp2 * RLLR + m1 * m1s * m2 * papb * RLRR - 2.*m1 * m2 * pap1 * pbp1 * RLRR +
                                                m1s * m2 * ml4 * papb * RRLR - 2.*m2 * ml4 * pap1 * pbp1 * RRLR -
                                                m1 * ml4 * p1p2 * papb * RRRR + m1 * ml4 * pap2 * pbp1 * RRRR + m1 * ml4 * pap1 * pbp2 * RRRR)));

    return std::real(me0);
}

double FI::MVtbo2t(int ieps)
{
    Npf * dd = new Npf(0., 0., m1s, m2s, 2.*papb, m1s - 2.*pap1,
                       m1l, m2l, m3l, m4l, mul, ieps);
    std::complex<double> la[7];

    dd->GetNpf(la, 2);
    std::complex<double> dd23 = la[2];
    std::complex<double> dd25 = la[4];
    std::complex<double> dd26 = la[5];
    std::complex<double> dd27 = la[6];
    dd->GetNpf(la, 1);
    std::complex<double> dd11 = la[0];
    std::complex<double> dd12 = la[1];
    std::complex<double> dd13 = la[2];
    dd->GetNpf(la, 0);
    std::complex<double> dd0 = la[0];
    delete dd;

    std::complex<double> me0 = 2.*ivt1s2 *
                               (dd23 * LR * LRLR * m1 * m1s * m2 * papb - dd0 * LL * LRLL * m1s * m2 * ml2 * papb +
                                dd0 * LL * LLLL * ml2 * ml4 * p1p2 * papb + 2.*dd0 * LL * LRLL * m2 * ml2 * pap1 * papb +
                                2.*dd12 * LL * LRLL * m2 * ml2 * pap1 * papb - 2.*dd23 * LR * LRLR * m1 * m2 * pap1 * pbp1 -
                                dd0 * LL * LLLL * ml2 * ml4 * pap2 * pbp1 + dd0 * LL * LLLL * ml2 * ml4 * pap1 * pbp2 -
                                dd23 * LRRL * m1s * p1p2 * papb * RL + 2.*dd26 * LRRL * p1p2 * pap1 * papb * RL +
                                dd23 * LRRL * m1s * pap2 * pbp1 * RL - 2.*dd26 * LRRL * pap1 * pap2 * pbp1 * RL +
                                dd23 * LRRL * m1s * pap1 * pbp2 * RL - 2.*dd26 * LRRL * std::pow(pap1, 2) * pbp2 * RL +
                                2.*dd25 * LRRL * m1s * papb * pbp2 * RL - 4.*dd25 * LRRL * pap1 * pbp1 * pbp2 * RL +
                                dd0 * LL * m1 * ml2 * p1p2 * papb * RLLL - dd0 * LL * m1 * ml2 * pap2 * pbp1 * RLLL +
                                dd0 * LL * m1 * ml2 * pap1 * pbp2 * RLLL - 2.*dd0 * LL * m1 * ml2 * papb * pbp2 * RLLL -
                                2.*dd11 * LL * m1 * ml2 * papb * pbp2 * RLLL - dd23 * LR * m1s * p1p2 * papb * RLLR +
                                2.*dd26 * LR * p1p2 * pap1 * papb * RLLR + dd23 * LR * m1s * pap2 * pbp1 * RLLR -
                                2.*dd26 * LR * pap1 * pap2 * pbp1 * RLLR + dd23 * LR * m1s * pap1 * pbp2 * RLLR -
                                2.*dd26 * LR * std::pow(pap1, 2) * pbp2 * RLLR + 2.*dd25 * LR * m1s * papb * pbp2 * RLLR -
                                4.*dd25 * LR * pap1 * pbp1 * pbp2 * RLLR + dd23 * m1 * m1s * m2 * papb * RL * RLRL -
                                2.*dd23 * m1 * m2 * pap1 * pbp1 * RL * RLRL +
                                2.*dd27 * (LR * LRLR * m1 * m2 * papb + 2.*LRRL * pap1 * pbp2 * RL + 2.*LR * pap1 * pbp2 * RLLR +
                                           m1 * m2 * papb * RL * RLRL) - dd0 * LLRR * m1 * m2 * ml2 * ml4 * papb * RR +
                                dd0 * LRRR * m1 * ml2 * p1p2 * papb * RR - dd0 * LRRR * m1 * ml2 * pap2 * pbp1 * RR +
                                dd0 * LRRR * m1 * ml2 * pap1 * pbp2 * RR - 2.*dd0 * LRRR * m1 * ml2 * papb * pbp2 * RR -
                                2.*dd11 * LRRR * m1 * ml2 * papb * pbp2 * RR - dd0 * m1s * m2 * ml2 * papb * RLRR * RR +
                                2.*dd0 * m2 * ml2 * pap1 * papb * RLRR * RR + 2.*dd12 * m2 * ml2 * pap1 * papb * RLRR * RR -
                                dd0 * LL * m1 * m2 * ml2 * ml4 * papb * RRLL +
                                dd13 * (LLRL * m1s * m2 * ml4 * papb * RL - LRRL * m1s * p1p2 * papb * RL +
                                        2.*LRRL * p1p2 * pap1 * papb * RL - 2.*LLRL * m2 * ml4 * pap1 * pbp1 * RL +
                                        LRRL * m1s * pap2 * pbp1 * RL - 2.*LRRL * pap1 * pap2 * pbp1 * RL +
                                        LRRL * m1s * pap1 * pbp2 * RL - 2.*LRRL * std::pow(pap1, 2) * pbp2 * RL +
                                        2.*LRRL * m1s * papb * pbp2 * RL -
                                        4.*LRRL * pap1 * pbp1 * pbp2 * RL + LL * ml2 * (-(LRLL * m1s * m2 * papb) +
                                                m1 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RLLL) + m1 * m1s * m2 * papb * RL * RLRL -
                                        2.*m1 * m2 * pap1 * pbp1 * RL * RLRL + LRRR * m1 * ml2 * p1p2 * papb * RR -
                                        LRRR * m1 * ml2 * pap2 * pbp1 * RR + LRRR * m1 * ml2 * pap1 * pbp2 * RR -
                                        m1s * m2 * ml2 * papb * RLRR * RR + LR * (LRLR * m1 * m2 * (m1s * papb - 2.*pap1 * pbp1) +
                                                LLLR * m1 * ml4 * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) -
                                                ((m1s - 2.*pap1) * (p1p2 * papb - pap2 * pbp1) -
                                                        (m1s * (pap1 + 2.*papb) - 2.*pap1 * (pap1 + 2.*pbp1)) * pbp2) * RLLR +
                                                m2 * ml4 * (m1s * papb - 2.*pap1 * pbp1) * RRLR) -
                                        m1 * ml4 * (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) * RL * RRRL) +
                                dd0 * ml2 * ml4 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RR * RRRR);

    return std::real(me0);
}

double FI::MVtbo2u(int ieps)
{
    Npf * dd = new Npf(0., 0., m1s, m2s, 2.*papb, m1s - 2.*pap1,
                       m1l, m2l, m3l, m4l, mul, ieps);
    std::complex<double> la[7];

    dd->GetNpf(la, 2);
    std::complex<double> dd23 = la[2];
    std::complex<double> dd25 = la[4];
    std::complex<double> dd26 = la[5];
    std::complex<double> dd27 = la[6];
    dd->GetNpf(la, 1);
    std::complex<double> dd11 = la[0];
    std::complex<double> dd12 = la[1];
    std::complex<double> dd13 = la[2];
    dd->GetNpf(la, 0);
    std::complex<double> dd0 = la[0];
    delete dd;

    std::complex<double> me0 = -2.*ivu2s2 *
                               (-(dd23 * LR * LRLR * m1s * p1p2 * papb) + 2.*dd26 * LR * LRLR * m1s * pap2 * papb +
                                4.*dd27 * LR * LRLR * pap2 * pbp1 + dd23 * LR * LRLR * m1s * pap2 * pbp1 -
                                4.*dd26 * LR * LRLR * pap1 * pap2 * pbp1 - 2.*dd11 * LL * LRLL * m2 * ml2 * papb * pbp1 +
                                2.*dd25 * LR * LRLR * p1p2 * papb * pbp1 - 2.*dd25 * LR * LRLR * pap2 * std::pow(pbp1, 2) +
                                dd23 * LR * LRLR * m1s * pap1 * pbp2 - 2.*dd25 * LR * LRLR * pap1 * pbp1 * pbp2 +
                                2.*dd27 * LRRL * m1 * m2 * papb * RL + dd23 * LRRL * m1 * m1s * m2 * papb * RL -
                                2.*dd23 * LRRL * m1 * m2 * pap1 * pbp1 * RL + 2.*dd12 * LL * m1 * ml2 * pap2 * papb * RLLL +
                                2.*dd27 * LR * m1 * m2 * papb * RLLR + dd23 * LR * m1 * m1s * m2 * papb * RLLR -
                                2.*dd23 * LR * m1 * m2 * pap1 * pbp1 * RLLR - dd23 * m1s * p1p2 * papb * RL * RLRL +
                                2.*dd26 * m1s * pap2 * papb * RL * RLRL + 4.*dd27 * pap2 * pbp1 * RL * RLRL +
                                dd23 * m1s * pap2 * pbp1 * RL * RLRL - 4.*dd26 * pap1 * pap2 * pbp1 * RL * RLRL +
                                2.*dd25 * p1p2 * papb * pbp1 * RL * RLRL - 2.*dd25 * pap2 * std::pow(pbp1, 2) * RL * RLRL +
                                dd23 * m1s * pap1 * pbp2 * RL * RLRL - 2.*dd25 * pap1 * pbp1 * pbp2 * RL * RLRL +
                                2.*ml2 * papb * (dd12 * LRRR * m1 * pap2 - dd11 * m2 * pbp1 * RLRR) * RR +
                                dd13 * (-(LR * LRLR * m1s * p1p2 * papb) + 2.*LR * LRLR * m1s * pap2 * papb +
                                        LR * LRLR * m1s * pap2 * pbp1 - 4.*LR * LRLR * pap1 * pap2 * pbp1 +
                                        2.*LR * LRLR * p1p2 * papb * pbp1 - 2.*LR * LRLR * pap2 * std::pow(pbp1, 2) +
                                        LLLR * LR * m2 * ml4 * (m1s * papb - 2.*pap1 * pbp1) + LR * LRLR * m1s * pap1 * pbp2 -
                                        2.*LR * LRLR * pap1 * pbp1 * pbp2 + LRRL * m1 * m1s * m2 * papb * RL -
                                        LLRL * m1 * ml4 * p1p2 * papb * RL - 2.*LRRL * m1 * m2 * pap1 * pbp1 * RL +
                                        LLRL * m1 * ml4 * pap2 * pbp1 * RL + LLRL * m1 * ml4 * pap1 * pbp2 * RL +
                                        LL * ml2 * (LRLL * m1s * m2 * papb - m1 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RLLL) +
                                        LR * m1 * m1s * m2 * papb * RLLR - 2.*LR * m1 * m2 * pap1 * pbp1 * RLLR -
                                        m1s * p1p2 * papb * RL * RLRL + 2.*m1s * pap2 * papb * RL * RLRL + m1s * pap2 * pbp1 * RL * RLRL -
                                        4.*pap1 * pap2 * pbp1 * RL * RLRL + 2.*p1p2 * papb * pbp1 * RL * RLRL -
                                        2.*pap2 * std::pow(pbp1, 2) * RL * RLRL + m1s * pap1 * pbp2 * RL * RLRL -
                                        2.*pap1 * pbp1 * pbp2 * RL * RLRL - LRRR * m1 * ml2 * p1p2 * papb * RR -
                                        LRRR * m1 * ml2 * pap2 * pbp1 * RR + LRRR * m1 * ml2 * pap1 * pbp2 * RR +
                                        m1s * m2 * ml2 * papb * RLRR * RR - LR * m1 * ml4 * p1p2 * papb * RRLR +
                                        LR * m1 * ml4 * pap2 * pbp1 * RRLR + LR * m1 * ml4 * pap1 * pbp2 * RRLR +
                                        m2 * ml4 * (m1s * papb - 2.*pap1 * pbp1) * RL * RRRL) +
                                dd0 * ml2 * (LL * (LRLL * m2 * papb * (m1s - 2.*pbp1) -
                                        LLLL * ml4 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) +
                                        m1 * (-(p1p2 * papb * RLLL) + 2.*pap2 * papb * RLLL - pap2 * pbp1 * RLLL +
                                                pap1 * pbp2 * RLLL + m2 * ml4 * papb * RRLL)) +
                                        RR * (LLRR * m1 * m2 * ml4 * papb + LRRR * m1 * (-(p1p2 * papb) + 2.*pap2 * papb -
                                                pap2 * pbp1 + pap1 * pbp2) + m2 * papb * (m1s - 2.*pbp1) * RLRR -
                                                ml4 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RRRR)));

    return std::real(me0);
}

double FI::MVubo2s(int ieps)
{
    Npf * dd = new Npf(0., 0., m2s, m1s, 2.*papb, m2s - 2.*pap2,
                       m1l, m2l, m3l, m4l, mul, ieps);
    std::complex<double> la[7];

    dd->GetNpf(la, 2);
    std::complex<double> dd23 = la[2];
    std::complex<double> dd25 = la[4];
    std::complex<double> dd26 = la[5];
    std::complex<double> dd27 = la[6];
    dd->GetNpf(la, 1);
    std::complex<double> dd13 = la[2];
    delete dd;

    std::complex<double> me0 = -4.*ivs3v2 *
                               (dd23 * LR * LRRR * m1 * m2 * m2s * papb - dd23 * LR * LRLR * m2s * p1p2 * papb +
                                2.*dd26 * LR * LRLR * m2s * pap1 * papb + dd23 * LR * LRLR * m2s * pap2 * pbp1 +
                                dd23 * LR * LRLR * m2s * pap1 * pbp2 - 2.*dd23 * LR * LRRR * m1 * m2 * pap2 * pbp2 -
                                4.*dd26 * LR * LRLR * pap1 * pap2 * pbp2 + 2.*dd25 * LR * LRLR * p1p2 * papb * pbp2 -
                                2.*dd25 * LR * LRLR * pap2 * pbp1 * pbp2 - 2.*dd25 * LR * LRLR * pap1 * std::pow(pbp2, 2) +
                                dd23 * LRRL * m1 * m2 * m2s * papb * RL - dd23 * LRLL * m2s * p1p2 * papb * RL +
                                2.*dd26 * LRLL * p1p2 * pap2 * papb * RL + dd23 * LRLL * m2s * pap2 * pbp1 * RL -
                                2.*dd26 * LRLL * std::pow(pap2, 2) * pbp1 * RL + 2.*dd25 * LRLL * m2s * papb * pbp1 * RL +
                                dd23 * LRLL * m2s * pap1 * pbp2 * RL - 2.*dd23 * LRRL * m1 * m2 * pap2 * pbp2 * RL -
                                2.*dd26 * LRLL * pap1 * pap2 * pbp2 * RL - 4.*dd25 * LRLL * pap2 * pbp1 * pbp2 * RL +
                                dd23 * m1 * m2 * m2s * papb * RL * RLLL - 2.*dd23 * m1 * m2 * pap2 * pbp2 * RL * RLLL +
                                dd23 * LR * m1 * m2 * m2s * papb * RLLR - 2.*dd23 * LR * m1 * m2 * pap2 * pbp2 * RLLR -
                                dd23 * m2s * p1p2 * papb * RL * RLRL + 2.*dd26 * m2s * pap1 * papb * RL * RLRL +
                                dd23 * m2s * pap2 * pbp1 * RL * RLRL + dd23 * m2s * pap1 * pbp2 * RL * RLRL -
                                4.*dd26 * pap1 * pap2 * pbp2 * RL * RLRL + 2.*dd25 * p1p2 * papb * pbp2 * RL * RLRL -
                                2.*dd25 * pap2 * pbp1 * pbp2 * RL * RLRL - 2.*dd25 * pap1 * std::pow(pbp2, 2) * RL * RLRL +
                                LR * (dd23 * m2s * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                                      2.*(-(dd26 * pap2 * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
                                          dd25 * pbp1 * (m2s * papb - 2.*pap2 * pbp2))) * RLRR +
                                2.*dd27 * (RL * (m1 * m2 * papb * (LRRL + RLLL) +
                                           2.*(LRLL * pap2 * pbp1 + pap1 * pbp2 * RLRL)) +
                                           LR * (m1 * m2 * papb * (LRRR + RLLR) + 2.*(LRLR * pap1 * pbp2 + pap2 * pbp1 * RLRR))) +
                                dd13 * (RL * (-(LRLL * m2s * p1p2 * papb) - LLRL * m2 * ml4 * p1p2 * papb +
                                        2.*LRLL * p1p2 * pap2 * papb + LRLL * m2s * pap2 * pbp1 + LLRL * m2 * ml4 * pap2 * pbp1 -
                                        2.*LRLL * std::pow(pap2, 2) * pbp1 + 2.*LRLL * m2s * papb * pbp1 +
                                        LRLL * m2s * pap1 * pbp2 +
                                        LLRL * m2 * ml4 * pap1 * pbp2 - 2.*LRLL * pap1 * pap2 * pbp2 - 4.*LRLL * pap2 * pbp1 * pbp2 +
                                        LRRL * m1 * m2 * (m2s * papb - 2.*pap2 * pbp2) + LLLL * m1 * ml4 *
                                        (m2s * papb - 2.*pap2 * pbp2) + m1 * m2 * m2s * papb * RLLL -
                                        2.*m1 * m2 * pap2 * pbp2 * RLLL - m2s * p1p2 * papb * RLRL + 2.*m2s * pap1 * papb * RLRL +
                                        m2s * pap2 * pbp1 * RLRL + m2s * pap1 * pbp2 * RLRL - 4.*pap1 * pap2 * pbp2 * RLRL +
                                        2.*p1p2 * papb * pbp2 * RLRL - 2.*pap2 * pbp1 * pbp2 * RLRL -
                                        2.*pap1 * std::pow(pbp2, 2) * RLRL -
                                        m2 * ml4 * p1p2 * papb * RRLL + m2 * ml4 * pap2 * pbp1 * RRLL + m2 * ml4 * pap1 * pbp2 * RRLL +
                                        m1 * ml4 * (m2s * papb - 2.*pap2 * pbp2) * RRRL) +
                                        LR * (-(LRLR * m2s * p1p2 * papb) - LLRR * m2 * ml4 * p1p2 * papb +
                                                2.*LRLR * m2s * pap1 * papb + LRLR * m2s * pap2 * pbp1 + LLRR * m2 * ml4 * pap2 * pbp1 +
                                                LRLR * m2s * pap1 * pbp2 + LLRR * m2 * ml4 * pap1 * pbp2 - 4.*LRLR * pap1 * pap2 * pbp2 +
                                                2.*LRLR * p1p2 * papb * pbp2 - 2.*LRLR * pap2 * pbp1 * pbp2 -
                                                2.*LRLR * pap1 * std::pow(pbp2, 2) + LRRR * m1 * m2 * (m2s * papb - 2.*pap2 * pbp2) +
                                                LLLR * m1 * ml4 * (m2s * papb - 2.*pap2 * pbp2) + m1 * m2 * m2s * papb * RLLR -
                                                2.*m1 * m2 * pap2 * pbp2 * RLLR - m2s * p1p2 * papb * RLRR + 2.*p1p2 * pap2 * papb * RLRR +
                                                m2s * pap2 * pbp1 * RLRR - 2.*std::pow(pap2, 2) * pbp1 * RLRR + 2.*m2s * papb * pbp1 * RLRR +
                                                m2s * pap1 * pbp2 * RLRR - 2.*pap1 * pap2 * pbp2 * RLRR - 4.*pap2 * pbp1 * pbp2 * RLRR -
                                                m2 * ml4 * p1p2 * papb * RRLR + m2 * ml4 * pap2 * pbp1 * RRLR + m2 * ml4 * pap1 * pbp2 * RRLR +
                                                m1 * ml4 * (m2s * papb - 2.*pap2 * pbp2) * RRRR)));

    return std::real(me0);
}

double FI::MVubo2t(int ieps)
{
    Npf * dd = new Npf(0., 0., m2s, m1s, 2.*papb, m2s - 2.*pap2,
                       m1l, m2l, m3l, m4l, mul, ieps);
    std::complex<double> la[7];

    dd->GetNpf(la, 2);
    std::complex<double> dd23 = la[2];
    std::complex<double> dd25 = la[4];
    std::complex<double> dd26 = la[5];
    std::complex<double> dd27 = la[6];
    dd->GetNpf(la, 1);
    std::complex<double> dd11 = la[0];
    std::complex<double> dd12 = la[1];
    std::complex<double> dd13 = la[2];
    dd->GetNpf(la, 0);
    std::complex<double> dd0 = la[0];
    delete dd;

    std::complex<double> me0 = -2.*ivt1s2 *
                               (-(dd23 * LR * LRLR * m2s * p1p2 * papb) + 2.*dd26 * LR * LRLR * m2s * pap1 * papb +
                                dd23 * LR * LRLR * m2s * pap2 * pbp1 + 4.*dd27 * LR * LRLR * pap1 * pbp2 +
                                dd23 * LR * LRLR * m2s * pap1 * pbp2 - 4.*dd26 * LR * LRLR * pap1 * pap2 * pbp2 -
                                2.*dd11 * LL * LRLL * m1 * ml2 * papb * pbp2 + 2.*dd25 * LR * LRLR * p1p2 * papb * pbp2 -
                                2.*dd25 * LR * LRLR * pap2 * pbp1 * pbp2 - 2.*dd25 * LR * LRLR * pap1 * std::pow(pbp2, 2) +
                                2.*dd27 * LRRL * m1 * m2 * papb * RL + dd23 * LRRL * m1 * m2 * m2s * papb * RL -
                                2.*dd23 * LRRL * m1 * m2 * pap2 * pbp2 * RL + 2.*dd12 * LL * m2 * ml2 * pap1 * papb * RLLL +
                                2.*dd27 * LR * m1 * m2 * papb * RLLR + dd23 * LR * m1 * m2 * m2s * papb * RLLR -
                                2.*dd23 * LR * m1 * m2 * pap2 * pbp2 * RLLR - dd23 * m2s * p1p2 * papb * RL * RLRL +
                                2.*dd26 * m2s * pap1 * papb * RL * RLRL + dd23 * m2s * pap2 * pbp1 * RL * RLRL +
                                4.*dd27 * pap1 * pbp2 * RL * RLRL + dd23 * m2s * pap1 * pbp2 * RL * RLRL -
                                4.*dd26 * pap1 * pap2 * pbp2 * RL * RLRL + 2.*dd25 * p1p2 * papb * pbp2 * RL * RLRL -
                                2.*dd25 * pap2 * pbp1 * pbp2 * RL * RLRL - 2.*dd25 * pap1 * std::pow(pbp2, 2) * RL * RLRL +
                                2.*ml2 * papb * (dd12 * LRRR * m2 * pap1 - dd11 * m1 * pbp2 * RLRR) * RR +
                                dd13 * (-(LR * LRLR * m2s * p1p2 * papb) + 2.*LR * LRLR * m2s * pap1 * papb +
                                        LR * LRLR * m2s * pap2 * pbp1 + LR * LRLR * m2s * pap1 * pbp2 -
                                        4.*LR * LRLR * pap1 * pap2 * pbp2 + 2.*LR * LRLR * p1p2 * papb * pbp2 -
                                        2.*LR * LRLR * pap2 * pbp1 * pbp2 - 2.*LR * LRLR * pap1 * std::pow(pbp2, 2) +
                                        LLLR * LR * m1 * ml4 * (m2s * papb - 2.*pap2 * pbp2) + LRRL * m1 * m2 * m2s * papb * RL -
                                        LLRL * m2 * ml4 * p1p2 * papb * RL + LLRL * m2 * ml4 * pap2 * pbp1 * RL +
                                        LLRL * m2 * ml4 * pap1 * pbp2 * RL - 2.*LRRL * m1 * m2 * pap2 * pbp2 * RL +
                                        LL * ml2 * (LRLL * m1 * m2s * papb - m2 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RLLL) +
                                        LR * m1 * m2 * m2s * papb * RLLR - 2.*LR * m1 * m2 * pap2 * pbp2 * RLLR -
                                        m2s * p1p2 * papb * RL * RLRL + 2.*m2s * pap1 * papb * RL * RLRL + m2s * pap2 * pbp1 * RL * RLRL +
                                        m2s * pap1 * pbp2 * RL * RLRL - 4.*pap1 * pap2 * pbp2 * RL * RLRL +
                                        2.*p1p2 * papb * pbp2 * RL * RLRL - 2.*pap2 * pbp1 * pbp2 * RL * RLRL -
                                        2.*pap1 * std::pow(pbp2, 2) * RL * RLRL - LRRR * m2 * ml2 * p1p2 * papb * RR +
                                        LRRR * m2 * ml2 * pap2 * pbp1 * RR - LRRR * m2 * ml2 * pap1 * pbp2 * RR +
                                        m1 * m2s * ml2 * papb * RLRR * RR - LR * m2 * ml4 * p1p2 * papb * RRLR +
                                        LR * m2 * ml4 * pap2 * pbp1 * RRLR + LR * m2 * ml4 * pap1 * pbp2 * RRLR +
                                        m1 * ml4 * (m2s * papb - 2.*pap2 * pbp2) * RL * RRRL) +
                                dd0 * ml2 * (LL * (LRLL * m1 * papb * (m2s - 2.*pbp2) -
                                        LLLL * ml4 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                                        m2 * (-(p1p2 * papb * RLLL) + 2.*pap1 * papb * RLLL + pap2 * pbp1 * RLLL -
                                                pap1 * pbp2 * RLLL + m1 * ml4 * papb * RRLL)) +
                                        RR * (LLRR * m1 * m2 * ml4 * papb + LRRR * m2 * (-(p1p2 * papb) + 2.*pap1 * papb +
                                                pap2 * pbp1 - pap1 * pbp2) + m1 * papb * (m2s - 2.*pbp2) * RLRR -
                                                ml4 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RRRR)));

    return std::real(me0);
}

double FI::MVubo2u(int ieps)
{
    Npf * dd = new Npf(0., 0., m2s, m1s, 2.*papb, m2s - 2.*pap2,
                       m1l, m2l, m3l, m4l, mul, ieps);
    std::complex<double> la[7];

    dd->GetNpf(la, 2);
    std::complex<double> dd23 = la[2];
    std::complex<double> dd25 = la[4];
    std::complex<double> dd26 = la[5];
    std::complex<double> dd27 = la[6];
    dd->GetNpf(la, 1);
    std::complex<double> dd11 = la[0];
    std::complex<double> dd12 = la[1];
    std::complex<double> dd13 = la[2];
    dd->GetNpf(la, 0);
    std::complex<double> dd0 = la[0];
    delete dd;

    std::complex<double> me0 = 2.*ivu2s2 *
                               (dd23 * LR * LRLR * m1 * m2 * m2s * papb - dd0 * LL * LRLL * m1 * m2s * ml2 * papb +
                                dd0 * LL * LLLL * ml2 * ml4 * p1p2 * papb + 2.*dd0 * LL * LRLL * m1 * ml2 * pap2 * papb +
                                2.*dd12 * LL * LRLL * m1 * ml2 * pap2 * papb + dd0 * LL * LLLL * ml2 * ml4 * pap2 * pbp1 -
                                dd0 * LL * LLLL * ml2 * ml4 * pap1 * pbp2 - 2.*dd23 * LR * LRLR * m1 * m2 * pap2 * pbp2 -
                                dd23 * LRRL * m2s * p1p2 * papb * RL + 2.*dd26 * LRRL * p1p2 * pap2 * papb * RL +
                                dd23 * LRRL * m2s * pap2 * pbp1 * RL - 2.*dd26 * LRRL * std::pow(pap2, 2) * pbp1 * RL +
                                2.*dd25 * LRRL * m2s * papb * pbp1 * RL + dd23 * LRRL * m2s * pap1 * pbp2 * RL -
                                2.*dd26 * LRRL * pap1 * pap2 * pbp2 * RL - 4.*dd25 * LRRL * pap2 * pbp1 * pbp2 * RL +
                                dd0 * LL * m2 * ml2 * p1p2 * papb * RLLL + dd0 * LL * m2 * ml2 * pap2 * pbp1 * RLLL -
                                2.*dd0 * LL * m2 * ml2 * papb * pbp1 * RLLL - 2.*dd11 * LL * m2 * ml2 * papb * pbp1 * RLLL -
                                dd0 * LL * m2 * ml2 * pap1 * pbp2 * RLLL - dd23 * LR * m2s * p1p2 * papb * RLLR +
                                2.*dd26 * LR * p1p2 * pap2 * papb * RLLR + dd23 * LR * m2s * pap2 * pbp1 * RLLR -
                                2.*dd26 * LR * std::pow(pap2, 2) * pbp1 * RLLR + 2.*dd25 * LR * m2s * papb * pbp1 * RLLR +
                                dd23 * LR * m2s * pap1 * pbp2 * RLLR - 2.*dd26 * LR * pap1 * pap2 * pbp2 * RLLR -
                                4.*dd25 * LR * pap2 * pbp1 * pbp2 * RLLR + dd23 * m1 * m2 * m2s * papb * RL * RLRL -
                                2.*dd23 * m1 * m2 * pap2 * pbp2 * RL * RLRL +
                                2.*dd27 * (LR * LRLR * m1 * m2 * papb + 2.*LRRL * pap2 * pbp1 * RL + 2.*LR * pap2 * pbp1 * RLLR +
                                           m1 * m2 * papb * RL * RLRL) - dd0 * LLRR * m1 * m2 * ml2 * ml4 * papb * RR +
                                dd0 * LRRR * m2 * ml2 * p1p2 * papb * RR + dd0 * LRRR * m2 * ml2 * pap2 * pbp1 * RR -
                                2.*dd0 * LRRR * m2 * ml2 * papb * pbp1 * RR - 2.*dd11 * LRRR * m2 * ml2 * papb * pbp1 * RR -
                                dd0 * LRRR * m2 * ml2 * pap1 * pbp2 * RR - dd0 * m1 * m2s * ml2 * papb * RLRR * RR +
                                2.*dd0 * m1 * ml2 * pap2 * papb * RLRR * RR + 2.*dd12 * m1 * ml2 * pap2 * papb * RLRR * RR -
                                dd0 * LL * m1 * m2 * ml2 * ml4 * papb * RRLL +
                                dd13 * (LLRL * m1 * m2s * ml4 * papb * RL - LRRL * m2s * p1p2 * papb * RL +
                                        2.*LRRL * p1p2 * pap2 * papb * RL + LRRL * m2s * pap2 * pbp1 * RL -
                                        2.*LRRL * std::pow(pap2, 2) * pbp1 * RL + 2.*LRRL * m2s * papb * pbp1 * RL +
                                        LRRL * m2s * pap1 * pbp2 * RL -
                                        2.*LLRL * m1 * ml4 * pap2 * pbp2 * RL - 2.*LRRL * pap1 * pap2 * pbp2 * RL -
                                        4.*LRRL * pap2 * pbp1 * pbp2 * RL + LL * ml2 * (-(LRLL * m1 * m2s * papb) +
                                                m2 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RLLL) + m1 * m2 * m2s * papb * RL * RLRL -
                                        2.*m1 * m2 * pap2 * pbp2 * RL * RLRL + LRRR * m2 * ml2 * p1p2 * papb * RR +
                                        LRRR * m2 * ml2 * pap2 * pbp1 * RR - LRRR * m2 * ml2 * pap1 * pbp2 * RR -
                                        m1 * m2s * ml2 * papb * RLRR * RR +
                                        LR * (LLLR * m2 * ml4 * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                                                LRLR * m1 * m2 * (m2s * papb - 2.*pap2 * pbp2) -
                                                (m2s * (p1p2 * papb - pap2 * pbp1 - 2.*papb * pbp1 - pap1 * pbp2) +
                                                        2.*pap2 * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2 + 2.*pbp1 * pbp2)) * RLLR +
                                                m1 * ml4 * (m2s * papb - 2.*pap2 * pbp2) * RRLR) -
                                        m2 * ml4 * (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) * RL * RRRL) +
                                dd0 * ml2 * ml4 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RR * RRRR);

    return std::real(me0);
}
