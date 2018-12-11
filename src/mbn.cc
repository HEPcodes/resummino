// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Born matrix elements.
//
// Notes on notation:
// - Squared Born matrix elements in D = 4 - 2 epsilon dimensions.
// - The incoming particle (antiparticle) has momentum pb (pa) and the
// - outgoing particle (antiparticle) has momentum p2 (p1).
// - Four general electroweak coupling coefficients referring to the four
//   vertices. Distinguishes between left- and right-handed.
// - General vector coupling: \Gamma_i = \gamma_mu (L_i P_L + R_i P_R )
//   e.g. LLLL = L_1 L_2 L_3^* L_4^* (vqq,vll,vqq,vll).
// - General scalar coupling analoge.
// - MBss:
//   - M = squared Matrix element
//   - B = Born level
//   - s = channel 1 is s-channel
//   - s = channel 2 is s-channel
// - Spin and color sum and average are implemented in hxs.cc.


#include <complex>
#include "fi.h"

// Process: q + \bar{q} -> f + \bar{f}.

// Squared Born matrix elements of O(epsilon^0).
double FI::MBss()
{
    return std::real(8.*ivs3v1 * ivs3v2 * (m1 * m2 * papb * (LLRL + LRRR + RLLL + RRLR) +
                                           2.*(pap1 * pbp2 * (LRLR + RLRL) + pap2 * pbp1 * (LLLL + RRRR))));
}


double FI::MBtt()
{
    return std::real(4.*ivt1s1 * ivt1s2 * pap1 * pbp2 * (LLLL + LRLR + RLRL + RRRR));
}
double FI::MBuu()
{
    return std::real(4.*ivu2s1 * ivu2s2 * pap2 * pbp1 * (LLLL + LRLR + RLRL + RRRR));
}
double FI::MBst()
{
    return std::real(4.*ivs3v1 * ivt1s2 * (2.*pap1 * pbp2 * (LRLR + RLRL)
                                           + m1 * m2 * papb * (LLRL + RRLR)));
}
double FI::MBsu()
{
    return std::real(-4.*ivs3v1 * ivu2s2 * (m1 * m2 * papb * (LRLR + RLRL)
                                            + 2.*pap2 * pbp1 * (LLRL + RRLR)));
}
double FI::MBtu()
{
    return std::real(-2.*ivt1s1 * ivu2s2 * ((LRLR + RLRL) * m1 * m2 * papb +
                                            (LLLL + RRRR) * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)));
}

// Squared Born matrix elements of O(epsilon^1).
double FI::MB1ss()
{
    return std::real(-8.*ivs3v1 * ivs3v2 *
                     (LLRL * m1 * m2 * papb + LRRR * m1 * m2 * papb + LLLL * p1p2 * papb +
                      LRLR * p1p2 * papb + 3.*LLLL * pap2 * pbp1 - 3.*LRLR * pap2 * pbp1 - 3.*LLLL * pap1 * pbp2 +
                      3.*LRLR * pap1 * pbp2 + m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - 3.*pap2 * pbp1 * RLRL +
                      3.*pap1 * pbp2 * RLRL + m1 * m2 * papb * RRLR + p1p2 * papb * RRRR + 3.*pap2 * pbp1 * RRRR -
                      3.*pap1 * pbp2 * RRRR));
}


double FI::MB1st()
{
    return std::real(-4.*ivs3v1 * ivt1s2 *
                     (LLRL * m1 * m2 * papb + (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) *
                      (LRLR + RLRL) + m1 * m2 * papb * RRLR));
}

double FI::MB1su()
{
    return std::real(4.*ivs3v1 * ivu2s2 *
                     (LRLR * m1 * m2 * papb + LLRL * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) +
                      m1 * m2 * papb * RLRL + p1p2 * papb * RRLR + pap2 * pbp1 * RRLR - pap1 * pbp2 * RRLR));
}

//Squared Born matrix elements of O(epsilon^2).
double FI::MB2ss()
{
    return std::real(16.*ivs3v1 * ivs3v2 * (pap2 * pbp1 - pap1 * pbp2) *
                     (LLLL - LRLR - RLRL + RRRR));
}

// Squared Born matrix element for slepton pair production.
// Process: q + \bar{q} -> sl + \bar{sl}.

// Squared Born matrix elements of O(epsilon^0).
double FI::MBssSL()
{
    return std::real(8.*ivs3v1 * ivs3v2 * (LRLR + LLLL) * (2.*pap1 * pap2 - pap1 * m2s - pap2 * m1s));

}
