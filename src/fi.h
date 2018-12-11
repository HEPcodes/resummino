// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Kinematics module.
//
// Notation: The incoming particle (antiparticle) has momentum pb (pa) and the
// outgoing particle (antiparticle) has momentum p2 (p1).

#ifndef FI_H_
#define FI_H_

#include "prm.h"

using namespace std;

class FI
{
public:
    // Constructors & Destructor
    FI() {};
    ~FI() {};

    // Sets kinematic variables for Born processes using masses of outgoing
    // particles mi(mj) and the common Mandelstam variables.
    void SetKinematic(const double mi, const double mj,
                      const double s, const double t);

    // Sets kinematic variables for real corrections using squared invariant mass of
    // the final state particles, squared transverse momentum pt2, scattering angle
    // theta, angle phi in the transverse plane and rapidity flag.
    void SetKinematic(const double mi, const double mj, const double s,
                      const double mi2, const double pt2, const double th,
                      const double ph, const int ys);

    // Sets kinematic variables for dipole term with parton pa as emitter.
    void SetDipKinematicA(double &x, double &fact);

    // Sets kinematic variables for dipole term with parton pb as emitter
    void SetDipKinematicB(double &x, double &fact);

    // Sets propagators using two mediator masses and the width.
    // Until now there is only the correct width in the s-channels and a small
    // width in the t channels with real quark emission.
    void SetPropagator(const double mass1, const double mass2, const double width1, const double width2);

    // Sets the maximal number of loop masses.
    // (Box diagrams -> 4 loop masses and 1 renormalization mass)
    void SetLoopMass(const double mls[5]);

    // Sets strong coupling coefficients for SUSY-QCD. Distinguishes between L- and R-type.
    void SetSCoupling(struct Coupling C[2]);

    // Sets general electroweak couplings. Four coupling coefficients for the
    // squared matrix element.
    // Distinguished between couplings to left- and right-handed fermions.
    // Complex conjugated coefficients in M*.
    void SetWCoupling(struct Coupling C[4]);

    // Born matrix elements
    double MBss();
    double MBtt();
    double MBuu();
    double MBst();
    double MBsu();
    double MBtu();
    double MB1ss();
    double MB1st();
    double MB1su();
    double MB2ss();
    double MBss_width();
    double MB1ss_width();
    double MB2ss_width();

    // Matrix elements for sleptons
    // Born
    double MBssSL();
    double MBssSL1();
    double MBssSL2();
    // Real Corrections
    double MGssSL();
    double MQssSL();
    double MQBssSL();
    // Self Energies
    double MVsbu1sSL(int eps);
    double MVsbu2sSL(int eps);
    // Vertex Corrections
    double MVstr1sSL(int eps);
    double MVstr2sSL(int eps);

    // Real matrix elements
    double MGss();
    double MGtt();
    double MGuu();
    double MGst();
    double MGsu();
    double MGtu();
    double MQss();
    double MQtt();
    double MQuu();
    double MQst();
    double MQsu();
    double MQtu();
    double MQBss();
    double MQBtt();
    double MQBuu();
    double MQBst();
    double MQBsu();
    double MQBtu();

    double MGss_width();
    double MQss_width();
    double MQBss_width();

    double MQttp();
    double MQuup();
    double MQBttp();
    double MQBuup();

    // Virtual matrix elements
    // Bubbles
    double MVsbu1s(int eps);
    double MVsbu1t(int eps);
    double MVsbu1u(int eps);
    double MVsbu2s(int eps);
    double MVsbu2t(int eps);
    double MVsbu2u(int eps);

    double MVtbu1s(int eps);
    double MVtbu1t(int eps);
    double MVtbu1u(int eps);
    double MVtbu2s(int eps);
    double MVtbu2t(int eps);
    double MVtbu2u(int eps);
    double MVtbu3s(int eps);
    double MVtbu3t(int eps);
    double MVtbu3u(int eps);
    double MVtbu4s(int eps);
    double MVtbu4t(int eps);
    double MVtbu4u(int eps);
    double MVtbu5s(int eps);
    double MVtbu5t(int eps);
    double MVtbu5u(int eps);

    double MVubu1s(int eps);
    double MVubu1t(int eps);
    double MVubu1u(int eps);
    double MVubu2s(int eps);
    double MVubu2t(int eps);
    double MVubu2u(int eps);
    double MVubu3s(int eps);
    double MVubu3t(int eps);
    double MVubu3u(int eps);
    double MVubu4s(int eps);
    double MVubu4t(int eps);
    double MVubu4u(int eps);
    double MVubu5s(int eps);
    double MVubu5t(int eps);
    double MVubu5u(int eps);

    double MVsbu1s_width(int eps);
    double MVsbu2s_width(int eps);

    // Triangles
    double MVstr1s(int eps);
    double MVstr1t(int eps);
    double MVstr1u(int eps);
    double MVstr2s(int eps);
    double MVstr2t(int eps);
    double MVstr2u(int eps);
    double MVstr3t(int eps);
    double MVstr3u(int eps);

    double MVttr1s(int eps);
    double MVttr1t(int eps);
    double MVttr1u(int eps);
    double MVttr2s(int eps);
    double MVttr2t(int eps);
    double MVttr2u(int eps);
    double MVttr3s(int eps);
    double MVttr3t(int eps);
    double MVttr3u(int eps);
    double MVttr4s(int eps);
    double MVttr4t(int eps);
    double MVttr4u(int eps);

    double MVutr1s(int eps);
    double MVutr1t(int eps);
    double MVutr1u(int eps);
    double MVutr2s(int eps);
    double MVutr2t(int eps);
    double MVutr2u(int eps);
    double MVutr3s(int eps);
    double MVutr3t(int eps);
    double MVutr3u(int eps);
    double MVutr4s(int eps);
    double MVutr4t(int eps);
    double MVutr4u(int eps);

    double MVstr1s_width(int eps);
    double MVstr2s_width(int eps);

    // Boxes
    double MVtbo1s(int eps);
    double MVtbo1t(int eps);
    double MVtbo1u(int eps);
    double MVtbo2s(int eps);
    double MVtbo2t(int eps);
    double MVtbo2u(int eps);

    double MVubo1s(int eps);
    double MVubo1t(int eps);
    double MVubo1u(int eps);
    double MVubo2s(int eps);
    double MVubo2t(int eps);
    double MVubo2u(int eps);

    // Counterterms
    double MVsct1s(int eps);
    double MVsct1t(int eps);
    double MVsct1u(int eps);
    double MVsct2s(int eps);
    double MVsct2t(int eps);
    double MVsct2u(int eps);

    double MVtct1s(int eps);
    double MVtct1t(int eps);
    double MVtct1u(int eps);
    double MVtct2s(int eps);
    double MVtct2t(int eps);
    double MVtct2u(int eps);
    double MVtct3s(int eps);
    double MVtct3t(int eps);
    double MVtct3u(int eps);
    double MVtct4s(int eps);
    double MVtct4t(int eps);
    double MVtct4u(int eps);
    double MVtct5s(int eps);
    double MVtct5t(int eps);
    double MVtct5u(int eps);

    double MVuct1s(int eps);
    double MVuct1t(int eps);
    double MVuct1u(int eps);
    double MVuct2s(int eps);
    double MVuct2t(int eps);
    double MVuct2u(int eps);
    double MVuct3s(int eps);
    double MVuct3t(int eps);
    double MVuct3u(int eps);
    double MVuct4s(int eps);
    double MVuct4t(int eps);
    double MVuct4u(int eps);
    double MVuct5s(int eps);
    double MVuct5t(int eps);
    double MVuct5u(int eps);


public:
    // Widths
    double prop_width[2];

private:
    // Kinematics
    double m1, m2, m1s, m2s;
    double papb, pap1, pap2, pap3, pbp1, pbp2, pbp3, p1p2, p1p3, p2p3;

    // Propagators
    double ivs;
    complex<double> ivs1s1, ivs1s2, ivs2s1, ivs2s2, ivs3v1, ivs3v2;
    double ivu3, ivt3, ivt1s1, ivu2s1, ivt1s2, ivt2s1, ivt2s2, ivu1s1, ivu1s2, ivu2s2;

    // Couplings
    // Strong Couplings
    complex<double> LL, LR, RR, RL;

    // Weak Couplings
    complex<double> LLLL, LLLR, LLRL, LRLL, RLLL, LLRR, LRLR, LRRL;
    complex<double> RRRR, RRRL, RRLR, RLRR, LRRR, RRLL, RLRL, RLLR;

    // Loops
    double mul, m1l, m2l, m3l, m4l;
    double ml1, ml2, ml3, ml4;
};

#endif
