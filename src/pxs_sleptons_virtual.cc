// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Computes the virtual part of the partonic cross section.

#include <cmath>
#include <complex>
#include <iostream>
#include <cstdlib>

#include "pxs.h"
#include "params.h"
#include "kinematics.h"
#include "utils.h"

#define IEPS 0
#define MUR2 params->murs

#define DIP // integrated dipole
#define VIR // virtual corrections

#define QSELF // external quark self-energy
#define TR1 // standard QCD vertex correction (quark-gluon-quark). IR divergent.
#define TR2 // SUSY-QCD vertex correction (squark-gluino-squark) IR safe.


// TR1 vertex correction.
double Vqcd_sleptons(const double S, const double T, Parameters *params) {
    double virt = 0.0;

    int aa = params->in1;
    int bb = params->in2;
    int ii = params->out1;
    int jj = params->out2;

    // For all s channel diagrams the same (leptons, sleptons, gauginos).
    int qs = iabs(aa / 3 - bb / 3);

    // Sets different Born kinematics
    FI *ff = new FI();
    ff->SetKinematic(params->mSL[ii - 10], params->mSL[jj - 10], S, T);


    // gluon-quark-quark coupling for strong corrections. (Basically g3s).
    struct Coupling Cs[2] = { params->gqq[bb][bb], params->gqq[aa][aa] };
    if (is_coupling_null(Cs, 2) == 1) {
        return virt;
    }

    ff->SetSCoupling(Cs);

    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
        for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
            struct Coupling Cw[4] = {
                params->vSLSL[i0][ii - 10][jj - 10],
                params->vqq[i0][bb][aa],
                params->vSLSL[i1][ii - 10][jj - 10],
                params->vqq[i1][bb][aa]
            };
            if (is_coupling_null(Cw, 4) == 1) {
                continue;
            }

            ff->SetPropagator(params->mv[i0],
                              params->mv[i1],
                              0.0,
                              0.0);
            ff->SetWCoupling(Cw);
            virt += ff->MVstr1sSL(0.0, 0.0, 2.0*ff->papb, 0.0, 0.0, 0.0, IEPS);
        }
    }
    delete ff;
    return virt;
}


// External quark self-energy at leg a.
double Vbu1_sleptons(const double S, const double T, Parameters *params) {
    double virt = 0.0;

    int aa = params->in1;
    int bb = params->in2;
    int ii = params->out1;
    int jj = params->out2;

    int qs = iabs(aa / 3 - bb / 3); // for all s channel diagrams the same (leptons, sleptons, gauginos)


    // Sets different Born kinematics
    FI *ff = new FI();
    ff->SetKinematic(params->mSL[ii - 10], params->mSL[jj - 10], S, T);

    for (int j0 = (aa / 3) * 6; j0 < (aa / 3 + 1) * 6; j0++) {
        struct Coupling Cs[2] = { params->GLqSQ[aa][j0], params->GLSQq[j0][aa] };
        if (is_coupling_null(Cs, 2) == 1) {
            continue;
        }
        ff->SetSCoupling(Cs);

        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                struct Coupling Cw[4] = {
                    params->vSLSL[i0][ii - 10][jj - 10], params->vqq[i0][bb][aa],
                    params->vSLSL[i1][ii - 10][jj - 10], params->vqq[i1][bb][aa]
                };
                if (is_coupling_null(Cw, 4) == 1) {
                    continue;
                }

                ff->SetPropagator(params->mv[i0], 
                                  params->mv[i1], 
                                  0.0, 
                                  0.0);
                ff->SetWCoupling(Cw);
                virt += ff->MVsbu1sSL(0.0, params->mGLs, params->mSQs[j0], IEPS);
            }
        }
    }

    delete ff;
    return virt;
}

// External quark self energy at leg b.
double Vbu2_sleptons(const double S, const double T, Parameters *params) {
    double virt = 0.0;

    int aa = params->in1;
    int bb = params->in2;
    int ii = params->out1;
    int jj = params->out2;

    int qs = iabs(aa / 3 - bb / 3); // for all s channel diagrams the same (leptons, sleptons, gauginos)


    // Sets different Born kinematics
    FI *ff = new FI();
    ff->SetKinematic(params->mSL[ii - 10], params->mSL[jj - 10], S, T);


    for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
        struct Coupling Cs[2] = { params->GLSQq[j0][bb], params->GLqSQ[bb][j0] };
        if (is_coupling_null(Cs, 2) == 1) {
            continue;
        }

        ff->SetSCoupling(Cs);

        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                struct Coupling Cw[4] = {
                    params->vSLSL[i0][ii - 10][jj - 10], params->vqq[i0][bb][aa],
                    params->vSLSL[i1][ii - 10][jj - 10], params->vqq[i1][bb][aa]
                };
                if (is_coupling_null(Cw, 4) == 0) {

                    ff->SetPropagator(params->mv[i0], 
                                      params->mv[i1], 
                                      0.0, 
                                      0.0);
                    ff->SetWCoupling(Cw);
                    virt += ff->MVsbu2sSL(0.0, params->mGLs, params->mSQs[j0], IEPS);
                }
            }
        }
    }

    delete ff;
    return virt;
}


// SUSY-QCD vertex correction.
double Vstr2_sleptons(const double S, const double T, Parameters *params) {
    double virt = 0.;

    int aa = params->in1;
    int bb = params->in2;
    int ii = params->out1;
    int jj = params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel

    FI *ff = new FI();
    ff->SetKinematic(params->mSL[ii - 10], params->mSL[jj - 10], S, T);

    for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
        for (int j1 = (aa / 3) * 6; j1 < (aa / 3 + 1) * 6; j1++) {
            struct Coupling Cs[2] = { params->GLqSQ[bb][j0], params->GLSQq[j1][aa] };
            if (is_coupling_null(Cs, 2) == 1) {
                continue;
            }
            ff->SetSCoupling(Cs);


            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                    struct Coupling Cw[4] = { params->vSLSL[i0][ii - 10][jj - 10], params->vqq[i0][bb][aa],
                                              params->vSLSL[i1][ii - 10][jj - 10], params->vqq[i1][bb][aa]
                    };
                    if (is_coupling_null(Cw, 4) == 0) {

                        ff->SetPropagator(params->mv[i0], 
                                          params->mv[i1], 
                                          0.0, 
                                          0.0);
                        ff->SetWCoupling(Cw);
                        virt += ff->MVstr2sSL(0.0, 0.0, 2.0*ff->papb, params->mGLs, 
                                              params->mSQs[j0], params->mSQs[j1], IEPS);
                    }
                }
            }
        }
    }
    delete ff;
    return virt;
}



double Virt_sleptons(const double S, const double T, Parameters *params) {

    const double g3s = std::norm(params->gqq[0][0].R);
    double result = 0.0;

#ifdef VIR

#ifdef TR1
    result += Vqcd_sleptons(S, T, params);
#endif

#ifdef QSELF
    result += 0.5 * Vbu1_sleptons(S, T, params) + 0.5 * Vbu2_sleptons(S, T, params);
#endif

#ifdef TR2
    result += Vstr2_sleptons(S, T, params);
#endif

#endif
    return result / g3s;
}

double Virt2_sleptons(const double S, const double T, Parameters *params) {
    const double g3s = std::norm(params->gqq[0][0].R);

    double result = 0.0;
    result += Vqcd_sleptons(S, T, params);

    return result / g3s;
}

// Integrated dipole.
double DipI_sleptons(const double S, const double T, Parameters *params) {
    const double lnmusq = std::log(params->murs / S);

#ifdef DIP    
    double murs = params->murs;
    double sab = S;
    
    if (IEPS == 2) {
        return 2.0 * born_sleptons(S, T, params);
    } else if (IEPS == 1) {
        return ((3.0 + 2.0 * lnmusq) * born_sleptons(S, T, params));
    } else if (IEPS == 0) {
      
      return born_sleptons(S, T, params)*(10.0 - pow(M_PI,2) + 3*log(murs/sab) + pow(log(murs/sab),2));
        // old: Jon has subtracted the (10 - pi^2) term from the coll.rem. 
       return ((3.0 + lnmusq) * lnmusq) * born_sleptons(S, T, params);
    }
    else {
        return 0.0;
    }
#endif
    return 0.0;
}
