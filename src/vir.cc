// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
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
#include "prm.h"
#include "fi.h"

#define IEPS 0
#define DRBAR
#define MUR2 Params->murs

#define SS
#define TT
#define UU
#define ST
#define SU
#define TU

#define SQUARKS // Include squarks in corrections.

static inline bool FastFlag(struct Coupling *C, int N)
{
    for (int i0 = 0; i0 < N; i0++) {
        if (C[i0].L == 0. && C[i0].R == 0.) {
            return 1;
        }
    }
    return 0;
}

static inline int iabs(int x)
{
    return (x < 0 ? -x : x);
}

double Vqcd(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
    int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

    double g3s = std::norm(Params->gqq[0][0].R);

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }

    struct Coupling Cs[2] = { Params->gqq[bb][bb], Params->gqq[aa][aa] };
    if (FastFlag(Cs, 2) == 1) {
        return virt;
    }
    ff->SetSCoupling(Cs);

    if (Params->out1 >= 10) {
#ifdef SS
        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                struct Coupling Cw[4] = { Params->vSLSL[i0][ii - 10][jj - 10], Params->vqq[i0][bb][aa],
                           Params->vSLSL[i1][ii - 10][jj - 10], Params->vqq[i1][bb][aa]
                };
                if (FastFlag(Cw, 4) == 1) {
                    continue;
                }

                ff->SetPropagator(Params->mv[i0], Params->mv[i1] , Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                double ml0[5] = { Params->murs, Params->mq[bb], 0., Params->mq[aa], 0. };
                ff->SetLoopMass(ml0);
                virt += ff->MVstr1sSL(IEPS);
            }
        }
#endif
    } else {
#ifdef SS
        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                           Params->vCHCH[i1][ii][jj], Params->vqq[i1][bb][aa]
                };
                if (FastFlag(Cw, 4) == 1) {
                    continue;
                }

                ff->SetPropagator(Params->mv[i0], Params->mv[i1] , Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                double ml0[5] = { Params->murs, Params->mq[bb], 0., Params->mq[aa], 0. };
                ff->SetLoopMass(ml0);
                virt += ff->MVstr1s(IEPS);
            }
        }
#endif
#ifdef TT
        if (qt >= 0) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    double ml1[5] = { Params->murs, Params->mq[bb], 0., Params->mq[aa], Params->mSQ[i0] };
                    ff->SetLoopMass(ml1);
                    virt += ff->MVtbo1t(IEPS);
                    double ml2[5] = { Params->murs, Params->mq[aa], 0., Params->mSQ[i0], 0. };
                    ff->SetLoopMass(ml2);
                    virt += ff->MVttr1t(IEPS);
                    double ml3[5] = { Params->murs, Params->mq[bb], 0., Params->mSQ[i0], 0. };
                    ff->SetLoopMass(ml3);
                    virt += ff->MVttr2t(IEPS);
                    double ml4[5] = { Params->murs, Params->mSQ[i0], 0., 0., 0. };
                    ff->SetLoopMass(ml4);
                    virt += ff->MVtbu3t(IEPS);
                    virt += ff->MVtct3t(IEPS);

#ifdef DRBAR
                    if (IEPS == 0) {
                        virt -= g3s * ff->MBtt(); // DRbar
                    }
#endif
                }
            }
        }
#endif
#ifdef UU
        if (qu <= 1) {
            for (int i0 = 6 * qu; i0 < 6 * qu + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[ii][bb][i0], Params->CHSQq[jj][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    double ml1[5] = { Params->murs, Params->mq[bb], 0., Params->mq[aa], Params->mSQ[i0] };
                    ff->SetLoopMass(ml1);
                    virt += ff->MVubo1u(IEPS);
                    double ml2[5] = { Params->murs, Params->mq[aa], 0., Params->mSQ[i0], 0. };
                    ff->SetLoopMass(ml2);
                    virt += ff->MVutr1u(IEPS);
                    double ml3[5] = { Params->murs, Params->mq[bb], 0., Params->mSQ[i0], 0. };
                    ff->SetLoopMass(ml3);
                    virt += ff->MVutr2u(IEPS);
                    double ml4[5] = { Params->murs, Params->mSQ[i0], 0., 0., 0. };
                    ff->SetLoopMass(ml4);
                    virt += ff->MVubu3u(IEPS);
                    virt += ff->MVuct3u(IEPS);

#ifdef DRBAR
                    if (IEPS == 0) {
                        virt -= g3s * ff->MBuu(); // DRbar
                    }
#endif
                }
            }
        }
#endif
#ifdef ST
        if (qt >= 0) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw1[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    struct Coupling Cw2[4] = { Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa],
                               Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa]
                    };
                    if (FastFlag(Cw1, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw1);
                    double ml0[5] = { Params->murs, Params->mq[bb], 0., Params->mq[aa], 0. };
                    ff->SetLoopMass(ml0);
                    virt += ff->MVstr1t(IEPS);

#ifdef DRBAR
                    if (IEPS == 0) {
                        virt -= g3s * ff->MBst(); // DRbar
                    }
#endif

                    ff->SetPropagator(Params->mSQ[i1], Params->mv[i0] , Params->mSQ[i1] * 1.e-2, Params->Gv[i0]);
                    ff->SetWCoupling(Cw2);
                    double ml1[5] = { Params->murs, Params->mq[bb], 0., Params->mq[aa], Params->mSQ[i1] };
                    ff->SetLoopMass(ml1);
                    virt += ff->MVtbo1s(IEPS);
                    double ml2[5] = { Params->murs, Params->mq[aa], 0., Params->mSQ[i1], 0. };
                    ff->SetLoopMass(ml2);
                    virt += ff->MVttr1s(IEPS);
                    double ml3[5] = { Params->murs, Params->mq[bb], 0., Params->mSQ[i1], 0. };
                    ff->SetLoopMass(ml3);
                    virt += ff->MVttr2s(IEPS);
                    double ml4[5] = { Params->murs, Params->mSQ[i1], 0., 0., 0. };
                    ff->SetLoopMass(ml4);
                    virt += ff->MVtbu3s(IEPS);
                    virt += ff->MVtct3s(IEPS);
                }
            }
        }
#endif
#ifdef SU
        if (qu <= 1) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw1[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    struct Coupling Cw2[4] = { Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa],
                               Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa]
                    };
                    if (FastFlag(Cw1, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw1);
                    double ml0[5] = { Params->murs, Params->mq[bb], 0., Params->mq[aa], 0. };
                    ff->SetLoopMass(ml0);
                    virt += ff->MVstr1u(IEPS);

#ifdef DRBAR
                    if (IEPS == 0) {
                        virt -= g3s * ff->MBsu(); // DRbar
                    }
#endif

                    ff->SetPropagator(Params->mSQ[i1], Params->mv[i0], Params->mSQ[i1] * 1.e-2, Params->Gv[i0]);
                    ff->SetWCoupling(Cw2);
                    double ml1[5] = { Params->murs, Params->mq[bb], 0., Params->mq[aa], Params->mSQ[i1] };
                    ff->SetLoopMass(ml1);
                    virt += ff->MVubo1s(IEPS);
                    double ml2[5] = { Params->murs, Params->mq[aa], 0., Params->mSQ[i1], 0. };
                    ff->SetLoopMass(ml2);
                    virt += ff->MVutr1s(IEPS);
                    double ml3[5] = { Params->murs, Params->mq[bb], 0., Params->mSQ[i1], 0. };
                    ff->SetLoopMass(ml3);
                    virt += ff->MVutr2s(IEPS);
                    double ml4[5] = { Params->murs, Params->mSQ[i1], 0., 0., 0. };
                    ff->SetLoopMass(ml4);
                    virt += ff->MVubu3s(IEPS);
                    virt += ff->MVuct3s(IEPS);
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw1[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    struct Coupling Cw2[4] = { Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa],
                               Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa]
                    };
                    if (FastFlag(Cw1, 4) == 1) {
                        continue;
                    }
                    //double mp1[2] = { Params->mSQ[i0], Params->mSQ[i1] };
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw1);
                    double mt1[5] = { Params->murs, Params->mq[bb], 0., Params->mq[aa], Params->mSQ[i0] };
                    ff->SetLoopMass(mt1);
                    virt += ff->MVtbo1u(IEPS);
                    double mt2[5] = { Params->murs, Params->mq[aa], 0., Params->mSQ[i0], 0. };
                    ff->SetLoopMass(mt2);
                    virt += ff->MVttr1u(IEPS);
                    double mt3[5] = { Params->murs, Params->mq[bb], 0., Params->mSQ[i0], 0. };
                    ff->SetLoopMass(mt3);
                    virt += ff->MVttr2u(IEPS);
                    double mt4[5] = { Params->murs, Params->mSQ[i0], 0., 0., 0. };
                    ff->SetLoopMass(mt4);
                    virt += ff->MVtbu3u(IEPS);
                    virt += ff->MVtct3u(IEPS);

#ifdef DRBAR
                    if (IEPS == 0) {
                        virt -= 2.*g3s * ff->MBtu(); // DRbar
                    }
#endif

                    ff->SetPropagator(Params->mSQ[i1], Params->mSQ[i0], Params->mSQ[i1] * 1.e-2, Params->mSQ[i0] * 1.e-2);
                    ff->SetWCoupling(Cw2);
                    double mu1[5] = { Params->murs, Params->mq[bb], 0., Params->mq[aa], Params->mSQ[i1] };
                    ff->SetLoopMass(mu1);
                    virt += ff->MVubo1t(IEPS);
                    double mu2[5] = { Params->murs, Params->mq[aa], 0., Params->mSQ[i1], 0. };
                    ff->SetLoopMass(mu2);
                    virt += ff->MVutr1t(IEPS);
                    double mu3[5] = { Params->murs, Params->mq[bb], 0., Params->mSQ[i1], 0. };
                    ff->SetLoopMass(mu3);
                    virt += ff->MVutr2t(IEPS);
                    double mu4[5] = { Params->murs, Params->mSQ[i1], 0., 0., 0. };
                    ff->SetLoopMass(mu4);
                    virt += ff->MVubu3t(IEPS);
                    virt += ff->MVuct3t(IEPS);
                }
            }
        }
#endif
    }


    delete ff;

    return virt;
}

// Gluino loops
double Vbu1(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    for (int j0 = (aa / 3) * 6; j0 < (aa / 3 + 1) * 6; j0++) {
        for (int j1 = (aa / 3) * 3; j1 < (aa / 3 + 1) * 3; j1++) {
            double mls[5] = { Params->murs, Params->mGL, Params->mSQ[j0], 0., 0. };
            struct Coupling Cs[2] = { Params->GLqSQ[j1][j0], Params->GLSQq[j0][aa] };
            if (FastFlag(Cs, 2) == 1) {
                continue;
            }
            ff->SetLoopMass(mls);
            ff->SetSCoupling(Cs);

            if (Params->out1 >= 10) {
#ifdef SS
                for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                        struct Coupling Cw[4] = { Params->vSLSL[i0][ii - 10][jj - 10], Params->vqq[i0][bb][aa],
                                   Params->vSLSL[i1][ii - 10][jj - 10], Params->vqq[i1][bb][aa]
                        };
                        if (FastFlag(Cw, 4) == 1) {
                            continue;
                        }

                        ff->SetPropagator(Params->mv[i0], Params->mv[i1] , Params->Gv[i0], Params->Gv[i1]);
                        ff->SetWCoupling(Cw);
                        virt += ff->MVsbu1sSL(IEPS);
                    }
                }
#endif
            } else {

#ifdef SS
                for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                        struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj],
                                   Params->vqq[i0][bb][j1],
                                   Params->vCHCH[i1][ii][jj],
                                   Params->vqq[i1][bb][aa]
                        };
                        if (FastFlag(Cw, 4) == 1) {
                            continue;
                        }
                        ff->SetPropagator(Params->mv[i0], Params->mv[i1] , Params->Gv[i0], Params->Gv[i1]);
                        ff->SetWCoupling(Cw);
                        virt += ff->MVsbu1s(IEPS);
                    }
                }
#endif
#ifdef TT
                if (qt >= 0) {
                    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                        for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0],
                                       Params->CHSQq[ii][i0][j1],
                                       Params->CHqSQ[jj][bb][i1],
                                       Params->CHSQq[ii][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 1) {
                                continue;
                            }
                            ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                            ff->SetWCoupling(Cw);
                            virt += ff->MVtbu1t(IEPS);
                            if (j1 == aa) {
                                virt += ff->MVtct1t(IEPS);
                            }
                        }
                    }
                }
#endif
#ifdef UU
                if (qu <= 1) {
                    for (int i0 = 6 * qu; i0 < 6 * qu + 6; i0++) {
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[ii][bb][i0],
                                       Params->CHSQq[jj][i0][j1],
                                       Params->CHqSQ[ii][bb][i1],
                                       Params->CHSQq[jj][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 1) {
                                continue;
                            }
                            ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                            ff->SetWCoupling(Cw);
                            virt += ff->MVubu1u(IEPS);
                            if (j1 == aa) {
                                virt += ff->MVuct1u(IEPS);
                            }
                        }
                    }
                }
#endif
#ifdef ST
                if (qt >= 0) {
                    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                        for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                            struct Coupling Cw1[4] = { Params->vCHCH[i0][ii][jj],
                                       Params->vqq[i0][bb][j1],
                                       Params->CHqSQ[jj][bb][i1],
                                       Params->CHSQq[ii][i1][aa]
                            };
                            struct Coupling Cw2[4] = { Params->CHqSQ[jj][bb][i1],
                                       Params->CHSQq[ii][i1][j1],
                                       Params->vCHCH[i0][ii][jj],
                                       Params->vqq[i0][bb][aa]
                            };
                            if (FastFlag(Cw1, 4) == 1) {
                                continue;
                            }
                            ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                            ff->SetWCoupling(Cw1);
                            virt += ff->MVsbu1t(IEPS);
                            if (j1 == aa) {
                                virt += ff->MVsct1t(IEPS);
                            }

                            ff->SetPropagator(Params->mSQ[i1], Params->mv[i0], Params->mSQ[i1] * 1.e-2, Params->Gv[i0]);
                            ff->SetWCoupling(Cw2);
                            virt += ff->MVtbu1s(IEPS);
                            if (j1 == aa) {
                                virt += ff->MVtct1s(IEPS);
                            }
                        }
                    }
                }
#endif
#ifdef SU
                if (qu <= 1) {
                    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw1[4] = { Params->vCHCH[i0][ii][jj],
                                       Params->vqq[i0][bb][j1],
                                       Params->CHqSQ[ii][bb][i1],
                                       Params->CHSQq[jj][i1][aa]
                            };
                            struct Coupling Cw2[4] = { Params->CHqSQ[ii][bb][i1],
                                       Params->CHSQq[jj][i1][j1],
                                       Params->vCHCH[i0][ii][jj],
                                       Params->vqq[i0][bb][aa]
                            };
                            if (FastFlag(Cw1, 4) == 1) {
                                continue;
                            }
                            ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], 0);
                            ff->SetWCoupling(Cw1);
                            virt += ff->MVsbu1u(IEPS);
                            if (j1 == aa) {
                                virt += ff->MVsct1u(IEPS);
                            }

                            ff->SetPropagator(Params->mSQ[i1], Params->mv[i0], Params->mSQ[i1] * 1.e-2, Params->Gv[i0]);
                            ff->SetWCoupling(Cw2);
                            virt += ff->MVubu1s(IEPS);
                            if (j1 == aa) {
                                virt += ff->MVuct1s(IEPS);
                            }
                        }
                    }
                }
#endif
#ifdef TU
                if (qt >= 0 && qu <= 1) {
                    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw1[4] = { Params->CHqSQ[jj][bb][i0],
                                       Params->CHSQq[ii][i0][j1],
                                       Params->CHqSQ[ii][bb][i1],
                                       Params->CHSQq[jj][i1][aa]
                            };
                            struct Coupling Cw2[4] = { Params->CHqSQ[ii][bb][i1],
                                       Params->CHSQq[jj][i1][j1],
                                       Params->CHqSQ[jj][bb][i0],
                                       Params->CHSQq[ii][i0][aa]
                            };
                            if (FastFlag(Cw1, 4) == 1) {
                                continue;
                            }
                            ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], 0, 0);
                            ff->SetWCoupling(Cw1);
                            virt += ff->MVtbu1u(IEPS);
                            if (j1 == aa) {
                                virt += ff->MVtct1u(IEPS);
                            }

                            ff->SetPropagator(Params->mSQ[i1], Params->mSQ[i0], Params->mSQ[i1] * 1.e-2, Params->mSQ[i0] * 1.e-2);
                            ff->SetWCoupling(Cw2);
                            virt += ff->MVubu1t(IEPS);
                            if (j1 == aa) {
                                virt += ff->MVuct1t(IEPS);
                            }
                        }
                    }
                }
#endif
            }

        }
    }

    delete ff;

    return virt;
}

double Vbu2(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
        for (int j1 = (bb / 3) * 3; j1 < (bb / 3 + 1) * 3; j1++) {
            double mls[5] = { Params->murs, Params->mGL, Params->mSQ[j0], 0., 0. };
            struct Coupling Cs[2] = { Params->GLqSQ[bb][j0], Params->GLSQq[j0][j1] };
            if (FastFlag(Cs, 2) == 1) {
                continue;
            }
            ff->SetLoopMass(mls);
            ff->SetSCoupling(Cs);

            if (Params->out1 >= 10) {
#ifdef SS
                for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                        struct Coupling Cw[4] = { Params->vSLSL[i0][ii - 10][jj - 10], Params->vqq[i0][bb][aa],
                                   Params->vSLSL[i1][ii - 10][jj - 10], Params->vqq[i1][bb][aa]
                        };
                        if (FastFlag(Cw, 4) == 0) {

                            ff->SetPropagator(Params->mv[i0], Params->mv[i1] , Params->Gv[i0], Params->Gv[i1]);
                            ff->SetWCoupling(Cw);
#ifdef SQUARKS
                            virt += ff->MVsbu2sSL(IEPS);
#endif
#ifdef DEC_SQUARKS
                            virt += 0;
#endif
                        }
                    }
                }
#endif
            } else {

#ifdef SS
                for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                        struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][j1][aa],
                                   Params->vCHCH[i1][ii][jj], Params->vqq[i1][bb][aa]
                        };
                        if (FastFlag(Cw, 4) == 0) {
                            ff->SetPropagator(Params->mv[i0], Params->mv[i1] , Params->Gv[i0], Params->Gv[i1]);
                            ff->SetWCoupling(Cw);
                            virt += ff->MVsbu2s(IEPS);
                        }
                    }
                }
#endif
#ifdef TT
                if (qt >= 0) {
                    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                        for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[jj][j1][i0], Params->CHSQq[ii][i0][aa],
                                       Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVtbu2t(IEPS);
                                if (j1 == bb) {
                                    virt += ff->MVtct2t(IEPS);
                                }
                            }
                        }
                    }
                }
#endif
#ifdef UU
                if (qu <= 1) {
                    for (int i0 = 6 * qu; i0 < 6 * qu + 6; i0++) {
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[ii][j1][i0], Params->CHSQq[jj][i0][aa],
                                       Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVubu2u(IEPS);
                                if (j1 == bb) {
                                    virt += ff->MVuct2u(IEPS);
                                }
                            }
                        }
                    }
                }
#endif
#ifdef ST
                if (qt >= 0) {
                    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                        for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                            struct Coupling Cw1[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][j1][aa],
                                       Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                            };
                            struct Coupling Cw2[4] = { Params->CHqSQ[jj][j1][i1], Params->CHSQq[ii][i1][aa],
                                       Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa]
                            };
                            if (FastFlag(Cw1, 4) == 0) {
                                ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw1);
                                virt += ff->MVsbu2t(IEPS);
                                if (j1 == bb) {
                                    virt += ff->MVsct2t(IEPS);
                                }

                                ff->SetPropagator(Params->mSQ[i1], Params->mv[i0], Params->mSQ[i1] * 1.e-2, Params->Gv[i0]);
                                ff->SetWCoupling(Cw2);
                                virt += ff->MVtbu2s(IEPS);
                                if (j1 == bb) {
                                    virt += ff->MVtct2s(IEPS);
                                }
                            }
                        }
                    }
                }
#endif
#ifdef SU
                if (qu <= 1) {
                    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw1[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][j1][aa],
                                       Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                            };
                            struct Coupling Cw2[4] = { Params->CHqSQ[ii][j1][i1], Params->CHSQq[jj][i1][aa],
                                       Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa]
                            };
                            if (FastFlag(Cw1, 4) == 0) {
                                ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw1);
                                virt += ff->MVsbu2u(IEPS);
                                if (j1 == bb) {
                                    virt += ff->MVsct2u(IEPS);
                                }

                                ff->SetPropagator(Params->mSQ[i1], Params->mv[i0], Params->mSQ[i1] * 1.e-2, Params->Gv[i0]);
                                ff->SetWCoupling(Cw2);
                                virt += ff->MVubu2s(IEPS);
                                if (j1 == bb) {
                                    virt += ff->MVuct2s(IEPS);
                                }
                            }
                        }
                    }
                }
#endif
#ifdef TU
                if (qt >= 0 && qu <= 1) {
                    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw1[4] = { Params->CHqSQ[jj][j1][i0], Params->CHSQq[ii][i0][aa],
                                       Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                            };
                            struct Coupling Cw2[4] = { Params->CHqSQ[ii][j1][i1], Params->CHSQq[jj][i1][aa],
                                       Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa]
                            };
                            if (FastFlag(Cw1, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw1);
                                virt += ff->MVtbu2u(IEPS);
                                if (j1 == bb) {
                                    virt += ff->MVtct2u(IEPS);
                                }

                                ff->SetPropagator(Params->mSQ[i1], Params->mSQ[i0] , Params->mSQ[i1] * 1.e-2, Params->mSQ[i0] * 1.e-2);
                                ff->SetWCoupling(Cw2);
                                virt += ff->MVubu2t(IEPS);
                                if (j1 == bb) {
                                    virt += ff->MVuct2t(IEPS);
                                }
                            }
                        }
                    }
                }
#endif
            }
        }
    }

    delete ff;

    return virt;
}

double Vbu4(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    if (qt >= 0) {
        for (int j0 = 6 * qt; j0 < 6 * qt + 6; j0++) {
            for (int j1 = 6 * qt; j1 < 6 * qt + 6; j1++) {
                for (int j2 = 3 * qt; j2 < 3 * qt + 3; j2++) {
                    double mls[5] = { Params->murs, Params->mq[j2], Params->mGL, Params->mSQ[j0], Params->mSQ[j1] };
                    struct Coupling Cs[2] = { Params->GLSQq[j0][j2], Params->GLqSQ[j2][j1] };
                    if (FastFlag(Cs, 2) == 1) {
                        continue;
                    }
                    ff->SetLoopMass(mls);
                    ff->SetSCoupling(Cs);

                    if (Params->out1 >= 10) {
                        // No s-channel
                    } else {

#ifdef TT
                        for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][j0], Params->CHSQq[ii][j1][aa],
                                       Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 1) {
                                continue;
                            }
                            ff->SetPropagator(Params->mSQ[i1], Params->mSQ[i1], Params->mSQ[i1] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                            ff->SetWCoupling(Cw);
                            virt += ff->MVtbu4t(IEPS);
                            if (j1 == j0) {
                                virt += ff->MVtct4t(IEPS);
                            }
                        }
#endif
#ifdef ST
                        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][j0], Params->CHSQq[ii][j1][aa],
                                       Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mv[i0], Params->mv[i0], Params->Gv[i0], Params->Gv[i0]);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVtbu4s(IEPS);
                                if (j1 == j0) {
                                    virt += ff->MVtct4s(IEPS);
                                }
                            }
                        }
#endif
#ifdef TU
                        if (qu <= 1) {
                            for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                                struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][j0], Params->CHSQq[ii][j1][aa],
                                           Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                                };
                                if (FastFlag(Cw, 4) == 0) {
                                    ff->SetPropagator(Params->mSQ[i1], Params->mSQ[i1], Params->mSQ[i1] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                    ff->SetWCoupling(Cw);
                                    virt += ff->MVtbu4u(IEPS);
                                    if (j1 == j0) {
                                        virt += ff->MVtct4u(IEPS);
                                    }
                                }
                            }
                        }
#endif
                    }
                }
            }
        }
        if (qu <= 1) {
            for (int j0 = 6 * qu; j0 < 6 * qu + 6; j0++) {
                for (int j1 = 6 * qu; j1 < 6 * qu + 6; j1++) {
                    for (int j2 = 3 * qu; j2 < 3 * qu + 3; j2++) {
                        double mls[5] = { Params->murs, Params->mq[j2], Params->mGL, Params->mSQ[j0], Params->mSQ[j1] };
                        struct Coupling Cs[2] = { Params->GLSQq[j0][j2], Params->GLqSQ[j2][j1] };
                        if (FastFlag(Cs, 2) == 1) {
                            continue;
                        }
                        ff->SetLoopMass(mls);
                        ff->SetSCoupling(Cs);
#ifdef UU
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[ii][bb][j0], Params->CHSQq[jj][j1][aa],
                                       Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[i1], Params->mSQ[i1], Params->mSQ[i1] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVubu4u(IEPS);
                                if (j1 == j0) {
                                    virt += ff->MVuct4u(IEPS);
                                }
                            }
                        }
#endif
#ifdef SU
                        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[ii][bb][j0], Params->CHSQq[jj][j1][aa],
                                       Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mv[i0], Params->mv[i0], Params->Gv[i0], Params->Gv[i0]);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVubu4s(IEPS);
                                if (j1 == j0) {
                                    virt += ff->MVuct4s(IEPS);
                                }
                            }
                        }
#endif
#ifdef TU
                        if (qt >= 0) {
                            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                                struct Coupling Cw[4] = { Params->CHqSQ[ii][bb][j0], Params->CHSQq[jj][j1][aa],
                                           Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa]
                                };
                                if (FastFlag(Cw, 4) == 0) {
                                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i0], Params->mSQ[i0] * 1.e-2, Params->mSQ[i0] * 1.e-2);
                                    ff->SetWCoupling(Cw);
                                    virt += ff->MVubu4t(IEPS);
                                    if (j1 == j0) {
                                        virt += ff->MVuct4t(IEPS);
                                    }
                                }
                            }
                        }
#endif
                    }
                }
            }
        }
    }

    delete ff;

    return virt;
}

double Vbu5(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    if (qt >= 0) {
        for (int j0 = 6 * qt; j0 < 6 * qt + 6; j0++) {
            for (int j1 = 6 * qt; j1 < 6 * qt + 6; j1++) {
                if (j0 == j1) {
                    continue;    // Counterterm: dm2 = diag(SelfEnergy)
                }
                for (int j2 = 6 * qt; j2 < 6 * qt + 6; j2++) {
                    double mls[5] = { Params->murs, Params->mSQ[j2], 0., Params->mSQ[j0], Params->mSQ[j1] };
                    struct Coupling DumC = {0., 1.};
                    struct Coupling Cs[2] = { Params->SQSQSQ[j0][j2][j1], DumC };
                    if (FastFlag(Cs, 2) == 1) {
                        continue;
                    }
                    ff->SetLoopMass(mls);
                    ff->SetSCoupling(Cs);

                    if (Params->out1 >= 10) {
                        // No s-channel
                    } else {

#ifdef TT
                        for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][j0], Params->CHSQq[ii][j1][aa],
                                       Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 1) {
                                continue;
                            }
                            ff->SetPropagator(Params->mSQ[i1], Params->mSQ[i1], 0, 0);
                            ff->SetWCoupling(Cw);
                            virt += ff->MVtbu5t(IEPS);
                        }
#endif
#ifdef ST
                        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][j0], Params->CHSQq[ii][j1][aa],
                                       Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mv[i0], Params->mv[i0], Params->Gv[i0], Params->Gv[i0]);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVtbu5s(IEPS);
                            }
                        }
#endif
#ifdef TU
                        if (qu <= 1) {
                            for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                                struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][j0], Params->CHSQq[ii][j1][aa],
                                           Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                                };
                                if (FastFlag(Cw, 4) == 0) {
                                    ff->SetPropagator(Params->mSQ[i1], Params->mSQ[i1], Params->mSQ[i1] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                    ff->SetWCoupling(Cw);
                                    virt += ff->MVtbu5u(IEPS);
                                }
                            }
                        }
#endif
                    }
                }
            }
        }
        if (qu <= 1) {
            for (int j0 = 6 * qu; j0 < 6 * qu + 6; j0++) {
                for (int j1 = 6 * qu; j1 < 6 * qu + 6; j1++) {
                    if (j0 == j1) {
                        continue;    // Counterterm: dm2 = diag(SelfEnergy)
                    }
                    for (int j2 = 6 * qu; j2 < 6 * qu + 6; j2++) {
                        double mls[5] = { Params->murs, Params->mSQ[j2], 0., Params->mSQ[j0], Params->mSQ[j1] };
                        struct Coupling DumC = {0., 1.};
                        struct Coupling Cs[2] = { Params->SQSQSQ[j0][j2][j1], DumC };
                        if (FastFlag(Cs, 2) == 1) {
                            continue;
                        }
                        ff->SetLoopMass(mls);
                        ff->SetSCoupling(Cs);
#ifdef UU
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[ii][bb][j0], Params->CHSQq[jj][j1][aa],
                                       Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[i1], Params->mSQ[i1], Params->mSQ[i1] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVubu5u(IEPS);
                            }
                        }
#endif
#ifdef SU
                        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[ii][bb][j0], Params->CHSQq[jj][j1][aa],
                                       Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mv[i0], Params->mv[i0], Params->Gv[i0], Params->Gv[i0]);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVubu5s(IEPS);
                            }
                        }
#endif
#ifdef TU
                        if (qt >= 0) {
                            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                                struct Coupling Cw[4] = { Params->CHqSQ[ii][bb][j0], Params->CHSQq[jj][j1][aa],
                                           Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa]
                                };
                                if (FastFlag(Cw, 4) == 0) {
                                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i0], Params->mSQ[i0] * 1.e-2, Params->mSQ[i0] * 1.e-2);
                                    ff->SetWCoupling(Cw);
                                    virt += ff->MVubu5t(IEPS);
                                }
                            }
                        }
#endif
                    }
                }
            }
        }
    }

    delete ff;

    return virt;
}

double Vstr2(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
        for (int j1 = (aa / 3) * 6; j1 < (aa / 3 + 1) * 6; j1++) {
            double mls[5] = { Params->murs, Params->mGL, Params->mSQ[j0], Params->mSQ[j1], 0. };
            struct Coupling Cs[2] = { Params->GLqSQ[bb][j0], Params->GLSQq[j1][aa] };
            if (FastFlag(Cs, 2) == 1) {
                continue;
            }
            ff->SetLoopMass(mls);
            ff->SetSCoupling(Cs);

            if (Params->out1 >= 10) {
#ifdef SS
                for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                        struct Coupling Cw[4] = { Params->vSLSL[i0][ii - 10][jj - 10], Params->vqq[i0][bb][aa],
                                   Params->vSLSL[i1][ii - 10][jj - 10], Params->vqq[i1][bb][aa]
                        };
                        if (FastFlag(Cw, 4) == 0) {

                            ff->SetPropagator(Params->mv[i0], Params->mv[i1] , Params->Gv[i0], Params->Gv[i1]);
                            ff->SetWCoupling(Cw);
#ifdef SQUARKS
                            virt += ff->MVstr2sSL(IEPS);
#endif
#ifdef DEC_SQUARKS
                            virt = +0;
#endif
                        }
                    }
                }
#endif

            } else {

#ifdef SS
                for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                        struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vSQSQ[i0][j0][j1],
                                   Params->vCHCH[i1][ii][jj], Params->vqq[i1][bb][aa]
                        };
                        if (FastFlag(Cw, 4) == 0) {
                            ff->SetPropagator(Params->mv[i0], Params->mv[i1] , Params->Gv[i0], Params->Gv[i1]);
                            ff->SetWCoupling(Cw);
                            virt += ff->MVstr2s(IEPS);
                        }
                    }
                }
#endif
#ifdef ST
                if (qt >= 0) {
                    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                        for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                            struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vSQSQ[i0][j0][j1],
                                       Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mv[i0], Params->mSQ[i1] , Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVstr2t(IEPS);
                            }
                        }
                    }
                }
#endif
#ifdef SU
                if (qu <= 1) {
                    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vSQSQ[i0][j0][j1],
                                       Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mv[i0], Params->mSQ[i1] , Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVstr2u(IEPS);
                            }
                        }
                    }
                }
#endif
            }
        }
    }

    delete ff;

    return virt;
}

double Vstr3(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
        for (int j1 = (aa / 3) * 6; j1 < (aa / 3 + 1) * 6; j1++) {
            double mls[5] = { Params->murs, Params->mGL, Params->mSQ[j0], Params->mSQ[j1], 0. };
            struct Coupling Cs[2] = { Params->GLqSQ[bb][j0], Params->GLSQq[j1][aa] };
            if (FastFlag(Cs, 2) == 1) {
                continue;
            }
            ff->SetLoopMass(mls);
            ff->SetSCoupling(Cs);

            if (Params->out1 >= 10) {
                // No s-channel
            } else {

#ifdef ST
                if (qt >= 0) {
                    for (int i0 = 4 * qs; i0 < 2 * qs + 4; i0++) {
                        for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                            struct Coupling Cw[4] = { Params->hCHCH[i0][ii][jj], Params->hSQSQ[i0][j0][j1],
                                       Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mh[i0], Params->mSQ[i1], Params->mh[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVstr3t(IEPS);
                            }
                        }
                    }
                }
#endif
#ifdef SU
                if (qu <= 1) {
                    for (int i0 = 4 * qs; i0 < 2 * qs + 4; i0++) {
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw[4] = { Params->hCHCH[i0][ii][jj], Params->hSQSQ[i0][j0][j1],
                                       Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mh[i0], Params->mSQ[i1], Params->mh[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVstr3u(IEPS);
                            }
                        }
                    }
                }
#endif
            }
        }
    }

    delete ff;

    return virt;
}

double Vttr3(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }

    if (qt >= 0) {
        for (int j0 = (aa / 3) * 6; j0 < (aa / 3 + 1) * 6; j0++) {
            for (int j1 = qt * 3; j1 < qt * 3 + 3; j1++) {
                for (int j2 = qt * 6; j2 < qt * 6 + 6; j2++) {
                    double mls[5] = { Params->murs, Params->mGL, Params->mSQ[j0], Params->mq[j1], 0. };
                    struct Coupling Cs[2] = { Params->GLSQq[j2][j1], Params->GLSQq[j0][aa] };
                    if (FastFlag(Cs, 2) == 1) {
                        continue;
                    }
                    ff->SetLoopMass(mls);
                    ff->SetSCoupling(Cs);

                    if (Params->out1 >= 10) {
                        // No s-channel
                    } else {

#ifdef TT
                        for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][j2], Params->CHqSQ[ii][j1][j0],
                                       Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[j2], Params->mSQ[i1], Params->mSQ[j2] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVttr3t(IEPS);
                            }
                        }
#endif
#ifdef ST
                        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][j2], Params->CHqSQ[ii][j1][j0],
                                       Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[j2], Params->mv[i0], Params->mSQ[j2] * 1.e-2, Params->Gv[i0]);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVttr3s(IEPS);
                            }
                        }
#endif
#ifdef TU
                        if (qu <= 1) {
                            for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                                struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][j2], Params->CHqSQ[ii][j1][j0],
                                           Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                                };
                                if (FastFlag(Cw, 4) == 0) {
                                    ff->SetPropagator(Params->mSQ[j2], Params->mSQ[i1], Params->mSQ[j2] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                    ff->SetWCoupling(Cw);
                                    virt += ff->MVttr3u(IEPS);
                                }
                            }
                        }
#endif
                    }
                }
            }
        }
    }
    delete ff;

    return virt;
}

double Vttr4(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    if (qt >= 0) {
        for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
            for (int j1 = qt * 3; j1 < qt * 3 + 3; j1++) {
                for (int j2 = qt * 6; j2 < qt * 6 + 6; j2++) {
                    double mls[5] = { Params->murs, Params->mGL, Params->mSQ[j0], Params->mq[j1], 0. };
                    struct Coupling Cs[2] = { Params->GLqSQ[j1][j2], Params->GLqSQ[bb][j0] };
                    if (FastFlag(Cs, 2) == 1) {
                        continue;
                    }
                    ff->SetLoopMass(mls);
                    ff->SetSCoupling(Cs);

                    if (Params->out1 >= 10) {
                        // No s-channel.
                    } else {

#ifdef TT
                        for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHSQq[jj][j0][j1], Params->CHSQq[ii][j2][aa],
                                       Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[j2], Params->mSQ[i1], Params->mSQ[j2] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVttr4t(IEPS);
                            }
                        }
#endif
#ifdef ST
                        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                            struct Coupling Cw[4] = { Params->CHSQq[jj][j0][j1], Params->CHSQq[ii][j2][aa],
                                       Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[j2], Params->mv[i0], Params->mSQ[j2] * 1.e-2, Params->Gv[i0]);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVttr4s(IEPS);
                            }
                        }
#endif
#ifdef TU
                        if (qu <= 1) {
                            for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                                struct Coupling Cw[4] = { Params->CHSQq[jj][j0][j1], Params->CHSQq[ii][j2][aa],
                                           Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                                };
                                if (FastFlag(Cw, 4) == 0) {
                                    ff->SetPropagator(Params->mSQ[j2], Params->mSQ[i1], Params->mSQ[j2] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                    ff->SetWCoupling(Cw);
                                    virt += ff->MVttr4u(IEPS);
                                }
                            }
                        }
#endif
                    }
                }
            }
        }
    }
    delete ff;

    return virt;
}

double Vutr3(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    if (qt >= 0) {
        for (int j0 = (aa / 3) * 6; j0 < (aa / 3 + 1) * 6; j0++) {
            for (int j1 = qu * 3; j1 < qu * 3 + 3; j1++) {
                for (int j2 = qu * 6; j2 < qu * 6 + 6; j2++) {
                    double mls[5] = { Params->murs, Params->mGL, Params->mSQ[j0], Params->mq[j1], 0. };
                    struct Coupling Cs[2] = { Params->GLSQq[j2][j1], Params->GLSQq[j0][aa] };
                    if (FastFlag(Cs, 2) == 1) {
                        continue;
                    }
                    ff->SetLoopMass(mls);
                    ff->SetSCoupling(Cs);

                    if (Params->out1 >= 10) {
                        // No s-channel.
                    } else {

#ifdef UU
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[ii][bb][j2], Params->CHqSQ[jj][j1][j0],
                                       Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[j2], Params->mSQ[i1], Params->mSQ[j2] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVutr3u(IEPS);
                            }
                        }
#endif
#ifdef SU
                        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[ii][bb][j2], Params->CHqSQ[jj][j1][j0],
                                       Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[j2], Params->mv[i0], Params->mSQ[j2] * 1.e-2, Params->Gv[i0]);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVutr3s(IEPS);
                            }
                        }
#endif
#ifdef TU
                        if (qt >= 0) {
                            for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                                struct Coupling Cw[4] = { Params->CHqSQ[ii][bb][j2], Params->CHqSQ[jj][j1][j0],
                                           Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                                };
                                if (FastFlag(Cw, 4) == 0) {
                                    ff->SetPropagator(Params->mSQ[j2], Params->mSQ[i1], Params->mSQ[j2] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                    ff->SetWCoupling(Cw);
                                    virt += ff->MVutr3t(IEPS);
                                }
                            }
                        }
#endif
                    }
                }
            }
        }
    }

    delete ff;

    return virt;
}

double Vutr4(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }

    if (qu <= 1) {
        for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
            for (int j1 = qu * 3; j1 < qu * 3 + 3; j1++) {
                for (int j2 = qu * 6; j2 < qu * 6 + 6; j2++) {
                    double mls[5] = { Params->murs, Params->mGL, Params->mSQ[j0], Params->mq[j1], 0. };
                    struct Coupling Cs[2] = { Params->GLqSQ[j1][j2], Params->GLqSQ[bb][j0] };
                    if (FastFlag(Cs, 2) == 1) {
                        continue;
                    }
                    ff->SetLoopMass(mls);
                    ff->SetSCoupling(Cs);

                    if (Params->out1 >= 10) {
                        // No s-channel.
                    } else {

#ifdef UU
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHSQq[ii][j0][j1], Params->CHSQq[jj][j2][aa],
                                       Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[j2], Params->mSQ[i1], Params->mSQ[j2] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVutr4u(IEPS);
                            }
                        }
#endif
#ifdef SU
                        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                            struct Coupling Cw[4] = { Params->CHSQq[ii][j0][j1], Params->CHSQq[jj][j2][aa],
                                       Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[j2], Params->mv[i0], Params->mSQ[j2] * 1.e-2, Params->Gv[i0]);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVutr4s(IEPS);
                            }
                        }
#endif
#ifdef TU
                        if (qt >= 0) {
                            for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                                struct Coupling Cw[4] = { Params->CHSQq[ii][j0][j1], Params->CHSQq[jj][j2][aa],
                                           Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                                };
                                if (FastFlag(Cw, 4) == 0) {
                                    ff->SetPropagator(Params->mSQ[j2], Params->mSQ[i1], Params->mSQ[j2] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                    ff->SetWCoupling(Cw);
                                    virt += ff->MVutr4t(IEPS);
                                }
                            }
                        }
#endif
                    }
                }
            }
        }
    }
    delete ff;

    return virt;
}

double Vtbo2(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    if (qt >= 0) {
        for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
            for (int j1 = (aa / 3) * 6; j1 < (aa / 3 + 1) * 6; j1++) {
                for (int j2 = qt * 3; j2 < 3 * qt + 3; j2++) {
                    double mls[5] = { Params->murs, Params->mSQ[j0], Params->mGL,
                                      Params->mSQ[j1], Params->mq[j2]
                                    };
                    struct Coupling Cs[2] = { Params->GLqSQ[bb][j0], Params->GLSQq[j1][aa] };
                    if (FastFlag(Cs, 2) == 1) {
                        continue;
                    }
                    ff->SetLoopMass(mls);
                    ff->SetSCoupling(Cs);

                    if (Params->out1 >= 10) {
                        // No s-channel.
                    } else {

#ifdef TT
                        for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[ii][j2][j1], Params->CHSQq[jj][j0][j2],
                                       Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[i1], Params->mSQ[i1], Params->mSQ[i1] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVtbo2t(IEPS);
                            }
                        }
#endif
#ifdef ST
                        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                            struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                                       Params->CHqSQ[ii][j2][j1], Params->CHSQq[jj][j0][j2]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mv[i0], Params->mv[i0], Params->Gv[i0], Params->Gv[i0]);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVtbo2s(IEPS);
                            }
                        }
#endif
#ifdef TU
                        if (qu <= 1) {
                            for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                                struct Coupling Cw[4] = { Params->CHqSQ[ii][j2][j1], Params->CHSQq[jj][j0][j2],
                                           Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                                };
                                if (FastFlag(Cw, 4) == 0) {
                                    ff->SetPropagator(Params->mSQ[i1], Params->mSQ[i1], Params->mSQ[i1] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                    ff->SetWCoupling(Cw);
                                    virt += ff->MVtbo2u(IEPS);
                                }
                            }
                        }
#endif
                    }
                }
            }
        }
    }
    delete ff;

    return virt;
}

double Vubo2(const double S, const double T, Parameters *Params)
{
    double virt = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    if (qu >= 0) {
        for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
            for (int j1 = (aa / 3) * 6; j1 < (aa / 3 + 1) * 6; j1++) {
                for (int j2 = qu * 3; j2 < 3 * qu + 3; j2++) {
                    double mls[5] = { Params->murs, Params->mSQ[j0], Params->mGL,
                                      Params->mSQ[j1], Params->mq[j2]
                                    };
                    struct Coupling Cs[2] = { Params->GLqSQ[bb][j0], Params->GLSQq[j1][aa] };
                    if (FastFlag(Cs, 2) == 1) {
                        continue;
                    }
                    ff->SetLoopMass(mls);
                    ff->SetSCoupling(Cs);

                    if (Params->out1 >= 10) {
                        // No s-channel.
                    } else {

#ifdef UU
                        for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                            struct Coupling Cw[4] = { Params->CHqSQ[jj][j2][j1], Params->CHSQq[ii][j0][j2],
                                       Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mSQ[i1], Params->mSQ[i1], Params->mSQ[i1] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVubo2u(IEPS);
                            }
                        }
#endif
#ifdef SU
                        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                            struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                                       Params->CHqSQ[jj][j2][j1], Params->CHSQq[ii][j0][j2]
                            };
                            if (FastFlag(Cw, 4) == 0) {
                                ff->SetPropagator(Params->mv[i0], Params->mv[i0], Params->Gv[i0], Params->Gv[i0]);
                                ff->SetWCoupling(Cw);
                                virt += ff->MVubo2s(IEPS);
                            }
                        }
#endif
#ifdef TU
                        if (qt >= 0) {
                            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                                struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                                           Params->CHqSQ[jj][j2][j1], Params->CHSQq[ii][j0][j2]
                                };
                                if (FastFlag(Cw, 4) == 0) {
                                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i0], Params->mSQ[i0] * 1.e-2, Params->mSQ[i0] * 1.e-2);
                                    ff->SetWCoupling(Cw);
                                    virt += ff->MVubo2t(IEPS);
                                }
                            }
                        }
#endif
                    }
                }
            }
        }
    }
    delete ff;

    return virt;
}

double Ctt3(const double S, const double T, Parameters *Params)
{
    double born = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    if (Params->out1 >= 10) {
        // No s-channel.
    } else {

#ifdef TT
        if (qt >= 0 && IEPS == 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->dCHSQq[ii][i0][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBtt();
                }
            }
        }
#endif
#ifdef ST
        if (qt >= 0) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[jj][bb][i1], Params->dCHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1] , Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    if (IEPS == 1) {
                        born += ff->MBst();
                    } else if (IEPS == 0) {
                        born += ff->MB1st();
                    }
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1 && IEPS == 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->dCHSQq[ii][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBtu();
                }
            }
        }
#endif
    }

    delete ff;

    return born;
}

double Ctt4(const double S, const double T, Parameters *Params)
{
    double born = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    if (Params->out1 >= 10) {
        // No s-channel.
    } else {

#ifdef TT
        if (qt >= 0 && IEPS == 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw[4] = { Params->dCHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBtt();
                }
            }
        }
#endif
#ifdef ST
        if (qt >= 0) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->dCHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1] , Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    if (IEPS == 1) {
                        born += ff->MBst();
                    } else if (IEPS == 0) {
                        born += ff->MB1st();
                    }
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1 && IEPS == 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->dCHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBtu();
                }
            }
        }
#endif
    }

    delete ff;

    return born;
}

double Ctu3(const double S, const double T, Parameters *Params)
{
    double born = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    if (Params->out1 >= 10) {
        // No s-channel.
    } else {

#ifdef UU
        if (qu <= 1 && IEPS == 1) {
            for (int i0 = 6 * qu; i0 < 6 * qu + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[ii][bb][i0], Params->dCHSQq[jj][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBuu();
                }
            }
        }
#endif
#ifdef SU
        if (qu <= 1) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[ii][bb][i1], Params->dCHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1] , Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    if (IEPS == 1) {
                        born += ff->MBsu();
                    } else if (IEPS == 0) {
                        born += ff->MB1su();
                    }
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1 && IEPS == 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->dCHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBtu();
                }
            }
        }
#endif
    }

    delete ff;

    return born;
}

double Ctu4(const double S, const double T, Parameters *Params)
{
    double born = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }

    if (Params->out1 >= 10) {
        // No s-channel.
    } else {

#ifdef UU
        if (qu <= 1 && IEPS == 1) {
            for (int i0 = 6 * qu; i0 < 6 * qu + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->dCHqSQ[ii][bb][i0], Params->CHSQq[jj][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBuu();
                }
            }
        }
#endif
#ifdef SU
        if (qu <= 1) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->dCHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1] , Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    if (IEPS == 1) {
                        born += ff->MBsu();
                    } else if (IEPS == 0) {
                        born += ff->MB1su();
                    }
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1 && IEPS == 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->dCHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1] , Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBtu();
                }
            }
        }
#endif
    }

    delete ff;

    return born;
}

double Virt(const double S, const double T, Parameters *Params)
{
    const double g3s = std::norm(Params->gqq[0][0].R);

    double result = 0.;
    result += Vqcd(S, T, Params); // Include bu3, tr1, tr2(t, u) and bo1(t, u)
    result += .5 * Vbu1(S, T, Params) + .5 * Vbu2(S, T, Params) + Vbu4(S, T, Params) + Vbu5(S, T, Params);
    result += Vstr2(S, T, Params) + (Vstr3(S, T, Params) + Vttr3(S, T, Params) +
                                     Vttr4(S, T, Params) + Vutr3(S, T, Params) + Vutr4(S, T, Params));
    result += (Vtbo2(S, T, Params) + Vubo2(S, T, Params));
    result += (Ctt3(S, T, Params) + Ctt4(S, T, Params) + Ctu3(S, T, Params) + Ctu4(S, T, Params));

    return result / g3s;
}

double Virt2(const double S, const double T, Parameters *Params)
{
    const double g3s = std::norm(Params->gqq[0][0].R);

    double result = 0.;
    result += Vqcd(S, T, Params); // Includes bu3, tr1, tr2(t, u) and bo1(t, u)

    return result / g3s;
}

double DipI(const double S, const double T, Parameters *Params)
{
    const double lnmusq = std::log(Params->murs / S);
    if (IEPS == 2) {
        return 2.*Born(S, T, Params);
    } else if (IEPS == 1) {
        return ((3. + 2.*lnmusq) * Born(S, T, Params) + 2.*Born1(S, T, Params));
    } else if (IEPS == 0) return ((3. + lnmusq) * lnmusq) * Born(S, T, Params)
                                     + (3 + 2.*lnmusq) * Born1(S, T, Params) + 2.*Born2(S, T, Params);
    else {
        return 0.;
    }
}
