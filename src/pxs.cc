// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Partonic cross section.

#include <cmath>
#include <complex>
#include <iostream>
#include <cstdlib>

#include "prm.h"
#include "fi.h"

#define ONSUB

#define SS
#define TT
#define UU
#define ST
#define SU
#define TU

inline bool FastFlag(struct Coupling *C, int N)
{
    for (int i0 = 0; i0 < N; i0++) {
        if (C[i0].L == 0. && C[i0].R == 0.) {
            return 1;
        }
    }
    return 0;
}

inline int iabs(int x)
{
    return x < 0 ? -x : x;
}

double Born(const double S, const double T, Parameters *Params)
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

    if (ii >= 10 && ii < 20) {
#ifdef SS
        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                struct Coupling Cw[4] = { Params->vSLSL[i0][ii - 10][jj - 10], Params->vqq[i0][bb][aa],
                           Params->vSLSL[i1][ii - 10][jj - 10], Params->vqq[i1][bb][aa]
                };
                if (FastFlag(Cw, 4) == 1) {
                    continue;
                }
                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);

                ff->SetWCoupling(Cw);
                born += ff->MBssSL();
            }
        }
#endif
    } else {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
#ifdef SS
        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                           Params->vCHCH[i1][ii][jj], Params->vqq[i1][bb][aa]
                };
                if (FastFlag(Cw, 4) == 1) {
                    continue;
                }
                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MBss();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBtt();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBuu();
                }
            }
        }
#endif
#ifdef ST
        if (qt >= 0) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], 0);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBst();
                }
            }
        }
#endif
#ifdef SU
        if (qu <= 1) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBsu();
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBtu();
                }
            }
        }
#endif
    }

    delete ff;

    return born;
}

double Born1(const double S, const double T, Parameters *Params)
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

    if (Params->out1 < 10) {

#ifdef SS
        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                           Params->vCHCH[i1][ii][jj], Params->vqq[i1][bb][aa]
                };
                if (FastFlag(Cw, 4) == 1) {
                    continue;
                }
                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MB1ss();
            }
        }
#endif
#ifdef ST
        if (qt >= 0) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MB1st();
                }
            }
        }
#endif
#ifdef SU
        if (qu <= 1) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MB1su();
                }
            }
        }
#endif
    }

    delete ff;

    return born;
}

double Born2(const double S, const double T, Parameters *Params)
{
    double born = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel

    FI *ff = new FI();
    //ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, T);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, T);
    }


    if (Params->out1 < 10) {

#ifdef SS
        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
                struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                           Params->vCHCH[i1][ii][jj], Params->vqq[i1][bb][aa]
                };
                if (FastFlag(Cw, 4) == 1) {
                    continue;
                }
                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MB2ss();
            }
        }
#endif

    }
    delete ff;

    return born;
}

double RealG(const double S, const double M2, const double PT2,
             const double TH, const double PH, const int YS, Parameters *Params)
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
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, M2, PT2, TH, PH, YS);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, M2, PT2, TH, PH, YS);
    }


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

                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MGssSL();
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
                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MGss();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2 , Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MGtt();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MGuu();
                }
            }
        }
#endif
#ifdef ST
        if (qt >= 0) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MGst();
                }
            }
        }
#endif
#ifdef SU
        if (qu <= 1) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MGsu();
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MGtu();
                }
            }
        }
#endif

    }

    delete ff;

    return born;
}

double RealQ(const double S, const double M2, const double PT2,
             const double TH, const double PH, const int YS, Parameters *Params)
{
    double born = 0.;

    const int aa = Params->in1;
    const int bb = Params->in2;
    const int ii = Params->out1;
    const int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, M2, PT2, TH, PH, YS);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, M2, PT2, TH, PH, YS);
    }


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

                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MQssSL();
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
                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MQss();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MQtt();
#ifdef ONSUB
                    if (i0 == i1) {
                        born -= ff->MQttp(); // On-shell substraction
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MQuu();
#ifdef ONSUB
                    if (i0 == i1) {
                        born -= ff->MQuup(); // On-shell substraction
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
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MQst();
                }
            }
        }
#endif
#ifdef SU
        if (qu <= 1) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MQsu();
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MQtu();
                }
            }
        }
#endif
    }

    delete ff;

    return born;
}

double RealQB(const double S, const double M2, const double PT2,
              const double TH, const double PH, const int YS, Parameters *Params)
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
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, M2, PT2, TH, PH, YS);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, M2, PT2, TH, PH, YS);
    }


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

                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MQBssSL();
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
                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MQBss();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MQBtt();
#ifdef ONSUB
                    if (i0 == i1) {
                        born -= ff->MQBttp(); // On-shell substraction
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MQBuu();
#ifdef ONSUB
                    if (i0 == i1) {
                        born -= ff->MQBuup(); // On-shell substraction
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
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MQBst();
                }
            }
        }
#endif
#ifdef SU
        if (qu <= 1) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MQBsu();
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MQBtu();
                }
            }
        }
#endif
    }

    delete ff;

    return born;
}

double DipGA(const double S, const double M2, const double PT2,
             const double TH, const double PH, const int YS, Parameters *Params)
{
    double born = 0.;
    double x = 0.;
    double fact = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, M2, PT2, TH, PH, YS);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, M2, PT2, TH, PH, YS);
    }

    ff->SetDipKinematicA(x, fact);

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

                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MBssSL();
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
                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MBss();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBtt();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBuu();
                }
            }
        }
#endif
#ifdef ST
        if (qt >= 0) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBst();
                }
            }
        }
#endif
#ifdef SU
        if (qu <= 1) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBsu();
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBtu();
                }
            }
        }
#endif
    }

    delete ff;

    return fact * (2. / (1. - x) - 1 - x) * born;
}

double DipGB(const double S, const double M2, const double PT2,
             const double TH, const double PH, const int YS, Parameters *Params)
{
    double born = 0.;
    double x = 0.;
    double fact = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, M2, PT2, TH, PH, YS);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, M2, PT2, TH, PH, YS);
    }


    ff->SetDipKinematicB(x, fact);

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
                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MBssSL();
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
                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MBss();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBtt();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBuu();
                }
            }
        }
#endif
#ifdef ST
        if (qt >= 0) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBst();
                }
            }
        }
#endif
#ifdef SU
        if (qu <= 1) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBsu();
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBtu();
                }
            }
        }
#endif
    }

    delete ff;

    return fact * (2. / (1. - x) - 1 - x) * born;
}
double DipQA(const double S, const double M2, const double PT2,
             const double TH, const double PH, const int YS, Parameters *Params)
{
    double born = 0.;
    double x = 0.;
    double fact = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, M2, PT2, TH, PH, YS);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, M2, PT2, TH, PH, YS);
    }


    ff->SetDipKinematicA(x, fact);


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

                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MBssSL();
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
                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MBss();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBtt();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBuu();
                }
            }
        }
#endif
#ifdef ST
        if (qt >= 0) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBst();
                }
            }
        }
#endif
#ifdef SU
        if (qu <= 1) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBsu();
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBtu();
                }
            }
        }
#endif
    }

    delete ff;

    return fact * (1. - 2.*x * (1 - x)) * born;
}

double DipQB(const double S, const double M2, const double PT2,
             const double TH, const double PH, const int YS, Parameters *Params)
{
    double born = 0.;
    double x = 0.;
    double fact = 0.;

    int aa = Params->in1;
    int bb = Params->in2;
    int ii = Params->out1;
    int jj = Params->out2;

    int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
    int qt = aa / 3 - ii / 4; // Squark charge for t-channel
    int qu = aa / 3 + jj / 4; // Squark charge for u-channel

    FI *ff = new FI();

    if (ii < 10) {
        ff->SetKinematic(Params->mCH[ii], Params->mCH[jj], S, M2, PT2, TH, PH, YS);
    } else if (ii >= 10 && ii < 20) {
        ff->SetKinematic(Params->mSL[ii - 10], Params->mSL[jj - 10], S, M2, PT2, TH, PH, YS);
    }


    ff->SetDipKinematicB(x, fact);

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

                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MBssSL();

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
                ff->SetPropagator(Params->mv[i0], Params->mv[i1], Params->Gv[i0], Params->Gv[i1]);
                ff->SetWCoupling(Cw);
                born += ff->MBss();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBtt();
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
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += ff->MBuu();
                }
            }
        }
#endif
#ifdef ST
        if (qt >= 0) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[jj][bb][i1], Params->CHSQq[ii][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBst();
                }
            }
        }
#endif
#ifdef SU
        if (qu <= 1) {
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->vCHCH[i0][ii][jj], Params->vqq[i0][bb][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mv[i0], Params->mSQ[i1], Params->Gv[i0], Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBsu();
                }
            }
        }
#endif
#ifdef TU
        if (qt >= 0 && qu <= 1) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
                for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                    struct Coupling Cw[4] = { Params->CHqSQ[jj][bb][i0], Params->CHSQq[ii][i0][aa],
                               Params->CHqSQ[ii][bb][i1], Params->CHSQq[jj][i1][aa]
                    };
                    if (FastFlag(Cw, 4) == 1) {
                        continue;
                    }
                    ff->SetPropagator(Params->mSQ[i0], Params->mSQ[i1], Params->mSQ[i0] * 1.e-2, Params->mSQ[i1] * 1.e-2);
                    ff->SetWCoupling(Cw);
                    born += 2.*ff->MBtu();
                }
            }
        }
#endif
    }

    delete ff;

    return fact * (1. - 2.*x * (1. - x)) * born;
}
