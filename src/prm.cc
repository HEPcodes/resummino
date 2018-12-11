// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Reads and sets the parameters.

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <cmath>
#include <string>
#include <string.h>
#include <complex>
#include "slhaea.h"
#include <map>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "prm.h"
#include "utils.h"

using namespace std;

#define SUSY

#define II complex<double>(0.0, 1.0)

map<string, string> Parameters::read_input_file(string filename)
{
    map<string, string> config;
    FILE *input_file;
    if (filename == "-") {
        input_file = stdin;
    } else {
        input_file = fopen(filename.c_str(), "r");
    }
    if (!input_file) {
        fprintf(stderr, "error: could not open file '%s'.\n", filename.c_str());
        exit(1);
    }


    char *variable = (char*)malloc(128 + 1);
    char *value = (char*)malloc(128 + 1);;
    for (;;) {
        char *input;
        size_t n;
        int chars_read;
        int elements_read;

        input = NULL;
        n = 0;
        chars_read = getline(&input, &n, input_file);
        if (chars_read == -1) {
            break;
        }
        if (chars_read < 2) {
            continue;
        }
        elements_read = sscanf(input, "%128[a-zA-Z0-9_-] = %128s\n",
                               variable, value);
        if (elements_read == 0) {
            // Blank line.
            continue;
        } else if (elements_read == 1 && variable[0] == '#') {
            // Comment.
            continue;
        } else if (elements_read == 1) {
            fprintf(stderr,
                    "error: while reading line '%s' of the input file.\n",
                    variable);
            exit(1);
        } else if (variable[0] == '#') {
            // Comment.
            continue;
        }
        config[variable] = value;
    }
    free(variable);
    free(value);
    fclose(input_file);
    return config;
}

void Parameters::init_couplings()
{
    // Higgs couplings
    for (int i0 = 0; i0 < nh; i0++) {
        for (int i1 = 0; i1 < nq; i1++) {
            for (int i2 = 0; i2 < nq; i2++) {
                hqq[i0][i1][i2].L = complex<double>(0., 0.);
                hqq[i0][i1][i2].R = complex<double>(0., 0.);
            }
        }
    }
    for (int i0 = 0; i0 < nh; i0++) {
        for (int i1 = 0; i1 < nSQ; i1++) {
            for (int i2 = 0; i2 < nSQ; i2++) {
                hSQSQ[i0][i1][i2].L = complex<double>(0., 0.);
                hSQSQ[i0][i1][i2].R = complex<double>(0., 0.);
            }
        }
    }
    for (int i0 = 0; i0 < nh; i0++) {
        for (int i1 = 0; i1 < nCH; i1++) {
            for (int i2 = 0; i2 < nCH; i2++) {
                hCHCH[i0][i1][i2].L = complex<double>(0., 0.);
                hCHCH[i0][i1][i2].R = complex<double>(0., 0.);
            }
        }
    }
    // Vector couplings to quarks.
    for (int i0 = 0; i0 < nv; i0++) {
        for (int i1 = 0; i1 < nq; i1++) {
            for (int i2 = 0; i2 < nq; i2++) {
                vqq[i0][i1][i2].L = complex<double>(0., 0.);
                vqq[i0][i1][i2].R = complex<double>(0., 0.);
            }
        }
    }
    // Vector couplings to squarks.
    for (int i0 = 0; i0 < nv; i0++) {
        for (int i1 = 0; i1 < nSQ; i1++) {
            for (int i2 = 0; i2 < nSQ; i2++) {
                vSQSQ[i0][i1][i2].L = complex<double>(0., 0.);
                vSQSQ[i0][i1][i2].R = complex<double>(0., 0.);
            }
        }
    }
    for (int i0 = 0; i0 < nv; i0++) {
        for (int i1 = 0; i1 < nCH; i1++) {
            for (int i2 = 0; i2 < nCH; i2++) {
                vCHCH[i0][i1][i2].L = complex<double>(0., 0.);
                vCHCH[i0][i1][i2].R = complex<double>(0., 0.);
            }
        }
    }
    // Gaugino couplings
    for (int i0 = 0; i0 < nCH; i0++) {
        for (int i1 = 0; i1 < nSQ; i1++) {
            for (int i2 = 0; i2 < nq; i2++) {
                CHSQq[i0][i1][i2].L = complex<double>(0., 0.);
                CHSQq[i0][i1][i2].R = complex<double>(0., 0.);
                CHqSQ[i0][i2][i1].L = complex<double>(0., 0.);
                CHqSQ[i0][i2][i1].R = complex<double>(0., 0.);

                dCHSQq[i0][i1][i2].L = complex<double>(0., 0.);
                dCHSQq[i0][i1][i2].R = complex<double>(0., 0.);
                dCHqSQ[i0][i2][i1].L = complex<double>(0., 0.);
                dCHqSQ[i0][i2][i1].R = complex<double>(0., 0.);
            }
        }
    }
    // Strong couplings
    for (int i0 = 0; i0 < nq; i0++) {
        for (int i1 = 0; i1 < nq; i1++) {
            gqq[i0][i1].L = complex<double>(0., 0.);
            gqq[i0][i1].R = complex<double>(0., 0.);
        }
    }
    for (int i0 = 0; i0 < nSQ; i0++) {
        for (int i1 = 0; i1 < nSQ; i1++) {
            gSQSQ[i0][i1].L = complex<double>(0., 0.);
            gSQSQ[i0][i1].R = complex<double>(0., 0.);
        }
    }
    for (int i0 = 0; i0 < nSQ; i0++) {
        for (int i1 = 0; i1 < nq; i1++) {
            GLSQq[i0][i1].L = complex<double>(0., 0.);
            GLSQq[i0][i1].R = complex<double>(0., 0.);
            GLqSQ[i1][i0].L = complex<double>(0., 0.);
            GLqSQ[i1][i0].R = complex<double>(0., 0.);
        }
    }
    for (int i0 = 0; i0 < nSQ; i0++) {
        for (int i1 = 0; i1 < nSQ; i1++) {
            for (int i2 = 0; i2 < nSQ; i2++) {
                SQSQSQ[i0][i1][i2].L = complex<double>(0., 0.);
                SQSQSQ[i0][i1][i2].R = complex<double>(0., 0.);
            }
        }
    }

    // Vector coupling to sleptons.
    for (int i0 = 0; i0 < nv; i0++) {
        for (int i1 = 0; i1 < nSL; i1++) {
            for (int i2 = 0; i2 < nSL; i2++) {
                vSLSL[i0][i1][i2].L = complex<double>(0., 0.);
                vSLSL[i0][i1][i2].R = complex<double>(0., 0.);
            }
        }
    }
}

void Parameters::set_couplings()
{
    // Sets the weak angle.
    xw = pow(g1, 2) / (pow(g1, 2) + pow(g2, 2));
    sw = sqrt(xw);
    cw = sqrt(1. - xw);

    // Sets quark charges.
    double eq[6], T3q[6];
    for (int i0 = 0; i0 < 3; i0++) {
        eq[i0] = -1.0 / 3.0;
        eq[i0 + 3] = 2.0 / 3.0;
        T3q[i0] = -0.5;
        T3q[i0 + 3] = 0.5;
    }

    init_couplings();

    // Sets electroweak couplings.
    // Notation: [0] = photon, [1] = Z, [2] = W.
    for (int i0 = 0; i0 < 6; i0++) {
        vqq[0][i0][i0].L = -g2 * sw * eq[i0];
        vqq[0][i0][i0].R = vqq[0][i0][i0].L;
        vqq[1][i0][i0].L = -g2 / cw * (T3q[i0] - eq[i0] * xw);
        vqq[1][i0][i0].R =  g2 / cw * eq[i0] * xw;
    }
    for (int i0 = 0; i0 < 3; i0++) {
        for (int i1 = 0; i1 < 3; i1++) {
            vqq[2][i0 + 3][i1].L = -g2 * M_SQRT1_2 * ckm[i0][i1];
            vqq[2][i0][i1 + 3].L = -g2 * M_SQRT1_2 * conj(ckm[i1][i0]);
        }
    }

    // Sets hard-coded widths.
    Gv[0] = 0.0; // y width
    Gv[1] = 2.41143316; // Z0 width
    Gv[2] = 2.00282196; // W width

#ifdef SUSY
    // Sets SUSY couplings.
    for (int i0 = 0; i0 < 12; i0++) {
        vSQSQ[0][i0][i0].L = vqq[0][i0 / 2][i0 / 2].L;
        vSQSQ[0][i0][i0].R = vqq[0][i0 / 2][i0 / 2].R;
    }

    for (int i0 = 0; i0 < 6; i0++) {
        for (int i1 = 0; i1 < 6; i1++) {
            for (int i2 = 0; i2 < 3; i2++) {
                vSQSQ[1][i0][i1].R +=
                    vqq[1][i0 / 2][i1 / 2].L * RSD[i0][i2] * conj(RSD[i1][i2]) +
                    vqq[1][i0 / 2][i1 / 2].R * RSD[i0][i2 + 3] * conj(RSD[i1][i2 + 3]);
                vSQSQ[1][i0 + 6][i1 + 6].R +=
                    vqq[1][i0 / 2 + 3][i1 / 2 + 3].L * RSD[i0][i2] * conj(RSD[i1][i2]) +
                    vqq[1][i0 / 2 + 3][i1 / 2 + 3].R * RSD[i0][i2 + 3] * conj(RSD[i1][i2 + 3]);
                for (int i3 = 0; i3 < 3; i3++) {
                    vSQSQ[2][i0 + 6][i1].R += vqq[2][i2 + 3][i3].L
                                              * RSU[i0][i2] * conj(RSD[i1][i3]);
                    vSQSQ[2][i1][i0 + 6].R += vqq[2][i2][i3 + 3].L
                                              * conj(RSU[i0][i2]) * RSD[i1][i3];
                }
            }
        }
    }

    for (int i0 = 0; i0 < 4; i0++) {
        for (int i1 = 0; i1 < 4; i1++) {
            vCHCH[1][i0][i1].L = g2 / cw * .5 *
                                 (-RN[i0][2] * conj(RN[i1][2]) + RN[i0][3] * conj(RN[i1][3]));
            vCHCH[1][i1][i0].R = -vCHCH[1][i0][i1].L;
        }
    }

    for (int i0 = 0; i0 < 4; i0++) {
        for (int i1 = 0; i1 < 2; i1++) {
            vCHCH[2][i0][i1 + 4].L = g2 * (RN[i0][1] * conj(RV[i1][0]) -
                                           RN[i0][3] * conj(RV[i1][1]) * M_SQRT1_2);
            vCHCH[2][i0][i1 + 4].R = g2 * (conj(RN[i0][1]) * RU[i1][0] +
                                           conj(RN[i0][2]) * RU[i1][1] * M_SQRT1_2);

            for (int i2 = 0; i2 < 3; i2++) {
                vCHCH[i2][i1 + 4][i0].L = conj(vCHCH[i2][i0][i1 + 4].L);
                vCHCH[i2][i1 + 4][i0].R = conj(vCHCH[i2][i0][i1 + 4].R);
            }
        }
    }

    for (int i0 = 0; i0 < 2; i0++) {
        for (int i1 = 0; i1 < 2; i1++) {
            double delta = (i0 == i1 ? 1. : 0.);
            vCHCH[0][i0 + 4][i1 + 4].L = -g2 * sw * delta;
            vCHCH[0][i0 + 4][i1 + 4].R = -g2 * sw * delta;
            vCHCH[1][i0 + 4][i1 + 4].L = g2 / cw * (-RV[i0][0] * conj(RV[i1][0])
                                                    - .5 * RV[i0][1] * conj(RV[i1][1])
                                                    + delta * xw);
            vCHCH[1][i0 + 4][i1 + 4].R = g2 / cw * (-conj(RU[i0][0]) * RU[i1][0]
                                                    - .5 * conj(RU[i0][1]) * RU[i1][1]
                                                    + delta * xw);
        }
    }

    // Slepton couplings.
    for (int i0 = 3; i0 < 9; i0++) {
        vSLSL[0][i0][i0].L = g2 * sw;
        vSLSL[0][i0][i0].R = vSLSL[0][i0][i0].L;
    }

    for (int i0 = 3; i0 < 9; i0++) {
        if (i0 % 2 != 0) {
            vSLSL[1][i0][i0].L = g2 / (2.0 * cw) * (1 - 2.0 * xw); // L-Type
            vSLSL[1][i0][i0].R = vSLSL[1][i0][i0].L;
        }

        if (i0 % 2 == 0) {
            vSLSL[1][i0][i0].L = -g2 / (2.0 * cw) * 2.0 * xw; // R-Type
            vSLSL[1][i0][i0].R = vSLSL[1][i0][i0].L;
        }
    }

    for (int i0 = 7; i0 < 9; i0++) {
        for (int i1 = 7; i1 < 9; i1++) {
            vSLSL[1][i0][i1].L = -g2 / (2.0 * cw) * (2.0 * xw * SSLM[i0 - 7][1] * SSLM[i1 - 7][1]
                                 - SSLM[i0 - 7][0] * SSLM[i1 - 7][0] * (1.0 - 2.0 * xw));
            vSLSL[1][i0][i1].R = vSLSL[1][i0][i1].L;
        }
    }

    for (int i0 = 0; i0 < 3; i0++) {
        vSLSL[1][i0][i0].L = -g2 / (2 * cw);
        vSLSL[1][i0][i0].R = vSLSL[1][i0][i0].L;
    }

    for (int i0 = 0; i0 < 2; i0++) {
        vSLSL[2][i0][i0 + (i0 + 3)].L = g2 * (1. / sqrt(2));
        vSLSL[2][i0][i0 + (i0 + 3)].R = g2 * (1. / sqrt(2));
        vSLSL[2][i0 + (i0 + 3)][i0].L = g2 * (1. / sqrt(2));
        vSLSL[2][i0 + (i0 + 3)][i0].R = g2 * (1. / sqrt(2));
    }

    vSLSL[2][7][2].L = g2 * (1. / sqrt(2)) * SSLM[0][1];
    vSLSL[2][7][2].R = g2 * (1. / sqrt(2)) * SSLM[0][1];
    vSLSL[2][8][2].L = g2 * (1. / sqrt(2)) * SSLM[1][1];
    vSLSL[2][8][2].R = g2 * (1. / sqrt(2)) * SSLM[1][1];
    vSLSL[2][2][7].L = g2 * (1. / sqrt(2)) * SSLM[0][1];
    vSLSL[2][2][7].R = g2 * (1. / sqrt(2)) * SSLM[0][1];
    vSLSL[2][2][8].L = g2 * (1. / sqrt(2)) * SSLM[1][1];
    vSLSL[2][2][8].R = g2 * (1. / sqrt(2)) * SSLM[1][1];

    for (int i0 = 0; i0 < 2; i0++) {
        for (int i1 = 0; i1 < 6; i1++) {
            for (int i2 = 0; i2 < 6; i2++) {
                for (int i3 = 0; i3 < 3; i3++) {
                    hSQSQ[i0][i1][i2].R += mv[1] * (RA[i0][0] * RB[0][1] - RA[i0][1] * RB[0][0]) *
                                           (vqq[1][i1 / 2][i2 / 2].L * RSD[i1][i3]  * conj(RSD[i2][i3]) -
                                            vqq[1][i1 / 2][i2 / 2].R * RSD[i1][i3 + 3] * conj(RSD[i2][i3 + 3]));
                    hSQSQ[i0][i1 + 6][i2 + 6].R += mv[1] * (RA[i0][0] * RB[0][1] - RA[i0][1] * RB[0][0]) *
                                                   (vqq[1][i1 / 2 + 3][i2 / 2 + 3].L * RSU[i1][i3]  * conj(RSU[i2][i3]) -
                                                    vqq[1][i1 / 2 + 3][i2 / 2 + 3].R * RSU[i1][i3 + 3] * conj(RSU[i2][i3 + 3]));
                }
            }
        }
    }

    for (int i0 = 0; i0 < 2; i0++) {
        for (int i1 = 0; i1 < 6; i1++) {
            for (int i2 = 0; i2 < 6; i2++) {
                for (int i3 = 0; i3 < 3; i3++) {
                    for (int i4 = 0; i4 < 3; i4++) {
                        for (int i5 = 0; i5 < 2; i5++) {
                            hSQSQ[i0 + 4][i1 + 6][i2].R += mv[2] * vqq[2][i3 + 3][i4].L *
                                                           RSU[i1][i3] * RSD[i2][i4] * RB[0][1 - i5] * RB[1][i5];
                        }
                    }
                }
                hSQSQ[i0 + 4][i2][i1 + 6].R = conj(hSQSQ[i0 + 4][i1 + 6][i2].R);
            }
        }
    }

    for (int i0 = 0; i0 < 2; i0++) {
        for (int i1 = 0; i1 < 4; i1++) {
            for (int i2 = 0; i2 < 4; i2++) {
                hCHCH[i0][i1][i2].R = .5 * g2 / cw *
                                      ((RA[i0][0] * RN[i1][2] - RA[i0][1] * RN[i1][3]) *
                                       (RN[i2][0] * sw - RN[i2][1] * cw) +
                                       (RA[i0][0] * RN[i2][2] - RA[i0][1] * RN[i2][3]) *
                                       (RN[i1][0] * sw - RN[i1][1] * cw));
                hCHCH[i0][i2][i1].L = conj(hCHCH[i0][i1][i2].R);
                hCHCH[i0 + 2][i1][i2].R = .5 * g2 / cw * II *
                                          ((RB[i0][0] * RN[i1][2] - RB[i0][1] * RN[i1][3]) *
                                           (RN[i2][0] * sw - RN[i2][1] * cw) +
                                           (RB[i0][0] * RN[i2][2] - RB[i0][1] * RN[i2][3]) *
                                           (RN[i1][0] * sw - RN[i1][1] * cw));
                hCHCH[i0 + 2][i2][i1].L = conj(hCHCH[i0 + 2][i1][i2].L);
            }
        }
    }

    for (int i0 = 0; i0 < 2; i0++) {
        for (int i1 = 0; i1 < 2; i1++) {
            for (int i2 = 0; i2 < 2; i2++) {
                hCHCH[i0][i1 + 4][i2 + 4].R = -g2 * M_SQRT1_2 *
                                              (RA[i0][0] * RV[i1][0] * RU[i2][1] +
                                               RA[i0][1] * RV[i1][1] * RU[i2][0]);
                hCHCH[i0][i2 + 4][i1 + 4].L = conj(hCHCH[i0][i1 + 4][i2 + 4].R);
                hCHCH[i0 + 2][i1 + 4][i2 + 4].R = -g2 * M_SQRT1_2 * II *
                                                  (RB[i0][0] * RV[i1][0] * RU[i2][1] +
                                                   RB[i0][1] * RV[i1][1] * RU[i2][0]);
                hCHCH[i0 + 2][i2 + 4][i1 + 4].L = conj(hCHCH[i0 + 2][i1 + 4][i2 + 4].R);
            }
        }
    }

    for (int i0 = 0; i0 < 2; i0++) {
        for (int i1 = 0; i1 < 2; i1++) {
            for (int i2 = 0; i2 < 4; i2++) {
                hCHCH[i0 + 4][i1 + 4][i2].L = g2 * RB[i0][0] *
                                              conj(M_SQRT1_2 * RU[i1][1] * (RN[i2][0] * sw / cw + RN[i2][1])
                                                   - RU[i1][0] * RN[i2][2]);
                hCHCH[i0 + 4][i1 + 4][i2].R = -g2 * RB[i0][1] *
                                              (M_SQRT1_2 * RV[i1][1] * (RN[i2][0] * sw / cw + RN[i2][1])
                                               + RV[i1][0] * RN[i2][3]);
                hCHCH[i0 + 4][i2][i1 + 4].L = conj(hCHCH[i0 + 4][i1 + 4][i2].R);
                hCHCH[i0 + 4][i2][i1 + 4].R = conj(hCHCH[i0 + 4][i1 + 4][i2].L);
            }
        }
    }

    for (int i0 = 0; i0 < 4; i0++) {
        for (int i1 = 0; i1 < 6; i1++) {
            for (int i2 = 0; i2 < 3; i2++) {
                CHSQq[i0][i1][i2].L = -g2 * M_SQRT2 * RSD[i1][i2] *
                                      ((eq[i2] - T3q[i2]) * sw / cw * conj(RN[i0][0])
                                       + T3q[i2] * conj(RN[i0][1]))
                                      - yq[i2] * conj(RN[i0][2]) * RSD[i1][i2 + 3];

                CHSQq[i0][i1][i2].R = g2 * M_SQRT2 * RSD[i1][i2 + 3] * eq[i2] * sw / cw * RN[i0][0]
                                      - yq[i2] * RN[i0][2] * RSD[i1][i2];

                CHSQq[i0][i1 + 6][i2 + 3].L = -g2 * M_SQRT2 * RSU[i1][i2]
                                              * ((eq[i2 + 3] - T3q[i2 + 3]) * sw / cw * conj(RN[i0][0])
                                                 + T3q[i2 + 3] * conj(RN[i0][1]))
                                              - yq[i2 + 3] * conj(RN[i0][3]) * RSU[i1][i2 + 3];

                CHSQq[i0][i1 + 6][i2 + 3].R = g2 * M_SQRT2 * RSU[i1][i2 + 3] * eq[i2 + 3] * sw / cw *
                                              RN[i0][0] - yq[i2 + 3] * RN[i0][3] * RSU[i1][i2];
            }
        }
    }

    for (int i0 = 0; i0 < 2; i0++) {
        for (int i1 = 0; i1 < 6; i1++) {
            for (int i2 = 0; i2 < 3; i2++) {
                for (int i3 = 0; i3 < 3; i3++) {
                    CHSQq[i0 + 4][i1 + 6][i2].L += (-g2 * conj(RV[i0][0]) * RSU[i1][i3]
                                                    + yq[i3 + 3] * conj(RV[i0][1])
                                                    * RSU[i1][i3 + 3]) * ckm[i3][i2];

                    CHSQq[i0 + 4][i1 + 6][i2].R += yq[i2] * RU[i0][1] * RSU[i1][i3] * ckm[i3][i2];

                    CHSQq[i0 + 4][i1][i2 + 3].L += (-g2 * conj(RU[i0][0]) * RSD[i1][i3] + yq[i3]
                                                    * conj(RU[i0][1]) * RSD[i1][i3 + 3]) *
                                                   conj(ckm[i2][i3]);

                    CHSQq[i0 + 4][i1][i2 + 3].R += yq[i2 + 3] * RV[i0][1] * RSD[i1][i3] *
                                                   conj(ckm[i2][i3]);
                }
            }
        }
    }

    for (int i0 = 0; i0 < 6; i0++) {
        for (int i1 = 0; i1 < 12; i1++) {
            for (int i2 = 0; i2 < 6; i2++) {
                CHqSQ[i0][i2][i1].L = conj(CHSQq[i0][i1][i2].R);
                CHqSQ[i0][i2][i1].R = conj(CHSQq[i0][i1][i2].L);
            }
        }
    }
#endif

    // Strong couplings.
    for (int i0 = 0; i0 < 6; i0++) {
        for (int i1 = 0; i1 < 6; i1++) {
            gqq[i0][i1].L = (i0 == i1 ? -g3 : 0.);
            gqq[i0][i1].R = (i0 == i1 ? -g3 : 0.);
        }
    }

#ifdef SUSY
    for (int i0 = 0; i0 < 12; i0++) {
        for (int i1 = 0; i1 < 12; i1++) {
            gSQSQ[i0][i1].L = (i0 == i1 ? -g3 : 0.);
            gSQSQ[i0][i1].R = (i0 == i1 ? -g3 : 0.);
        }
    }

    for (int i0 = 0; i0 < 6; i0++) {
        for (int i1 = 0; i1 < 3; i1++) {
            GLSQq[i0][i1].L = -g3 * M_SQRT2 * RSD[i0][i1];
            GLSQq[i0][i1].R = g3 * M_SQRT2 * RSD[i0][i1 + 3];
            GLSQq[i0][i1 + 3].L = 0.;
            gSQSQ[i0][i1 + 3].R = 0.;
            GLSQq[i0 + 6][i1].L = 0.;
            gSQSQ[i0 + 6][i1].R = 0.;
            GLSQq[i0 + 6][i1 + 3].L = -g3 * M_SQRT2 * RSU[i0][i1];
            GLSQq[i0 + 6][i1 + 3].R = g3 * M_SQRT2 * RSU[i0][i1 + 3];
        }
    }

    for (int i0 = 0; i0 < 12; i0++) {
        for (int i1 = 0; i1 < 6; i1++) {
            GLqSQ[i1][i0].L = conj(GLSQq[i0][i1].R);
            GLqSQ[i1][i0].R = conj(GLSQq[i0][i1].L);
        }
    }

    // Effective vertex for SQ-SQloop-SQ.
    complex<double> SD[6][6], SU[6][6];
    for (int i0 = 0; i0 < 6; i0++) {
        for (int i1 = 0; i1 < 6; i1++) {
            for (int i2 = 0; i2 < 3; i2++) {
                SD[i0][i1] += conj(RSD[i0][i2]) * RSD[i1][i2]
                              - conj(RSD[i0][i2 + 3]) * RSD[i1][i2 + 3];
                SU[i0][i1] += conj(RSU[i0][i2]) * RSU[i1][i2]
                              - conj(RSU[i0][i2 + 3]) * RSU[i1][i2 + 3];
            }
        }
    }

    for (int i0 = 0; i0 < 12; i0++) {
        for (int i1 = 0; i1 < 12; i1++) {
            for (int i2 = 0; i2 < 12; i2++) {
                SQSQSQ[i0][i1][i2].L = complex<double>(0., 0.);
                SQSQSQ[i0][i1][i2].R = complex<double>(0., 0.);
            }
        }
    }

    for (int i0 = 0; i0 < 6; i0++) {
        for (int i1 = 0; i1 < 6; i1++) {
            for (int i2 = 0; i2 < 6; i2++) {
                SQSQSQ[i0][i1][i2].R = -pow(g3, 2) * SD[i1][i0] * conj(SD[i1][i2]);
                SQSQSQ[i0 + 6][i1 + 6][i2 + 6].R = -pow(g3, 2) * SU[i1][i0] *
                                                   conj(SU[i1][i2]);
            }
        }
    }

    // Coupling Counterterms.
    for (int i0 = 0; i0 < 6; i0++) { // CH
        for (int i1 = 0; i1 < 12; i1++) { // SQ
            for (int i2 = 0; i2 < 6; i2++) { // q
                for (int i3 = 0; i3 < 12; i3++) { // SQ
                    for (int i4 = 0; i4 < 12; i4++) { // SQ
                        if (abs(mSQ[i1] - mSQ[i3]) / abs(mSQ[i1] + mSQ[i3]) > 1.e-7) {
                            dCHSQq[i0][i1][i2].L += -pow(mSQ[i4], 2) * SQSQSQ[i1][i4][i3].R
                                                    / (pow(mSQ[i1], 2) - pow(mSQ[i3], 2)) * CHSQq[i0][i3][i2].L;
                            dCHSQq[i0][i1][i2].R += -pow(mSQ[i4], 2) * SQSQSQ[i1][i4][i3].R
                                                    / (pow(mSQ[i1], 2) - pow(mSQ[i3], 2)) * CHSQq[i0][i3][i2].R;
                        }
                    }
                }
            }
        }
    }

    for (int i0 = 0; i0 < 6; i0++) {
        for (int i1 = 0; i1 < 12; i1++) {
            for (int i2 = 0; i2 < 6; i2++) {
                dCHqSQ[i0][i2][i1].L = conj(dCHSQq[i0][i1][i2].R);
                dCHqSQ[i0][i2][i1].R = conj(dCHSQq[i0][i1][i2].L);
            }
        }
    }

#endif

}

void Parameters::read_slha(const char *file)
{
    using namespace std;
    using namespace SLHAea;

    ifstream *ifs;

    if (strcmp(file, "-") == 0) {
        ifs = new ifstream("/dev/stdin");
    } else {
        ifs = new ifstream(file);
    }

    if (!(*ifs).good()) {
        fprintf(stderr, "error: Failed to open '%s'\n", file);
        exit(1);
    }

    Coll input(*ifs);
    try {
        mv[1] = to<double>(input.at("SMINPUTS").at("4").at(1)); // Z mass
        // mv[1] = to<double>(input.at("MASS").at("23").at(1)); // Z mass
        mv[2] = to<double>(input.at("MASS").at("24").at(1)); // W mass
    } catch (...) {
        fprintf(stderr, "warning: No vector boson masses defined.\n");
        // exit(0);
    }
    try {
        mq[5] = to<double>(input.at("SMINPUTS").at("6").at(1)); // top mass
    } catch (...) {
        fprintf(stderr, "warning: No top mass defined.\n");
    }
    try {
        mh[1] = to<double>(input.at("MASS").at("25").at(1)); // h mass
        mh[0] = to<double>(input.at("MASS").at("35").at(1)); // H mass
        mh[2] = to<double>(input.at("MASS").at("36").at(1)); // A mass
        mh[4] = to<double>(input.at("MASS").at("37").at(1)); // H+ mass
    } catch (...) {
        fprintf(stderr, "warning: No Higgs masses defined.\n");
    }
    try {
        mGL = to<double>(input.at("MASS").at("1000021").at(1)); // gluino mass
    } catch (...) {
        fprintf(stderr, "warning: No gluino masses defined.\n");
    }
    // Squark masses
    try {
        mSQ[0] = to<double>(input.at("MASS").at("1000001").at(1));
        mSQ[1] = to<double>(input.at("MASS").at("1000003").at(1));
        mSQ[2] = to<double>(input.at("MASS").at("1000005").at(1));
        mSQ[3] = to<double>(input.at("MASS").at("2000001").at(1));
        mSQ[4] = to<double>(input.at("MASS").at("2000003").at(1));
        mSQ[5] = to<double>(input.at("MASS").at("2000005").at(1));
        mSQ[6] = to<double>(input.at("MASS").at("1000002").at(1));
        mSQ[7] = to<double>(input.at("MASS").at("1000004").at(1));
        mSQ[8] = to<double>(input.at("MASS").at("1000006").at(1));
        mSQ[9] = to<double>(input.at("MASS").at("2000002").at(1));
        mSQ[10] = to<double>(input.at("MASS").at("2000004").at(1));
        mSQ[11] = to<double>(input.at("MASS").at("2000006").at(1));
    } catch (...) {
        fprintf(stderr, "warning: No squark masses defined.\n");
    }
    try {
        // Neutralino masses
        mCH[0] = to<double>(input.at("MASS").at("1000022").at(1));
        mCH[1] = to<double>(input.at("MASS").at("1000023").at(1));
        mCH[2] = to<double>(input.at("MASS").at("1000025").at(1));
        mCH[3] = to<double>(input.at("MASS").at("1000035").at(1));
        // Chargino masses
        mCH[4] = to<double>(input.at("MASS").at("1000024").at(1));
        mCH[5] = to<double>(input.at("MASS").at("1000037").at(1));
    } catch (...) {
        fprintf(stderr, "warning: No gaugino masses defined.\n");
    }
    try {
        // Sneutrino masses
        mSL[0] = to<double>(input.at("MASS").at("1000012").at(1));
        mSL[1] = to<double>(input.at("MASS").at("1000014").at(1));
        mSL[2] = to<double>(input.at("MASS").at("1000016").at(1));
        // Charged Sleptons masses
        mSL[3] = to<double>(input.at("MASS").at("1000011").at(1));
        mSL[4] = to<double>(input.at("MASS").at("2000011").at(1));
        mSL[5] = to<double>(input.at("MASS").at("1000013").at(1));
        mSL[6] = to<double>(input.at("MASS").at("2000013").at(1));
        mSL[7] = to<double>(input.at("MASS").at("1000015").at(1));
        mSL[8] = to<double>(input.at("MASS").at("2000015").at(1));
    } catch (...) {
        fprintf(stderr, "warning: No slepton masses defined.\n");
    }

    // Mixings

    // Neutralino mixing
    try {
        for (int i0 = 0; i0 < 4; i0++) {
            for (int i1 = 0; i1 < 4; i1++) {
                RN[i0][i1] = to<double>(input.at("NMIX").at(i0 + 1, i1 + 1).at(2));
            }
        }

        for (int i0 = 0; i0 < 4; i0++) {
            if (mCH[i0] < 0.) {
                mCH[i0] *= -1.;
                for (int i1 = 0; i1 < 4; i1++) {
                    RN[i0][i1] *= complex<double>(0., 1.);
                }
            }
        }
    } catch (...) {
        fprintf(stderr, "warning: No neutralino mixing defined.\n");
    }

    // Chargino mixing matrix U
    try {
        for (int i0 = 0; i0 < 2; i0++) {
            for (int i1 = 0; i1 < 2; i1++) {
                RU[i0][i1] = to<double>(input.at("UMIX").at(i0 + 1, i1 + 1).at(2));
            }
        }

        // Chargino mixing matrix V
        for (int i0 = 0; i0 < 2; i0++) {
            for (int i1 = 0; i1 < 2; i1++) {
                RV[i0][i1] = to<double>(input.at("VMIX").at(i0 + 1, i1 + 1).at(2));
            }
        }
    } catch (...) {
        fprintf(stderr, "warning: No chargino mixing defined.\n");
    }


    // Sbottom mixing
    try {
        complex<double> RSB[2][2];
        for (int i0 = 0; i0 < 2; i0++) {
            for (int i1 = 0; i1 < 2; i1++) {
                RSB[i0][i1] = to<double>(input.at("SBOTMIX").at(i0 + 1, i1 + 1).at(2));
            }
        }

        for (int i0 = 0; i0 < 6; i0++) {
            for (int i1 = 0; i1 < 6; i1++) {
                if (i0 % 3 == 2 && i1 % 3 == 2) {
                    RSD[i0][i1] = RSB[i0 / 3][i1 / 3];
                } else {
                    RSD[i0][i1] = (i0 == i1 ? 1. : 0.);
                }
            }
        }
    } catch (...) {
        fprintf(stderr, "warning: No sbottom mixing defined.\n");
    }

    // Stop mixing
    try {
        complex<double> RST[2][2];
        for (int i0 = 0; i0 < 2; i0++) {
            for (int i1 = 0; i1 < 2; i1++) {
                RST[i0][i1] = to<double>(input.at("STOPMIX").at(i0 + 1, i1 + 1).at(2));
            }
        }

        // Heavy squark mixing matrix
        for (int i0 = 0; i0 < 6; i0++) {
            for (int i1 = 0; i1 < 6; i1++) {
                if (i0 % 3 == 2 && i1 % 3 == 2) {
                    RSU[i0][i1] = RST[i0 / 3][i1 / 3];
                } else {
                    RSU[i0][i1] = (i0 == i1 ? 1. : 0.);
                }
            }
        }
    } catch (...) {
        fprintf(stderr, "warning: No stop mixing defined.\n");
    }


    // Alpha
    try {
        double alpha = 0.0;
        for (Block::const_iterator line = input.at("ALPHA").begin();
                line != input.at("ALPHA").end(); ++line) {
            if (!line->is_data_line()) {
                continue;
            }

            alpha += to<double>(line->at(0));
        }

        RA[0][0] = cos(alpha);
        RA[0][1] = sin(alpha);
        RA[1][0] = -sin(alpha);
        RA[1][1] = cos(alpha);
    } catch (...) {
        fprintf(stderr, "warning: No block ALPHA defined.\n");
    }

    // Higgisino mixing
    try {
        double beta;
        beta = to<double>(input.at("HMIX").at("2").at(1));
        beta = atan(beta);
        RB[0][0] = sin(beta);
        RB[0][1] = cos(beta);
        RB[1][0] = -cos(beta);
        RB[1][1] = sin(beta);
    } catch (...) {
        fprintf(stderr, "warning: No Higgsino mixing defined.\n");
    }


    // Gauge
    double g[3];
    try {
        g[0] = to<double>(input.at("GAUGE").at("1").at(1));
        g[1] = to<double>(input.at("GAUGE").at("2").at(1));
        g[2] = to<double>(input.at("GAUGE").at("3").at(1));
        g1 = g[0];
        g2 = g[1];
        g3 = g[2];
    } catch (...) {
        double G_F = to<double>(input.at("SMINPUTS").at("2").at(1));
        double a_s = to<double>(input.at("SMINPUTS").at("3").at(1));
        xw = 1. - std::pow(mv[2] / mv[1], 2);
        sw = std::sqrt(xw);
        cw = std::sqrt(1. - xw);
        g2 = 2.0 * mv[2] * std::sqrt(M_SQRT2 * G_F);
        g3 = sqrt(4 * M_PI * a_s);
        g1 = sw / cw * g2; // not used
        fprintf(stderr, "warning: No gauge couplings defined. Using default values.\n");
    }

    // Stau mixing
    try {
        complex<double> STA[2][2];
        for (int i0 = 0; i0 < 2; i0++) {
            for (int i1 = 0; i1 < 2; i1++) {
                STA[i0][i1] = to<double>(input.at("STAUMIX").at(i0 + 1, i1 + 1).at(2));
            }
        }
        SSLM[0][0] = STA[0][0];
        SSLM[1][0] = STA[1][0];
        SSLM[0][1] = STA[0][1];
        SSLM[1][1] = STA[1][1];
    } catch (...) {
        fprintf(stderr, "warning: No stau mixing defined.\n");
    }


    // Yukawa Yd and Yu
    try {
        yq[0] =  to<double>(input.at("YD").at(1, 1).at(2));
    } catch (...) {
        fprintf(stderr, "warning: No Yukawa coupling Yd defined.\n");
    }
    try {
        yq[1] =  to<double>(input.at("YD").at(2, 2).at(2));
    } catch (...) {
        fprintf(stderr, "warning: No Yukawa coupling Ys defined.\n");
    }
    try {
        yq[2] =  to<double>(input.at("YD").at(3, 3).at(2));
    } catch (...) {
        fprintf(stderr, "warning: No Yukawa coupling Yb defined.\n");
    }
    try {
        yq[3] =  to<double>(input.at("YU").at(1, 1).at(2));
    } catch (...) {
        fprintf(stderr, "warning: No Yukawa coupling Yu defined.\n");
    }
    try {
        yq[4] =  to<double>(input.at("YU").at(2, 2).at(2));
    } catch (...) {
        fprintf(stderr, "warning: No Yukawa coupling Yc defined.\n");
    }
    try {
        yq[5] =  to<double>(input.at("YU").at(3, 3).at(2));
    } catch (...) {
        fprintf(stderr, "warning: No Yukawa coupling Yt defined.\n");
    }

}

void Parameters::write_log(const char *file)
{
    ofstream fout;
    fout.open(file);
    fout.precision(4);

    fout << "Parameters log\n";

    fout << "----------------------\n"
         << "DRb couplings\n"
         << "----------------------\n"
         << "g1 =\t" << g1    << "\n"
         << "g2 =\t" << g2    << "\n"
         << "g3 =\t" << g3    << "\n"
         << "xw =\t" << xw    << "\n"
         << "Yd =\t" << yq[0] << "\n"
         << "Ys =\t" << yq[1] << "\n"
         << "Yb =\t" << yq[2] << "\n"
         << "Yu =\t" << yq[3] << "\n"
         << "Yc =\t" << yq[4] << "\n"
         << "Yt =\t" << yq[5] << "\n";

    fout << "----------------------\n"
         << "Higgs sector\n"
         << "----------------------\n"
         << "mh =";
    for (int i0 = 0; i0 < 6; i0++) {
        fout << "\t" << mh[i0];
    }
    fout << "\nAlpha =";
    for (int i0 = 0; i0 < 2; i0++) {
        for (int i1 = 0; i1 < 2; i1++) {
            fout << "\t" << RA[i0][i1];
        }
        fout << "\n";
    }
    fout << "Beta  =";
    for (int i0 = 0; i0 < 2; i0++) {
        for (int i1 = 0; i1 < 2; i1++) {
            fout << "\t" << RB[i0][i1];
        }
        fout << "\n";
    }

    fout << "----------------------\n"
         << "Neutralino sector\n"
         << "----------------------\n"
         << "mN =";
    for (int i0 = 0; i0 < 4; i0++) {
        fout << "\t" << mCH[i0];
    }
    fout << "\nN  =";
    for (int i0 = 0; i0 < 4; i0++) {
        for (int i1 = 0; i1 < 4; i1++) {
            fout << "\t" << RN[i0][i1];
        }
        fout << "\n";
    }

    fout << "----------------------\n"
         << "Chargino sector\n"
         << "----------------------\n"
         << "mC =";
    for (int i0 = 0; i0 < 2; i0++) {
        fout << "\t" << mCH[i0 + 4];
    }
    fout << "\nU  =";
    for (int i0 = 0; i0 < 2; i0++) {
        for (int i1 = 0; i1 < 2; i1++) {
            fout << "\t" << RU[i0][i1];
        }
        fout << "\n";
    }
    fout << "V  =";
    for (int i0 = 0; i0 < 2; i0++) {
        for (int i1 = 0; i1 < 2; i1++) {
            fout << "\t" << RV[i0][i1];
        }
        fout << "\n";
    }

    fout << "----------------------\n"
         << "Gluino sector\n"
         << "----------------------\n"
         << "mG =\t" << mGL << "\n";

    fout << "----------------------\n"
         << "Down-squark sector\n"
         << "----------------------\n"
         << "mSD =";
    for (int i0 = 0; i0 < 6; i0++) {
        fout << "\t" << mSQ[i0];
    }
    fout << "\nRSD =";
    for (int i0 = 0; i0 < 6; i0++) {
        for (int i1 = 0; i1 < 6; i1++) {
            fout << "\t" << RSD[i0][i1];
        }
        fout << "\n";
    }

    fout << "----------------------\n"
         << "Up-squark sector\n"
         << "----------------------\n"
         << "mSU =";
    for (int i0 = 0; i0 < 6; i0++) {
        fout << "\t" << mSQ[i0 + 6];
    }
    fout << "\nRSU =";
    for (int i0 = 0; i0 < 6; i0++) {
        for (int i1 = 0; i1 < 6; i1++) {
            fout << "\t" << RSU[i0][i1];
        }
        fout << "\n";
    }

    fout << "----------------------\n"
         << "Sneutrino sector\n"
         << "----------------------\n"
         << "mSN =";
    for (int i0 = 0; i0 < 3; i0++) {
        fout << "\t" << mSL[i0];
    }
    fout << "\n";

    fout << "----------------------\n"
         << "Charged Slepton sector\n"
         << "----------------------\n"
         << "mCSL =";
    for (int i0 = 3; i0 < 9; i0++) {
        fout << "\t" << mSL[i0];
    }
    fout << "\nSSLM =";
    for (int i0 = 0; i0 < 2; i0++) {
        for (int i1 = 0; i1 < 2; i1++) {
            fout << "\t" << SSLM[i0][i1];
        }
        fout << "\n";
    }

    fout.close();
}
