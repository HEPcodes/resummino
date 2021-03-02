// This File will include all the Matrix Elements for gaugino squark pair
// production

#include "kinematics.h"
#include "npf.h"
#include "utils.h"
#include <complex>
#include <iostream>
using namespace std;

#define cA 3.0
#define cF (4.0 / 3.0)
#define NA 8.0
#define I2R 0.5
#define NC 3.0
#define MPIsI 1.0 / pow2(M_PI)

// TODO: LASSE

// Born
double FI::Mss_SQGA1() {
  return real(32.0 * pow2(ivs) * papb * pbp2 * (LL + RR));
}

double FI::Muu_SQGA1() {
  return real(-8.0 * ivu2s1 * ivu2s2 * pap2 *
              (m1s + m2s - 2.0 * (p1p2 - pap1 + pap2)) * (LL + RR));
}

double FI::Msu_SQGA1() {
  return real(-8.0 * ivs * ivu2s2 *
              (-2.0 * pow2(pap2) + (m2s - p1p2) * papb +
               pap2 * (pbp1 - 2.0 * pbp2) + pap1 * (2.0 * pap2 + pbp2)) *
              (LL + RR));
}

double FI::Mss_SQGA2() {
  return real(32.0 * pow2(ivs) * pap1 * papb * (LL + RR));
}

double FI::Muu_SQGA2() {
  return real(-8.0 * ivu1s1 * ivu1s2 * pbp1 *
              (m1s + m2s - 2.0 * (p1p2 + pbp1 - pbp2)) * (LL + RR));
}

double FI::Msu_SQGA2() {
  return real(-8.0 * ivs * ivu1s2 *
              (m1s * papb - p1p2 * papb - 2.0 * pap1 * pbp1 + pap2 * pbp1 -
               2.0 * pow2(pbp1) + pap1 * pbp2 + 2.0 * pbp1 * pbp2) *
              (LL + RR));
}

double FI::MG_SQGA1() { return 0.0; }

// squared matrix element (SS) w/ qqg vertex (vertex on the left) correction
// first vertex correction (quark-gluon-gluon loop)
// arguments: p1s = pa^2; p2s = pb^2; p3^2 = (pa + pb)^2
// ml1s = ml2s = ml3s = 0
double FI::Mss_qqg1_SQGA(double p1s, double p2s, double p3s, double ml1s,
                         double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 = -16.0 * papb * (c22 * pap2 * papb - 5.0 * c00 * pbp2 -
                                        2.0 * (-c1 + c12 + c22) * papb * pbp2 +
                                        c2 * papb * (pap2 + pbp2)) *
                        (LL + RR);

  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      16.0 * papb *
      (c22 * pap2 * papb - 2.0 * (4.0 * c00 + (-c1 + c12 + c22) * papb) * pbp2 +
       c2 * papb * (pap2 + pbp2)) *
      (LL + RR);

  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 = 48.0 * c00 * papb * pbp2 * (LL + RR);

  return real(pow2(cA) * cF / (32.0 * pow2(M_PI)) * pow2(ivs) *
              (me0 + me1 + me2));
}

// squared matrix element (SS) w/ qqg vertex (vertex on the left) correction
// second vertex correction (quark-gluon-gluon loop)
// --> Later (have to sum over squarks?)
// m1ls = msquark, ml2s=m3ls=mgluino
double FI::Mss_qqg2_SQGA(double p1s, double p2s, double p3s, double ml1s,
                         double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      -8.0 * papb *
      (LL * LsRs * (6.0 * c00 + c0 * ml2s + 2.0 * c12 * papb) * pbp2 +
       c0 * LsRs * m2 * ml2 * papb * RL +
       (c0 * LR * m2 * ml2 * papb +
        (6.0 * c00 + c0 * ml2s + 2.0 * c12 * papb) * pbp2 * RR) *
           RsLs -
       2.0 * c22 * papb * (pap2 - pbp2) * (LL * LsRs + RR * RsLs) -
       c2 * papb * (2.0 * LL * LsRs * (pap2 + pbp2) + LsRs * m2 * ml2 * RL +
                    (LR * m2 * ml2 + 2.0 * (pap2 + pbp2) * RR) * RsLs));
  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      8.0 * papb *
      (8.0 * c00 * pbp2 * (LL * LsRs + RR * RsLs) +
       2.0 * (c12 - c2 + c22) * papb * pbp2 * (LL * LsRs + RR * RsLs) +
       c0 * (LL * LsRs * ml2s * pbp2 + LsRs * m2 * ml2 * papb * RL +
             LR * m2 * ml2 * papb * RsLs + ml2s * pbp2 * RR * RsLs));

  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 = 16.0 * c00 * pbp2 * papb * RR * RsLs -
                        16.0 * c00 * pbp2 * papb * LL * LsRs;

  return real(pow2(cA) * cF / (32.0 * pow2(M_PI)) * pow2(ivs) *
              (me0 + me1 + me2));
}

// squared matrix element (SS) w/ qqg vertex (vertex on the left) correction
// third vertex correction (quark-gluon-gluon loop)
// arguments: p1s = pa^2; p2s = pb^2; p3^2 = (pa + pb)^2
// ml1s = ml2s = ml3s = 0
double FI::Mss_qqg3_SQGA(double p1s, double p2s, double p3s, double ml1s,
                         double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 = -32.0 * papb *
                        ((c2 + c22) * pap2 * papb - 3.0 * c00 * pbp2 +
                         (c0 + c1 - c12 + 2.0 * c2 - c22) * papb * pbp2) *
                        (LL + RR);

  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      32.0 * papb *
      ((c2 + c22) * pap2 * papb - 7.0 * c00 * pbp2 +
       (c0 + c1 - 2.0 * c12 + 3.0 * c2 - 2.0 * c22) * papb * pbp2) *
      (LL + RR);

  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 =
      32.0 * papb * (5.0 * c00 + (c12 - c2 + c22) * papb) * pbp2 * (LL + RR);

  return real(-cF / (32.0 * pow2(M_PI)) * pow2(ivs) * (me0 + me1 + me2));
}

// squared matrix element (SS) w/ qqg vertex (vertex on the left) correction
// fourth vertex correction (quark-gluon-gluon loop)
// arguments: p1s = pa^2; p2s = pb^2; p3^2 = (pa + pb)^2
// ml1s = mgluino, ml2s=ml3s=msquark
// Later sum over squarks
double FI::Mss_qqg4_SQGA(double p1s, double p2s, double p3s, double ml1s,
                         double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      4.0 * papb *
      (4.0 * c00 * LL * LsRs * pbp2 + c0 * LsRs * m1l * m2 * papb * RL +
       c0 * LR * m1l * m2 * papb * RsLs + 4.0 * c00 * pbp2 * RR * RsLs -
       4.0 * c22 * pap2 * papb * (LL * LsRs + RR * RsLs) +
       2.0 * c2 * papb * (-(LL * LsRs * pap2) + LsRs * m1l * m2 * RL +
                          LR * m1l * m2 * RsLs - pap2 * RR * RsLs));
  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = -16.0 * c00 * papb * pbp2 * (LL * LsRs + RR * RsLs);
  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 = 0;

  return real(cF / (32.0 * pow2(M_PI)) * pow2(ivs) * (me0 + me1 + me2));
}

// squared matrix element (SS) w/ qsqga vertex (vertex on the right) correction
// first vertex correction (quark-gluon-squark-loop)
// arguments: p1s = k^2; p2s = (-p1)^2; p3^2 = (-p2)^2
// ml1s = mgluino, ml2s=ml3s=msquark
// Later sum over squarks
double FI::Mss_qsqga1_SQGA(double p1s, double p2s, double p3s, double ml1s,
                           double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      -8.0 * papb * ((2.0 * c1 + c2) * (p1p2 * papb - pap2 * pbp1) +
                     (4.0 * c00 + 2.0 * c2 * m1s + c22 * m1s - 2.0 * c1 * pap1 -
                      2.0 * c12 * pap1 - 3.0 * c2 * pap1 - 2.0 * c22 * pap1 +
                      2.0 * c1 * papb + 2.0 * c11 * papb + 4.0 * c12 * papb +
                      2.0 * c2 * papb + 2.0 * c22 * papb -
                      2.0 * (2.0 * c1 + c12 + 2.0 * c2 + c22) * pbp1) *
                         pbp2) *
      (LL + RR);

  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      8 * papb * ((2.0 * c1 + c2) * (p1p2 * papb - pap2 * pbp1) +
                  (6.0 * c00 + 2.0 * c2 * m1s + c22 * m1s - 2.0 * c1 * pap1 -
                   2.0 * c12 * pap1 - 3.0 * c2 * pap1 - 2.0 * c22 * pap1 +
                   2.0 * c1 * papb + 2.0 * c11 * papb + 4.0 * c12 * papb +
                   2.0 * c2 * papb + 2.0 * c22 * papb -
                   2.0 * (2.0 * c1 + c12 + 2.0 * c2 + c22) * pbp1) *
                      pbp2) *
      (LL + RR);
  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 = -16.0 * c00 * papb * pbp2 * (LL + RR);

  return real(-cA * pow2(cF) / (16.0 * pow2(M_PI)) * pow2(ivs) *
              (me0 + me1 + me2));
}

// how to define couplings like L*Lsp or L*L or Rsp*Rsp??

// squared matrix element (SS) w/ qsqga vertex (vertex on the right) correction
// second vertex correction (quark-gluon-squark-loop)
// arguments: p1s = k^2; p2s = (-p1)^2; p3^2 = (-p2)^2
// ml1s = msquark, ml2s=mgluino, ml3s=msquark
// Later sum over squarks and quarks(?)
/*
double FI::Mss_qsqga2_SQGA(double p1s, double p2s, double p3s,
                         double ml1s, double ml2s, double ml3s, int ieps) {

  complex<double> c0,c00,c1,c2,c11,c12,c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s,  p3s,  ml1s,  ml2s,  ml3s, ieps, &c0, &c00, &c1, &c2, &c11,
&c12, &c22);

  complex<double> me0 = 8.0*papb*(-Lsp*Lsp*m2*ml2*
                                  ((c0 + c1 + c2)*papb - (c0 + c2)*pbp1)* RR +
                                2.0*LLsp*((c0 + c1 + c2)*(p1p2*papb - pap2*pbp1)
+
                                         (-4.0*c00 - c2*m1s - c22*m1s + c0*pap1
+
                                                     c1*pap1 + 2.0*c12*pap1 +
3.0*c2*pap1 +
                                                     2.0*c22*pap1 - 2.0*c0*papb
- 4.0*c1*papb -
                                                     2.0*c11*papb - 4.0*c12*papb
- 4.0*c2*papb -
                                          2.0*c22*papb + 2.0*(c12 + c2 +
c22)*pbp1)* pbp2)*RRsp -
                                  LL*m2*ml2*((c0 + c1 + c2)*papb - (c0 +
c2)*pbp1)*RspRsp;

    // I1
  SetC(p1s, p2s,  p3s,  ml1s,  ml2s,  ml3s, ieps+1, &c0, &c00, &c1, &c2, &c11,
&c12, &c22);

  complex<double> me1 = 8*papb*(Power(Lsp,2)*m2*ml2*((c0 + c1 + c2)*papb - (c0 +
c2)*pbp1)*Power(R,2) - 2*L*Lsp*
                                ((c0 + c1 + c2)*(p1p2*papb - pap2*pbp1) +
(-6*c00 - c2*m1s - c22*m1s + c0*pap1 +
                                                                             c1*pap1
+ 2*c12*pap1 + 3*c2*pap1 +
                                                                             2*c22*pap1
- 2*c0*papb - 4*c1*papb -
                                                                             2*c11*papb
- 4*c12*papb - 4*c2*papb -
                                                                  2*c22*papb +
2*(c12 + c2 + c22)*pbp1)*
                                                                 pbp2)*R*Rsp +
                                                                Power(L,2)*m2*ml2*
                                                                ((c0 + c1 +
c2)*papb - (c0 + c2)*pbp1)*
                                Power(Rsp,2));
   // I2
  SetC(p1s, p2s,  p3s,  ml1s,  ml2s,  ml3s, ieps+2, &c0, &c00, &c1, &c2, &c11,
&c12, &c22);

  complex<double> me2 = - 32.0*c00*pbp2*papb*L*R*Rsp*Lsp;

  return real( cA*pow2(cF)  / (16.0 * pow2(M_PI) ) *  pow2(ivs) * ( me0 + me1 +
me2));

}*/

// UU-channel
// squared matrix element (uu) w/ sqaurk-sqaurk-gaugino vertex (top)
// first correction

double FI::Muu_qsqga1_SQGA(double p1s, double p2s, double p3s, double ml1s,
                           double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      2.0 * (-((2.0 * c1 + c2) * m2s * pap1) +
             (4.0 * c00 + c11 * m1s + 2.0 * c12 * m1s + c2 * m1s + c22 * m1s +
              2.0 * c2 * m2s + c22 * m2s + 2.0 * c12 * p1p2 + 4.0 * c2 * p1p2 +
              2.0 * c22 * p1p2 + c1 * (m1s + 4.0 * p1p2) -
              2.0 * (c12 + c2 + c22) * pap1) *
                 pap2 -
             2.0 * (2.0 * c2 + c22) * pow2(pap2)) *
      (m1s + m2s - 2.0 * (p1p2 - pap1 + pap2)) * (LL + RR);

  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      -4.0 * c00 * pap2 * (m1s + m2s - 2.0 * (p1p2 - pap1 + pap2)) * (LL + RR);

  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 = 0;

  return real(-cA * pow2(cF) / (16.0 * pow2(M_PI)) * ivu1s1 * ivu1s2 *
              (me0 + me1 + me2));
}

// squared matrix element (uu) w/ sqaurk-sqaurk-gaugino vertex (top)
// second correction (squark-gluino-quark-loop)
// yet to come -> Couplings

// squared matrix element (uu) w/ gluon-squark-squark vertex (bottom) correction
// first correction
double FI::Muu_gsqsq1_SQGA(double p1s, double p2s, double p3s, double ml1s,
                           double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      -2.0 * pap2 *
      (-(c1 * p1p2 * papb) + 2.0 * c11 * p1p2 * papb + 3.0 * c2 * p1p2 * papb -
       2.0 * c22 * p1p2 * papb + c1 * pap1 * papb - 2.0 * c11 * pap1 * papb -
       3.0 * c2 * pap1 * papb + 2.0 * c22 * pap1 * papb +
       2.0 * c12 * pow(papb, 2) + 2.0 * c2 * pow(papb, 2) -
       2.0 * c22 * pow(papb, 2) - c1 * m2s * pbp1 + 2.0 * c11 * m2s * pbp1 +
       3.0 * c2 * m2s * pbp1 - 2.0 * c22 * m2s * pbp1 + 2.0 * c1 * pap2 * pbp1 -
       4.0 * c11 * pap2 * pbp1 - 6.0 * c2 * pap2 * pbp1 +
       4.0 * c22 * pap2 * pbp1 + 2.0 * c12 * papb * pbp1 +
       2.0 * c2 * papb * pbp1 - 2.0 * c22 * papb * pbp1 +
       6.0 * c00 * (2.0 * m2s - 2.0 * p1p2 + 2.0 * pap1 - 4.0 * pap2 + papb +
                    pbp1 - pbp2) +
       (c1 * p1p2 - 2.0 * c11 * p1p2 - 3.0 * c2 * p1p2 + 2.0 * c22 * p1p2 -
        c1 * pap1 + 2.0 * c11 * pap1 + 3.0 * c2 * pap1 - 2.0 * c22 * pap1 -
        4.0 * c12 * papb - 4.0 * c2 * papb + 4.0 * c22 * papb -
        2.0 * (c12 + c2 - c22) * pbp1) *
           pbp2 +
       2.0 * (c12 + c2 - c22) * pow(pbp2, 2) +
       c0 * (pap1 * papb - m2s * pbp1 + 2.0 * pap2 * pbp1 - pap1 * pbp2 +
             p1p2 * (-papb + pbp2))) *
      (LL + RR);

  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      8.0 * c00 * pap2 *
      (2.0 * m2s - 2.0 * p1p2 + 2.0 * pap1 - 4.0 * pap2 + papb + pbp1 - pbp2) *
      (LL + RR);

  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 = 0;

  return real(cF * pow2(cA) / (32.0 * pow2(M_PI)) * ivu1s1 * ivu1s2 *
              (me0 + me1 + me2));
}

// squared matrix element (uu) w/ gluon-squark-squark vertex (bottom) correction
// second correction
double FI::Muu_gsqsq2_SQGA(double p1s, double p2s, double p3s, double ml1s,
                           double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22, c001, c002, c003, c111, c112,
      c122, c222;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22, &c001, &c002, &c111, &c112, &c122, &c222);
  complex<double> c121 = c112;

  complex<double> me0 =
      -2.0 * pap2 *
      (6.0 * c002 * m2s + c1 * pow(m2s, 2) + 2.0 * c11 * pow(m2s, 2) +
       c111 * pow(m2s, 2) + 2.0 * c112 * pow(m2s, 2) + 4.0 * c12 * pow(m2s, 2) +
       c121 * pow(m2s, 2) + 3.0 * c122 * pow(m2s, 2) + c2 * pow(m2s, 2) +
       2.0 * c22 * pow(m2s, 2) + c222 * pow(m2s, 2) + c1 * m2s * pow(ml2, 2) +
       c2 * m2s * pow(ml2, 2) - 6.0 * c002 * p1p2 - c1 * m2s * p1p2 -
       2.0 * c11 * m2s * p1p2 - c111 * m2s * p1p2 - 2.0 * c112 * m2s * p1p2 -
       4.0 * c12 * m2s * p1p2 - c121 * m2s * p1p2 - 3.0 * c122 * m2s * p1p2 -
       c2 * m2s * p1p2 - 2.0 * c22 * m2s * p1p2 - c222 * m2s * p1p2 -
       c1 * pow(ml2, 2) * p1p2 - c2 * pow(ml2, 2) * p1p2 + 6.0 * c002 * pap1 +
       c1 * m2s * pap1 + 2.0 * c11 * m2s * pap1 + c111 * m2s * pap1 +
       2.0 * c112 * m2s * pap1 + 4.0 * c12 * m2s * pap1 + c121 * m2s * pap1 +
       3.0 * c122 * m2s * pap1 + c2 * m2s * pap1 + 2.0 * c22 * m2s * pap1 +
       c222 * m2s * pap1 + c1 * pow(ml2, 2) * pap1 + c2 * pow(ml2, 2) * pap1 +
       6.0 * c001 * (m2s - p1p2 + pap1 - 2.0 * pap2) - 12.0 * c002 * pap2 -
       4.0 * c1 * m2s * pap2 - 8.0 * c11 * m2s * pap2 -
       4.0 * c111 * m2s * pap2 - 8.0 * c112 * m2s * pap2 -
       16.0 * c12 * m2s * pap2 - 4.0 * c121 * m2s * pap2 -
       12.0 * c122 * m2s * pap2 - 4.0 * c2 * m2s * pap2 -
       8.0 * c22 * m2s * pap2 - 4.0 * c222 * m2s * pap2 -
       2.0 * c1 * pow(ml2, 2) * pap2 - 2.0 * c2 * pow(ml2, 2) * pap2 +
       2.0 * c1 * p1p2 * pap2 + 4.0 * c11 * p1p2 * pap2 +
       2.0 * c111 * p1p2 * pap2 + 4.0 * c112 * p1p2 * pap2 +
       8.0 * c12 * p1p2 * pap2 + 2.0 * c121 * p1p2 * pap2 +
       6.0 * c122 * p1p2 * pap2 + 2.0 * c2 * p1p2 * pap2 +
       4.0 * c22 * p1p2 * pap2 + 2.0 * c222 * p1p2 * pap2 -
       2.0 * c1 * pap1 * pap2 - 4.0 * c11 * pap1 * pap2 -
       2.0 * c111 * pap1 * pap2 - 4.0 * c112 * pap1 * pap2 -
       8.0 * c12 * pap1 * pap2 - 2.0 * c121 * pap1 * pap2 -
       6.0 * c122 * pap1 * pap2 - 2.0 * c2 * pap1 * pap2 -
       4.0 * c22 * pap1 * pap2 - 2.0 * c222 * pap1 * pap2 +
       4.0 * c1 * pow(pap2, 2) + 8.0 * c11 * pow(pap2, 2) +
       4.0 * c111 * pow(pap2, 2) + 8.0 * c112 * pow(pap2, 2) +
       16.0 * c12 * pow(pap2, 2) + 4.0 * c121 * pow(pap2, 2) +
       12.0 * c122 * pow(pap2, 2) + 4.0 * c2 * pow(pap2, 2) +
       8.0 * c22 * pow(pap2, 2) + 4.0 * c222 * pow(pap2, 2) +
       6.0 * c002 * papb + c1 * m2s * papb + c11 * m2s * papb +
       2.0 * c112 * m2s * papb + 6.0 * c12 * m2s * papb + c121 * m2s * papb +
       6.0 * c122 * m2s * papb + 2.0 * c2 * m2s * papb +
       5.0 * c22 * m2s * papb + 3.0 * c222 * m2s * papb +
       c2 * pow(ml2, 2) * papb - c112 * p1p2 * papb - 4.0 * c12 * p1p2 * papb -
       c121 * p1p2 * papb - 4.0 * c122 * p1p2 * papb - 2.0 * c2 * p1p2 * papb -
       4.0 * c22 * p1p2 * papb - 2.0 * c222 * p1p2 * papb + c112 * pap1 * papb +
       4.0 * c12 * pap1 * papb + c121 * pap1 * papb + 4.0 * c122 * pap1 * papb +
       2.0 * c2 * pap1 * papb + 4.0 * c22 * pap1 * papb +
       2.0 * c222 * pap1 * papb - 2.0 * c1 * pap2 * papb -
       2.0 * c11 * pap2 * papb - 4.0 * c112 * pap2 * papb -
       12.0 * c12 * pap2 * papb - 2.0 * c121 * pap2 * papb -
       12.0 * c122 * pap2 * papb - 4.0 * c2 * pap2 * papb -
       10.0 * c22 * pap2 * papb - 6.0 * c222 * pap2 * papb +
       2.0 * c12 * pow(papb, 2) + 2.0 * c122 * pow(papb, 2) +
       2.0 * c22 * pow(papb, 2) + 2.0 * c222 * pow(papb, 2) +
       6.0 * c002 * pbp1 + c1 * m2s * pbp1 + c11 * m2s * pbp1 +
       c112 * m2s * pbp1 + 2.0 * c12 * m2s * pbp1 + 2.0 * c122 * m2s * pbp1 +
       c22 * m2s * pbp1 + c222 * m2s * pbp1 + c2 * pow(ml2, 2) * pbp1 -
       2.0 * c1 * pap2 * pbp1 - 2.0 * c11 * pap2 * pbp1 -
       2.0 * c112 * pap2 * pbp1 - 4.0 * c12 * pap2 * pbp1 -
       4.0 * c122 * pap2 * pbp1 - 2.0 * c22 * pap2 * pbp1 -
       2.0 * c222 * pap2 * pbp1 + 2.0 * c12 * papb * pbp1 +
       2.0 * c122 * papb * pbp1 + 2.0 * c22 * papb * pbp1 +
       2.0 * c222 * papb * pbp1 +
       4.0 * c00 * (2.0 * m2s - 2.0 * p1p2 + 2.0 * pap1 - 4.0 * pap2 + papb +
                    pbp1 - pbp2) -
       (6.0 * c002 + c1 * m2s + c11 * m2s + 2.0 * c112 * m2s + 6.0 * c12 * m2s +
        c121 * m2s + 6.0 * c122 * m2s + 2.0 * c2 * m2s + 5.0 * c22 * m2s +
        3.0 * c222 * m2s + c2 * pow(ml2, 2) - c112 * p1p2 - 4.0 * c12 * p1p2 -
        c121 * p1p2 - 4.0 * c122 * p1p2 - 2.0 * c2 * p1p2 - 4.0 * c22 * p1p2 -
        2.0 * c222 * p1p2 + c112 * pap1 + 4.0 * c12 * pap1 + c121 * pap1 +
        4.0 * c122 * pap1 + 2.0 * c2 * pap1 + 4.0 * c22 * pap1 +
        2.0 * c222 * pap1 - 2.0 * c1 * pap2 - 2.0 * c11 * pap2 -
        4.0 * c112 * pap2 - 12.0 * c12 * pap2 - 2.0 * c121 * pap2 -
        12.0 * c122 * pap2 - 4.0 * c2 * pap2 - 10.0 * c22 * pap2 -
        6.0 * c222 * pap2 + 4.0 * c12 * papb + 4.0 * c122 * papb +
        4.0 * c22 * papb + 4.0 * c222 * papb +
        2.0 * (c12 + c122 + c22 + c222) * pbp1) *
           pbp2 +
       2.0 * (c12 + c122 + c22 + c222) * pow(pbp2, 2)) *
      (LL + RR) * (LsRs + RsLs);

  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22, &c001, &c002, &c111, &c112, &c122, &c222);
  complex<double> me1 =
      4.0 * pap2 *
      (c001 * (m2s - p1p2 + pap1 - 2.0 * pap2) +
       c00 * (2.0 * m2s - 2.0 * p1p2 + 2.0 * pap1 - 4.0 * pap2 + papb + pbp1 -
              pbp2) +
       c002 * (m2s - p1p2 + pap1 - 2.0 * pap2 + papb + pbp1 - pbp2)) *
      (LL + RR) * (LsRs + RsLs);

  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22, &c001, &c002, &c111, &c112, &c122, &c222);

  complex<double> me2 = 0;

  return real(cF * pow2(cA) / (32.0 * pow2(M_PI)) * ivu1s1 * ivu1s2 *
              (me0 + me1 + me2));
}
