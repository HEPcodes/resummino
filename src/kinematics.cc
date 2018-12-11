// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

#include <stdio.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <cstdlib>

#include "kinematics.h"
#include "utils.h"

using namespace std;

// Set Born kinematics
// uses mandelstam variables s and t to set the scalar products.
void FI::SetKinematic(const double mi, const double mj,
                      const double s,  const double t) {
    m1 = mi; // Mass of outgoing particle
    m2 = mj; // Mass of outgoing antiparticle
    m1s = pow2(m1);
    m2s = pow2(m2);
    papb = 0.5 * s;
    p1p2 = 0.5 * (s - m1s - m2s);
    const double u = m1s + m2s - s - t;
    pap1 = 0.5 * (m1s - t);
    pbp2 = 0.5 * (m2s - t);
    pbp1 = 0.5 * (m1s - u);
    pap2 = 0.5 * (m2s - u);
}

// Set 2->3 kinematics using the helicity frame 
// (see particle kinematics - Byckling and Kajantie for further information).
void FI::SetKinematic_HelicityFrame(const double mi,  const double mj,  const double s,
                       const double s2, const double t1,
                       const double s1, const double phi, int flag) {

    m1 = mi;
    m1s = pow2(m1);
    m2 = mj;
    m2s = pow2(m2);
    double cosphi = cos(phi);
    double t2;
    double detM;

    // Helicity frame for dps3_23
    if (flag == 1) {
        // express t in terms of the helicity angle phi (Eq. 8.9 on p. 136)
        detM = -((s - s2)*(m2s*s - s1*s2)) + (s*(-2*m2s - s + s1) 
                                              + (s + s1)*s2)*t1 + m1s*(m2s*s - 2*s*s2 + s1*s2 + s*t1 - s1*t1);
        t2 = 1.0/(Rkln(s,s2,m1s)) * (-detM + 2.0 * sqrt(G(s,t1,s2,0,0,m1s) * G(s1,s2,s,m2s,m1s,0)) * cosphi);

        papb = 0.5 * s;
        p1p2 = 0.5 * (s1 - m1s - m2s);
        pap1 = 0.5 * (m1s -t1);
        pap3 = 0.5 * (s - s1 + t2); 
        pbp1 = 0.5 * (s - s2 + t1); 
        pbp3 = 0.5 * (-t2); 
        pap2 = 0.5 * (s1 + t1 - t2 - m1s);
        pbp2 = 0.5 * (s2 + t2 - t1); 
        p1p3 = 0.5 * ( s -s1 - s2 + m2s); 
        p2p3 = 0.5 * (s2 - m2s ); 

    }

    // Helicity frame for dps3_13 (p1 <-> p2)
    if (flag == 2) {
        detM = -((s - s2)*(m1s*s - s1*s2)) + (s*(-2*m1s - s + s1) 
                                              + (s + s1)*s2)*t1 + m2s*(m1s*s - 2*s*s2 + s1*s2 + s*t1 - s1*t1);
        t2 = 1.0/(Rkln(s,s2,m2s)) * (-detM + 2.0 * sqrt(G(s,t1,s2,0,0,m2s) * G(s1,s2,s,m1s,m2s,0)) * cosphi);

        papb = 0.5 * s;
        p1p2 = 0.5 * (s1 - m1s - m2s);
        pap2 = 0.5 * (m2s -t1);
        pap3 = 0.5 * (s - s1 + t2);
        pbp2 = 0.5 * (s - s2 + t1);
        pbp3 = 0.5 * (-t2);
        pap1 = 0.5 * (s1 + t1 - t2 - m2s);
        pbp1 = 0.5 * (s2 + t2 - t1);
        p2p3 = 0.5 * ( s -s1 - s2 + m1s);
        p1p3 = 0.5 * (s2 - m1s );
    }
}

// Sets on-shell kinematics s23.
void FI::SetKinematicRES23(const double mSQs, double &factor) {

    // Jacobian before reshuffling.
    double Jac =  kln(2.0*papb, 2.0 * p2p3 + m2s,m1s) 
        * kln(2.0 * p2p3 + m2s,m2s,0) / (2.0 * p2p3 + m2s); 

    /* The reshuffling */

    // (this can maybe be improved using dipole kinematics
    //  for final emitter and spectator ( for more information see arxiv:1305.4061v2) )

    // temporary variables for on shell kinematics
    double Qs = m1s + m2s + 2.* ( p1p3 + p2p3 + p1p2);
    double p1Q = m1s + p1p2 + p1p3;
    double paQ = pap1 + pap2 + pap3;
    double pbQ = pbp1 + pbp2 + pbp3;
    
    pap1 = kln(Qs,mSQs,m1s) / kln(Qs,2.0*p2p3 + m2s,m1s) * (pap1 - p1Q/Qs * paQ) 
        + (Qs + m1s - mSQs) / (2*Qs) * paQ;

    // set pbp1 to on shell kinematics
    pbp1 = kln(Qs,mSQs,m1s) / kln(Qs,2.0*p2p3 + m2s,m1s) * (pbp1 - p1Q/Qs * pbQ) 
        + (Qs + m1s - mSQs) / (2*Qs) * pbQ;

    // set p2p3 to on shell kinematics
    p2p3 = 0.5 * (mSQs - m2s);

    // new mandelstams
    double t1 = - 2.0 * pap1 + m1s;
    double s2 = 2.0 * p2p3 + m2s;
    double t2 = -2.0 * pbp3;
    double s1 = 2.0 * p1p2 + m1s + m2s;
    double s = 2.0 * papb;

    // new momenta
    pap2 = 0.5 * (s1 + t1 - t2 - m1s);
    pbp2 = 0.5 * (s2 + t2 - t1); 
    p1p3 = 0.5 * ( s -s1 - s2 + m2s); 
        
    // Jacobian after reshuffling.
    double Jac_tilde = kln(2.0*papb, 2.0 * p2p3 + m2s,m1s) 
        * kln(2.0 * p2p3 + m2s,m2s,0) / (2.0 * p2p3 + m2s);
  
    // set factor to rescale the full dP3s to the restricted phase space.
     factor = Jac_tilde/Jac;    
}

// Sets on-shell kinematics for s13 (similar to RES23)
void FI::SetKinematicRES13(const double mSQs, double &factor) {

    double Jac = kln(2.0*papb, 2.0 * p1p3 + m1s,m2s) 
        * kln(2.0 * p1p3 + m1s,m1s,0) / (2.0 * p1p3 + m1s);

    double paQ = pap1 + pap2 + pap3;
    double pbQ = pbp1 + pbp2 + pbp3;
    double Qs = m1s + m2s + 2.* ( p1p3 + p2p3 + p1p2);
    double p2Q = m2s + p1p2 + p2p3;

    pap2 = kln(Qs,mSQs,m2s) / kln(Qs,2.0*p1p3 + m1s,m2s) * (pap2 - p2Q/Qs * paQ) 
        + (Qs + m2s - mSQs) / (2*Qs) * paQ;

    pbp2 = kln(Qs,mSQs,m2s) / kln(Qs,2.0*p1p3 + m1s,m2s) * (pbp2 - p2Q/Qs * pbQ) 
        + (Qs + m2s - mSQs) / (2*Qs) * pbQ;
        
    p1p3 = 0.5 * (mSQs - m1s); 

        // new mandelstams
    double t1 = - 2.0 * pap2 + m2s;
    double s2 = 2.0 * p1p3 + m1s;
    double t2 = -2.0 * pbp3;
    double s1 = 2.0 * p1p2 + m1s + m2s;
    double s = 2.0 * papb;
    
    // new momenta
    pap1 = 0.5 * (s1 + t1 - t2 - m2s);
    pbp1 = 0.5 * (s2 + t2 - t1);
    p2p3 = 0.5 * ( s -s1 - s2 + m1s);
  
    double Jac_tilde = kln(2.0*papb, 2.0 * p1p3 + m1s,m2s) 
        * kln(2.0 * p1p3 + m1s,m1s,0) / (2.0 * p1p3 + m1s); 

    factor = Jac_tilde/Jac;
}

// set 2->3 kinematics 
void FI::SetKinematic(const double mi,  const double mj,  const double s,
                      const double invariant_mass2, const double pt2,
                      const double th, const double ph, const int ys) {

  // set the five linear independent scalar products (papb,p1p2,pap3,pbp1,pap1) 
  // in terms of the new invariants used for the integration (inv_mass2, pt2, th, ph,s).
    m1 = mi;
    m1s = pow2(m1);
    m2 = mj;
    m2s = pow2(m2);
    papb = 0.5 * s; //1
    p1p2 = 0.5 * (invariant_mass2 - m1s - m2s); //2

    if (ys == 0) {
      pap3 = 0.25 * (s - invariant_mass2 + sqrt(pow2(s - invariant_mass2) - 4.0 * s * pt2)); //3
//      pbp3 = 0.25 * (s - invariant_mass2 - sqrt(pow2(s - invariant_mass2) - 4.0 * s * pt2));
    } else if (ys == 1) {
      pap3 = 0.25 * (s - invariant_mass2 - sqrt(pow(s - invariant_mass2, 2) - 4.0 * s * pt2)); //3
//      pbp3 = 0.25 * (s - invariant_mass2 + sqrt(pow(s - invariant_mass2, 2) - 4.0 * s * pt2));
    } else {
        cout << "Rapidity flag = " << ys << "\n";
        exit(0);
    }

    double kM12, cA, sA;
    kM12 = kln(invariant_mass2, m1s, m2s);
    pap1 = .25 * (invariant_mass2 + 2.0 * (-0.5 * (m1s + m2s) - p1p2 - pap3 + papb)) / invariant_mass2 
      * (invariant_mass2 + m1s - m2s - kM12 * cos(th)); //4

    cA = 1.0 - 2.0 * s * invariant_mass2 / (invariant_mass2 + 2.0 * pap3) 
      / (invariant_mass2 + 2.0 * (-0.5 * (m1s + m2s) - p1p2 - pap3 + papb));
    sA = sqrt(1.0 - cA * cA);
    
    pbp1 = .25 * (invariant_mass2 + 2.0 * pap3) / invariant_mass2 *
      (invariant_mass2 + m1s - m2s - kM12 * (cA * cos(th) + sA * sin(th) * cos(ph))); //5

    // set the remaining 5 dependent scalar products
    pbp3 = -0.5 * (m1s + m2s) - p1p2 - pap3 + papb;
    pap2 = -pap1 - pap3 + papb;
    pbp2 = 0.5 * (m1s + m2s) + p1p2 + pap3 - pbp1;
    p1p3 = -m1s - p1p2 + pap1 + pbp1;
    p2p3 = 0.5*(m1s - m2s) - pap1 + papb - pbp1;

}

// Dipole kinematics for Emitter & Spectator
// (see Catani Seymour Dipole papers: hep-ph/0201036 and hep-ph/9605323 )

//Emitter A & Spectator B
void FI::SetDipKinematicAB(double &x) {
  x = 1.0 - (pap3 + pbp3) / papb;


  const double papb_temp = x * papb;    

  const double pap1_temp = (x*(p1p3*pap3*(pap3 - 2*papb + pbp3) - 
                               pbp1*(2*pow(pap3,2) + pap3*(2*pbp3 + papb*(-5 + x)) - 
                                     2*papb*(pbp3 + papb*(-1 + x))) + 
                               p1p3*(pap3 - 2*papb)*papb*x + 
                               pap1*(-(pap3*(papb - 2*pbp3 + papb*x)) + 
                                     2*(pow(papb - pbp3,2) + pow(papb,2)*x))))/
    ((pap3 - papb + pbp3)*(pap3 - 2*papb + 2*pbp3 + pap3*x - 
                           2*papb*x));

  const double pap2_temp = (x*(p2p3*pap3*(pap3 - 2*papb + pbp3) - 
                               pbp2*(2*pow(pap3,2) + pap3*(2*pbp3 + papb*(-5 + x)) - 
                                     2*papb*(pbp3 + papb*(-1 + x))) + 
                               p2p3*(pap3 - 2*papb)*papb*x + 
                               pap2*(-(pap3*(papb - 2*pbp3 + papb*x)) + 
                                     2*(pow(papb - pbp3,2) + pow(papb,2)*x))))/
    ((pap3 - papb + pbp3)*(pap3 - 2*papb + 2*pbp3 + pap3*x - 
                           2*papb*x));

  const double pbp1_temp = ((pap3*pbp1 - (p1p3 - pap1)*(papb - pbp3))*(pap3 - papb + pbp3) + 
                            (pap1*(papb - pbp3)*(pap3 + pbp3) + p1p3*papb*(-papb + pbp3) + 
                             pbp1*(pow(pap3,2) + 2*papb*(papb - pbp3) + 
                                   pap3*(-2*papb + pbp3)))*x + 
                            papb*(p1p3*(pap3 - 2*papb) - pap3*pbp1 + 2*papb*pbp1 + 
                                  pap1*(papb + pbp3))*pow(x,2))/
    ((pap3 - papb + pbp3)*(pap3 - 2*papb + 2*pbp3 + pap3*x - 
                           2*papb*x));


  const double pbp2_temp = ((pap3*pbp2 - (p2p3 - pap2)*(papb - pbp3))*(pap3 - papb + pbp3) + 
                            (pap2*(papb - pbp3)*(pap3 + pbp3) + p2p3*papb*(-papb + pbp3) + 
                             pbp2*(pow(pap3,2) + 2*papb*(papb - pbp3) + 
                                   pap3*(-2*papb + pbp3)))*x + 
                            papb*(p2p3*(pap3 - 2*papb) - pap3*pbp2 + 2*papb*pbp2 + 
                                  pap2*(papb + pbp3))*pow(x,2))/
    ((pap3 - papb + pbp3)*(pap3 - 2*papb + 2*pbp3 + pap3*x - 
                           2*papb*x));

  papb = papb_temp;
  pap1 = pap1_temp;
  pap2 = pap2_temp;
  pbp1 = pbp1_temp;
  pbp2 = pbp2_temp;

}

//Emitter B & Spectator A
void FI::SetDipKinematicBA(double &x) {
  x = 1.0 - (pap3 + pbp3) / papb;
 
  const double papb_temp = x * papb;    

  const double pap1_temp =((pap3 - papb + pbp3)*((pap3 - papb)*(p1p3 - pbp1) + pap1*pbp3) 
                           + ((pap3 - papb)*(p1p3*papb - 2*pap1*papb - pap3*pbp1) 
                              + (pap1*(pap3 - 2*papb) + (-pap3 + papb)*pbp1)*pbp3 + pap1*pow(pbp3,2))*x 
                           + papb*((pap3 + papb)*pbp1 + pap1*(2*papb - pbp3) + p1p3*(-2*papb + pbp3))*pow(x,2))/
    ((pap3 - papb + pbp3)*(2*pap3 - 2*papb + pbp3 - 2*papb*x + pbp3*x));

  const double pap2_temp = ((pap3 - papb + pbp3)*((pap3 - papb)*(p2p3 - pbp2) + pap2*pbp3) 
                            + ((pap3 - papb)*(p2p3*papb - 2*pap2*papb - pap3*pbp2) 
                               + (pap2*(pap3 - 2*papb) + (-pap3 + papb)*pbp2)*pbp3 + pap2*pow(pbp3,2))*x + 
                            papb*((pap3 + papb)*pbp2 + pap2*(2*papb - pbp3) + p2p3*(-2*papb + pbp3))*pow(x,2))/
    ((pap3 - papb + pbp3)*(2*pap3 - 2*papb + pbp3 - 2*papb*x + pbp3*x));

  const double pbp1_temp = (x*(2*pow(pap3,2)*pbp1 + pap3*(-4*papb*pbp1 + p1p3*pbp3 + 2*pbp1*pbp3) 
                               + pap1*(2*pap3*(papb - pbp3) + (2*papb - pbp3)*(2*pbp3 + papb*(-1 + x))) + 
                               (2*papb - pbp3)*(papb*pbp1 - p1p3*pbp3 + papb*(-p1p3 + pbp1)*x)))/
    ((pap3 - papb + pbp3)*(2*pap3 - 2*papb + pbp3 - 2*papb*x + pbp3*x));


  const double pbp2_temp = (x*(2*pow(pap3,2)*pbp2 + pap3*(-4*papb*pbp2 + p2p3*pbp3 + 2*pbp2*pbp3) 
                               + pap2*(2*pap3*(papb - pbp3) + (2*papb - pbp3)*(2*pbp3 + papb*(-1 + x))) + 
                               (2*papb - pbp3)*(papb*pbp2 - p2p3*pbp3 + papb*(-p2p3 + pbp2)*x)))/
    ((pap3 - papb + pbp3)*(2*pap3 - 2*papb + pbp3 - 2*papb*x + pbp3*x));

  papb = papb_temp;
  pap1 = pap1_temp;
  pap2 = pap2_temp;
  pbp1 = pbp1_temp;
  pbp2 = pbp2_temp;

}


// hep-ph/0201036 p.29 eq 5.73-5.80
//Emitter A & Spectator 1
void FI::SetDipKinematicA1(double &x, double &zi, double &zj) {

  x = 1.0 - p1p3/(pap1 + pap3);
  zi = pap3/(pap3 + pap1);
  zj = pap1/(pap3 + pap1);

  
  const double papb_temp = papb*x;
  const double pap1_temp = (pap1 + pap3)*x;
  const double pap2_temp = pap2*x;                         
  const double pbp1_temp = pbp1 + pbp3 + papb*(-1.0 + x);
  const double pbp2_temp = pbp2;

  const double pa_old_p1_new = pap3 + pap1;

  papb = papb_temp;
  pap1 = pap1_temp;
  pap2 = pap2_temp;
  pbp1 = pbp1_temp;
  pbp2 = pbp2_temp;

}

//Emitter B & Spectator 1
void FI::SetDipKinematicB1(double &x, double &zi, double &zj) {

  x = 1.0 - p1p3/(pbp1 + pbp3);
  zi = pbp3/(pbp3 + pbp1);
  zj = pbp1/(pbp3 + pbp1);


  const double papb_temp = papb*x;
  const double pap1_temp = pap1 + pap3 + papb*(-1.0 + x);
  const double pap2_temp = pap2;                         
  const double pbp1_temp = (pbp1 + pbp3)*x;
  const double pbp2_temp = pbp2*x;

  papb = papb_temp;
  pap1 = pap1_temp;
  pap2 = pap2_temp;
  pbp1 = pbp1_temp;
  pbp2 = pbp2_temp;

}

// hep-ph/0201036 p.25 eq 5.42-5.49
//Emitter 1 & Spectator A
void FI::SetDipKinematic1A(double &x, double &zi, double &zj, double &zplus, double &zminus, double mi, double mj, double mij) {

  x = 1 - (pow(mi,2) - pow(mij,2) + pow(mj,2) + 2*p1p3)/(2.*(pap1 + pap3));
  zi = pap3/(pap3 + pap1);
  zj = pap1/(pap3 + pap1);

  const double papb_temp = papb*x;
  const double pap1_temp = (pap1 + pap3)*x;
  const double pap2_temp = pap2*x;                         
  const double pbp1_temp = pbp1 + pbp3 + papb*(-1 + x);
  const double pbp2_temp = pbp2;

  const double mui = mi/sqrt(2.0 * (pap1 + pap3));
  const double muj = mj/sqrt(2.0 * (pap1 + pap3));
  const double muij = mij/sqrt(2.0 * (pap1 + pap3));


  papb = papb_temp;
  pap1 = pap1_temp;
  pap2 = pap2_temp;
  pbp1 = pbp1_temp;
  pbp2 = pbp2_temp;

  zplus = (1 - x + pow2(muij) + pow2(mui) - pow2(muj) + sqrt( pow2(1 - x + pow2(muij) - pow2(mui)) - 4.0 * pow2(mui) * pow2(muj) )) 
    / (2.0 * (1 - x + pow2(muij)));

  zminus = (1 - x + pow2(muij) + pow2(mui) - pow2(muj) - sqrt( pow2(1 - x + pow2(muij) - pow2(mui)) - 4.0 * pow2(mui) * pow2(muj) )) 
    / (2.0 * (1 - x + pow2(muij)));

}

//Emitter 1 & Spectator B
void FI::SetDipKinematic1B(double &x, double &zi, double &zj, double &zplus, double &zminus, double mi, double mj, double mij) {

  x = 1 - (pow(mi,2) - pow(mij,2) + pow(mj,2) + 2*p1p3)/(2.*(pbp1 + pbp3));   
  zi = pbp3/(pbp3 + pbp1);
  zj = pbp1/(pbp3 + pbp1);


  const double papb_temp = papb*x;
  const double pap1_temp = pap1 + pap3 + papb*(-1 + x);
  const double pap2_temp = pap2;                         
  const double pbp1_temp = (pbp1 + pbp3)*x;
  const double pbp2_temp = pbp2*x;

  const double mui = mi/sqrt(2.0 * (pbp1 + pbp3));
  const double muj = mj/sqrt(2.0 * (pbp1 + pbp3));
  const double muij = mij/sqrt(2.0 * (pbp1 + pbp3));


  papb = papb_temp;
  pap1 = pap1_temp;
  pap2 = pap2_temp;
  pbp1 = pbp1_temp;
  pbp2 = pbp2_temp;

  zplus = (1 - x + pow2(muij) + pow2(mui) - pow2(muj) + sqrt( pow2(1 - x + pow2(muij) - pow2(mui)) - 4.0 * pow2(mui) * pow2(muj) ))
    / (2.0 * (1 - x + pow2(muij)));

  zminus = (1 - x + pow2(muij) + pow2(mui) - pow2(muj) - sqrt( pow2(1 - x + pow2(muij) - pow2(mui)) - 4.0 * pow2(mui) * pow2(muj) )) 
    / (2.0 * (1 - x + pow2(muij)));

}


// still used for gauginos and leptons (same as above, but more simplifications used)
// (and harder to read)
void FI::SetDipKinematicA(double &x, double &fact) {
    x = 1.0 - (pap3 + pbp3) / papb;
    fact = 1.0 / x / pap3;

    papb *= x;
    const double L = 4.0 * papb + (1.0 - x) * pap3;
    pap1 = (pap3 * m1s + pap3 * p1p2 - 2.0 * papb * p1p3 + 2.0 * papb * pap1 + x * pap3 * m1s
            + x * pap3 * p1p3 + x * pap3 * p1p2 + 2.0 * x * papb * pap1) / L;
    pap2 = (pap3 * m2s + pap3 * p1p2 - 2.0 * papb * p2p3 + 2.0 * papb * pap2 + x * pap3 * m2s
            + x * pap3 * p2p3 + x * pap3 * p1p2 + 2.0 * x * papb * pap2) / L;
    pbp1 = (2.0 * papb * m1s + 2.0 * papb * p1p2 + 2.0 * papb * pbp1 - 2.0 * x * pap3 * m1s
            - x * pap3 * p1p3 - 2.0 * x * pap3 * p1p2 - 2.0 * x * papb * m1s - 2.0 * x * papb * p1p3
            - 2.0 * x * papb * p1p2 + 2.0 * x * papb * pbp1) / L;
    pbp2 = (2.0 * papb * m2s + 2.0 * papb * p1p2 + 2.0 * papb * pbp2 - 2.0 * x * pap3 * m2s
            - x * pap3 * p2p3 - 2.0 * x * pap3 * p1p2 - 2.0 * x * papb * m2s - 2.0 * x * papb * p2p3
            - 2.0 * x * papb * p1p2 + 2.0 * x * papb * pbp2) / L;

}


void FI::SetDipKinematicB(double &x, double &fact) {
    x = 1.0 - (pap3 + pbp3) / papb;
    fact = 1.0 / x / pbp3;

    papb *= x;
    const double L = 4.0 * papb + (1.0 - x) * pbp3;
    pbp1 = (pbp3 * m1s + pbp3 * p1p2 - 2.0 * papb * p1p3 + 2.0 * papb * pbp1 + x * pbp3 * m1s
            + x * pbp3 * p1p3 + x * pbp3 * p1p2 + 2.0 * x * papb * pbp1) / L;
    pbp2 = (pbp3 * m2s + pbp3 * p1p2 - 2.0 * papb * p2p3 + 2.0 * papb * pbp2 + x * pbp3 * m2s
            + x * pbp3 * p2p3 + x * pbp3 * p1p2 + 2.0 * x * papb * pbp2) / L;
    pap1 = (2.0 * papb * m1s + 2.0 * papb * p1p2 + 2.0 * papb * pap1 - 2.0 * x * pbp3 * m1s
            - x * pbp3 * p1p3 - 2.0 * x * pbp3 * p1p2 - 2.0 * x * papb * m1s - 2.0 * x * papb * p1p3
            - 2.0 * x * papb * p1p2 + 2.0 * x * papb * pap1) / L;
    pap2 = (2.0 * papb * m2s + 2.0 * papb * p1p2 + 2.0 * papb * pap2 - 2.0 * x * pbp3 * m2s
            - x * pbp3 * p2p3 - 2.0 * x * pbp3 * p1p2 - 2.0 * x * papb * m2s - 2.0 * x * papb * p2p3
            - 2.0 * x * papb * p1p2 + 2.0 * x * papb * pap2) / L;
}

// set the propagators
void FI::SetPropagator(const double mass1, const double mass2,
                       const double width1, const double width2) {
    static const complex<double> III(0.0, 1.0);

    // Precompute values which are used several times.
    double mass1s = pow2(mass1);
    double mass2s = pow2(mass2);
    complex<double> imass1width1 = III * mass1 * width1;
    complex<double> imass2width2 = III * mass2 * width2;

    // Electroweak s-channel propagators.
    ivs3v1 = 1.0 / (m1s + 2.0 * p1p2 + m2s - mass1s - imass1width1);
    ivs3v2 = 1.0 / (m1s + 2.0 * p1p2 + m2s - mass2s + imass2width2);

    // Electroweak t-channel propagators
    //Propagators without width. 
    ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mass1s);
    ivt2s1 = 1.0 / (m2s - 2.0 * pbp2 - mass1s);
    ivt1s2 = 1.0 / (m1s - 2.0 * pap1 - mass2s);
    ivt2s2 = 1.0 / (m2s - 2.0 * pbp2 - mass2s);


    // Electroweak u-channel propagators.
    // Propagators without width.
    ivu1s1 = 1.0 / (m1s - 2.0 * pbp1 - mass1s);
    ivu1s2 = 1.0 / (m1s - 2.0 * pbp1 - mass2s);
    ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mass1s);
    ivu2s2 = 1.0 / (m2s - 2.0 * pap2 - mass2s);

    // Strong quark propagators for real corrections.
    ivt3 = -0.5 / pap3;
    ivu3 = -0.5 / pbp3;
    ivs = 0.5 / papb;

    // Electroweak s-channel propagators in t-channel (and u-channel?) diagrams with real quark
    // emissions.
    ivs1s1 = 1.0 / (m2s + 2.0 * p2p3 - mass1s -  imass1width1);
    ivs2s1 = 1.0 / (m1s + 2.0 * p1p3 - mass1s -  imass1width1);
    ivs1s2 = 1.0 / (m2s + 2.0 * p2p3 - mass2s +  imass2width2);
    ivs2s2 = 1.0 / (m1s + 2.0 * p1p3 - mass2s +  imass2width2);

    // new propagatos for gaugl and final state emission
    iv13 = 0.5/p1p3;
    iv23 = 0.5/p2p3;
    
    iv13s1 = 1.0/(m1s - mass1s + 2.0*p1p3 - imass1width1);
    iv13s2 = 1.0/(m1s - mass2s + 2.0*p1p3 + imass2width2);

    iv23s1 = 1.0/(m2s - mass1s + 2.0*p2p3 - imass1width1);
    iv23s2 = 1.0/(m2s - mass2s + 2.0*p2p3 + imass2width2);

    ivb1 = -0.5/(pbp1);
    iva1 = -0.5/(pap1);

    iv12 = 1.0/(m1s + m2s + 2.0 * p1p2);

    // Electroweak s-channel propagators in t-channel diagrams with real quark
    // emissions.
    ivs1s1 = 1. / (m2s + 2.*p2p3 - std::pow(mass1, 2) -  mass1 * width1 * III);
    ivs2s1 = 1. / (m1s + 2.*p1p3 - std::pow(mass1, 2) -  mass1 * width1 * III);
    ivs1s2 = 1. / (m2s + 2.*p2p3 - std::pow(mass2, 2) +  mass2 * width2 * III);
    ivs2s2 = 1. / (m1s + 2.*p1p3 - std::pow(mass2, 2) +  mass2 * width2 * III);

}


// Breit Wigner form (squared Breit Wigner propagator)
// used for on-shell subtraction.
void FI::SetBreitWigner(double ps,double ms,double DecayWidth) {
    BreitWigner = ((ms * pow2(DecayWidth)) ) / ( pow2(ps - ms) + ms * pow2(DecayWidth));
}

// two couplings for the strong corrections
void FI::SetSCoupling(struct Coupling C[2]) {
    LL = C[0].L * C[1].L;
    LR = C[0].L * C[1].R;
    RR = C[0].R * C[1].R;
    RL = C[0].R * C[1].L;
}

// for gaugino-squark production at NLO
// couplings Ls,Lsp,Rs,Rsp, i.e. mathcal L(') etc
void FI::SetLsRsCoupling(struct Coupling C[2]) {
  LsLs = C[0].L * C[1].L; // Ls*Lsp
  LsRs = C[0].L * C[1].R; // Ls*Rsp
  RsRs = C[0].R * C[1].R; // Rs*Rsp
  RsLs = C[0].R * C[1].L; // Rs*Lsp
}

// for branching ratio I need the second one complex conj
void FI::SetBCoupling(struct Coupling C[2]) {
  LL = C[0].L * conj(C[1].L);
  LR = C[0].L * conj(C[1].R);
  RR = C[0].R * conj(C[1].R);
  RL = C[0].R * conj(C[1].L);
}



// Born couplings
// name is misleading for the new processes. They contain not only electroweak couplings.
void FI::SetWCoupling(struct Coupling C[4]) {
    // Several tests show that this function is critical in the execution time
    // of the program, so it is necessary to optimize it as much as possible.

    const complex<double> C2Lc = conj(C[2].L);
    const complex<double> C2Rc = conj(C[2].R);
    const complex<double> C3Lc = conj(C[3].L);
    const complex<double> C3Rc = conj(C[3].R);

    const complex<double> LL_1 = C[0].L * C[1].L;
    const complex<double> LR_1 = C[0].L * C[1].R;
    const complex<double> RL_1 = C[0].R * C[1].L;
    const complex<double> RR_1 = C[0].R * C[1].R;

    const complex<double> LL_2 = C2Lc * C3Lc;
    const complex<double> LR_2 = C2Lc * C3Rc;
    const complex<double> RL_2 = C2Rc * C3Lc;
    const complex<double> RR_2 = C2Rc * C3Rc;

    LLLL = LL_1 * LL_2;
    LLLR = LL_1 * LR_2;
    LLRL = LL_1 * RL_2;
    LRLL = LR_1 * LL_2;
    RLLL = RL_1 * LL_2;
    LLRR = LL_1 * RR_2;
    LRLR = LR_1 * LR_2;
    LRRL = LR_1 * RL_2;

    RRRR = RR_1 * RR_2;
    RRRL = RR_1 * RL_2;
    RRLR = RR_1 * LR_2;
    RLRR = RL_1 * RR_2;
    LRRR = LR_1 * RR_2;
    RRLL = RR_1 * LL_2;
    RLRL = RL_1 * RL_2;
    RLLR = RL_1 * LR_2;

}
