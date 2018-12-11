// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2013 David R. Lamprea.
// Copyright 2011-2013 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// N-point functions module.

#include <complex>

using namespace std;

// External computation of scalar one-loop integrals.
extern "C" {
    void qlinit_();

    // TODO: complex<double> is in principle compatible with Fortran complex
    // but this can be rewritten in C-style.
    complex<double> qli1_(double &m1, double &mu2, int &ep);
    complex<double> qli2_(double &p1, double &m1, double &m2,
                          double &mu2, int &ep);
    complex<double> qli3_(double &p1, double &p2, double &p3, double &m1,
                          double &m2, double &m3, double &mu2, int &ep);
    complex<double> qli4_(double &p1, double &p2, double &p3,
                          double &p4, double &s12, double &s23,
                          double &m1, double &m2, double &m3,
                          double &m4, double &mu2, int &ep);
}

// Declaration of class arguments and functions.
class Npf
{
public:
    // One-point function constructor
    Npf(double m1, double Scale, int Index) {
        qlinit_();
        npt = 10;
        mu = Scale;
        ieps = Index;
        SetArgument(m1);
    };
    // Two-point function constructor
    Npf(double p1, double m1, double m2, double Scale, int Index) {
        qlinit_();
        npt = 20;
        mu = Scale;
        ieps = Index;
        SetArgument(p1, m1, m2);
    };
    // Three-point function constructor
    Npf(double p1, double p2, double p3, double m1, double m2, double m3,
        double Scale, int Index) {
        qlinit_();
        npt = 30;
        mu = Scale;
        ieps = Index;
        SetArgument(p1, p2, p3, m1, m2, m3);
    };
    // Four-point function constructor
    Npf(double p1, double p2, double p3, double p4, double s12, double s23,
        double m1, double m2, double m3, double m4 , double Scale, int Index) {
        qlinit_();
        npt = 40;
        mu = Scale;
        ieps = Index;
        SetArgument(p1, p2, p3, p4, s12, s23, m1, m2, m3, m4);
    };
    // N-point function destructor
    ~Npf() {};

    // Allow to change the Laurent coeff. you get
    void SetIeps(int Index);
    // Get what you ask for!
    void GetNpf(complex<double> Array[7], int Index);

private:
    int npt, ieps;
    double mu;

    double m1a;
    double p1b, m1b, m2b;
    double p1c, p2c, p3c, m1c, m2c, m3c;
    double p1d, p2d, p3d, p4d, s12d, s23d, m1d, m2d, m3d, m4d;

    double f1b;
    double f1c, f2c;
    double f1d, f2d, f3d;

    double Xc[2][2];
    double Xd[3][3];

    void SetArgument(const double m1);
    void SetArgument(const double p1, const double m1, const double m2);
    void SetArgument(const double p1, const double p2, const double p3,
                     const double m1, const double m2, const double m3);
    void SetArgument(const double p1, const double p2, const double p3,
                     const double p4, const double s12, const double s23,
                     const double m1, const double m2, const double m3,
                     const double m4);

    void GetA0(complex<double> &a0);
    void GetB0(complex<double> &b0);
    void GetB1(complex<double> *b1);
    void GetB2(complex<double> *b2);
    void GetC0(complex<double> &c0);
    void GetC1(complex<double> *c1);
    void GetC2(complex<double> *c2);
    void GetD0(complex<double> &d0);
    void GetD1(complex<double> *d1);
    void GetD2(complex<double> *d2);
};
