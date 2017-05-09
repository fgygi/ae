////////////////////////////////////////////////////////////////////////////////
//
//  vae.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef VAE_H
#define VAE_H
#include <cmath>
double h(int Z, double a, double b, double r);
double hp(int Z, double a, double b, double r);
double hpp(int Z, double a, double b, double r);
double phi(int Z, double a, double b, double r);
double v(int Z, double a, double b, double r);
double vsetb(int Z, double a, double& b);
#endif
