////////////////////////////////////////////////////////////////////////////////
//
//  vae.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef VAE_H
#define VAE_H
#include <cmath>
double czab(int Z, double a, double b);
double h(int Z, double a, double b, double c, double r);
double hp(int Z, double a, double b, double c, double r);
double hpp(int Z, double a, double b, double c, double r);
double phi(int Z, double a, double b, double c, double r);
double v(int Z, double a, double b, double c, double r);
void vsetb(int Z, double a, double& b);
double psnorm2(int Z, double a, double b);
#endif
