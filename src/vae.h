////////////////////////////////////////////////////////////////////////////////
//
//  vae.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef VAE_H
#define VAE_H
const char *const version = "v4.1";
#include <cmath>
double cab(double a, double b);
double h1(double a, double b, double c, double r);
double h1p(double a, double b, double c, double r);
double h1pp(double a, double b, double c, double r);
double phi(int Z, double a, double b, double c, double r);
double phi1(double a, double b, double c, double r);
double h(int Z, double a, double b, double c, double r);
double hp(int Z, double a, double b, double c, double r);
double hpp(int Z, double a, double b, double c, double r);
double phi(int Z, double a, double b, double c, double r);
double v1(double a, double b, double c, double r);
double v(int Z, double a, double b, double c, double r);
void vsetb(double a, double& b);
double psnorm2(double a, double b);
void vsetbc(double a, double& b, double& c);
double psnorm2_c0(double a, double b);
#endif
