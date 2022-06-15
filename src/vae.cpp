////////////////////////////////////////////////////////////////////////////////
//
// vae.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "vae.h"
#include "RootFinder.h"
#include <cmath>
#include <cassert>
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
double czab(int Z, double a, double b)
{
  // find value of c that enforces zero-curvature of V at r=0
  const double x = Z/(a*sqrt(M_PI));
  // compute both roots cm, cp
  const double s = sqrt( 25.0 + 10.0 * b + 80.0 * x / 3.0 ) ;
  const double cm = a * a * ( 2.5 + b + 2.0 * x - 0.5 * s );
  const double cp = a * a * ( 2.5 + b + 2.0 * x + 0.5 * s );
  assert( fabs(cm) < fabs(cp) );
  return cm;
}

////////////////////////////////////////////////////////////////////////////////
double h(int Z, double a, double b, double c, double r)
{
  return -Z*r*erf(a*r) + ( b + c * r*r ) * exp(-a*a*r*r);
}

////////////////////////////////////////////////////////////////////////////////
double hp(int Z, double a, double b, double c, double r)
{
  const double a2r2 = a*a*r*r;
  return -Z*erf(a*r) + 2.0*( c*(1.0-a2r2) - a*a*b - a*Z/sqrt(M_PI) ) *
          r * exp(-a2r2);
}

////////////////////////////////////////////////////////////////////////////////
double hpp(int Z, double a, double b, double c, double r)
{
  double a2 = a*a;
  double a4 = a2*a2;
  double a2r2 = a2*r*r;
  const double x = Z/(a*sqrt(M_PI));
  return ( 2.0*c - 2.0*a2*b - 4.0*a2*x +
           ( 4.0*a4*b-10.0*a2*c+4.0*a4*x )*r*r +
           4.0*a2r2*a2r2*c ) * exp(-a2r2);
}

////////////////////////////////////////////////////////////////////////////////
double phi(int Z, double a, double b, double c, double r)
{
  return sqrt(Z*Z*Z/M_PI)*exp(h(Z,a,b,c,r));
}

////////////////////////////////////////////////////////////////////////////////
double v(int Z, double a, double b, double c, double r)
{
  if ( r == 0.0 )
  {
    return -0.5*Z*Z + 3.0*c - 3.0*a*a*b - 6.0*a*Z/sqrt(M_PI);
  }
  else
  {
    double hpval = hp(Z,a,b,c,r);
    return -0.5*Z*Z + hpval/r + 0.5*hpval*hpval + 0.5*hpp(Z,a,b,c,r);
  }
}

////////////////////////////////////////////////////////////////////////////////
void vsetb(int Z, double a, double& b)
{
  // Adjust b to enforce norm conservation
  RootFinder finder(0.1);
  int iter = 0;
  const double epsilon = 1.e-6;
  double pseudo_norm2 = psnorm2(Z,a,b);
  while ( fabs(pseudo_norm2-1.0) > epsilon )
  {
    iter++;
    assert(iter<200);
    b = finder.next(b,pseudo_norm2-1.0);
    pseudo_norm2 = psnorm2(Z,a,b);
  }
}

////////////////////////////////////////////////////////////////////////////////
// pseudofunction 2-norm
double psnorm2(int Z, double a, double b)
{
  double c = czab(Z,a,b);
  int np = 501;
  const double dr = 0.002/Z;
  int npout = 20 * np;
  double sum = 0.0;
  for ( int i = 0; i < npout; i++ )
  {
    double r = dr * i;
    double t = r * phi(Z,a,b,c,r);
    sum += t*t*dr;
  }
  return 4 * M_PI * sum;
}
