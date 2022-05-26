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

const double xi = 1.0;
// M_1_SQRTPI = 1/sqrt(pi)

////////////////////////////////////////////////////////////////////////////////
double h(int Z, double a, double b, double r)
{
  double a2r2 = a*a*r*r;
  return -Z*r*erf(a*r) - xi*Z/(a*sqrt(M_PI)) * exp(-a2r2)
         + b*r*r*exp(-a2r2);
}

////////////////////////////////////////////////////////////////////////////////
double hp(int Z, double a, double b, double r)
{
  double a2r2 = a*a*r*r;
  return -Z*erf(a*r) - 2.0*(a*Z*r/sqrt(M_PI)) * (1.0 - xi) * exp(-a2r2)
         + (2.0*b*r*(1.0-a2r2)) * exp(-a2r2);
}

////////////////////////////////////////////////////////////////////////////////
double hpp(int Z, double a, double b, double r)
{
  double a2r2 = a*a*r*r;
  return   -2*Z*(a/sqrt(M_PI)) * ( 2.0 + 2.0*a2r2*(xi-1.0) -xi ) * exp(-a2r2)
         + (2.0*b*(1.0-5.0*a2r2+2.0*a2r2*a2r2)) * exp(-a2r2);
}

////////////////////////////////////////////////////////////////////////////////
double phi(int Z, double a, double b, double r)
{
  return exp(h(Z,a,b,r));
}

////////////////////////////////////////////////////////////////////////////////
double v(int Z, double a, double b, double r)
{
  if ( r == 0.0 )
  {
    return -0.5*Z*Z + 3.0*a*Z/sqrt(M_PI) * (xi - 2.0) + 3.0*b;
  }
  else
  {
    double hpval = hp(Z,a,b,r);
    return -0.5*Z*Z + hpval/r + 0.5*hpval*hpval + 0.5*hpp(Z,a,b,r);
  }
}

////////////////////////////////////////////////////////////////////////////////
// adjust b to enforce norm conservation in [1/Z, Infinity]
// return the 2-norm of the phi_10 function for parameters Z,a,b
double vsetb(int Z, double a, double& b)
{
  b = 0.0;
  const double dr = 0.002/Z;
  int np = 501;
  // adjust np so that v(r) = -Z/r at r=(np-)*dr
  while ( fabs(v(Z,a,b,(np-1)*dr)+Z/((np-1)*dr)) > 1.e-10 )
  {
    np += 100;
    assert(np<=10001);
  }

  // compute norm of exact solution outside of r=1/Z
  // Normalized exact solution is Z^(3/2)/sqrt(pi) exp(-Z*r)
  double npn = 1.0/(Z*dr);
  double exact_norm2 = 0.0;
  for ( int i = npn; i < 10*np; i++ )
  {
    double r = dr * i;
    double t = r * (sqrt(Z*Z*Z)/sqrt(M_PI)) * exp(-Z*r);
    exact_norm2 += t*t*dr;
  }
  exact_norm2 *= 4 * M_PI;
  cerr << "Exact  norm2 in [1/Z,infinity]: " << exact_norm2 << endl;

  // Adjust b to enforce norm conservation at r=1/Z
  RootFinder finder(0.5);
  double pseudo_norm2=0,fac;
  int iter = 0;
  const double epsilon = 1.e-6;
  while ( fabs(pseudo_norm2-exact_norm2) > epsilon )
  {
    iter++;
    assert(iter<200);
    // compute normalization factor of phi(r)
    double sum = 0.0;
    for ( int i = 0; i < 10*np; i++ )
    {
      double r = dr * i;
      double t = r * phi(Z,a,b,r);
      sum += t*t*dr;
    }
    fac = 1.0 / sqrt(4.0*M_PI*sum);

    // compute norm2 of normalized phi outside of r=1/Z
    pseudo_norm2 = 0.0;
    for ( int i = npn; i < 10*np; i++ )
    {
      double r = dr * i;
      double t = r * fac * phi(Z,a,b,r);
      pseudo_norm2 += t*t*dr;
    }
    pseudo_norm2 *= 4 * M_PI;
    if ( fabs(pseudo_norm2-exact_norm2) > epsilon )
      b = finder.next(b,pseudo_norm2-exact_norm2);
  }
  cerr << "Pseudo norm2 in [1/Z,infinity]: " << pseudo_norm2
       << " b=" << b << endl;
  return fac;
}
