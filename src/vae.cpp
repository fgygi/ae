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
  const double x = a*Z/sqrt(M_PI);
  return 2.5*a*a + a*a*b + 2.0*x -
         0.5*sqrt( (75.0*a*a*a*a + 30.0*a*a*a*a*b + 80.0*a*a*x ) / 3.0 );
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
  const double x = a*Z/sqrt(M_PI);
  return -Z*erf(a*r) + 2.0*( c*(1.0-a2r2) - a*a*b - a*Z/sqrt(M_PI) ) *
          r * exp(-a2r2);
}

////////////////////////////////////////////////////////////////////////////////
double hpp(int Z, double a, double b, double c, double r)
{
  double a2r2 = a*a*r*r;
  const double x = a*Z/sqrt(M_PI);
  return ( 2.0*c - 2.0*a*a*b - 4.0*x +
           ( 4.0*a*a*a*a*b-10.0*a*a*c+4.0*a*a*x )*r*r +
           4.0*a2r2*a2r2*c ) * exp(-a2r2);
}

////////////////////////////////////////////////////////////////////////////////
double phi(int Z, double a, double b, double c, double r)
{
  return exp(h(Z,a,b,c,r));
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
// adjust b to enforce norm conservation in [1/Z, Infinity]
// return the 2-norm of the phi_10 function for parameters Z,a,b
double vsetb(int Z, double a, double& b)
{
  b = -Z/(a*sqrt(M_PI));
  double c = czab(Z,a,b);
  const double dr = 0.002/Z;
  int np = 501;
  // adjust np so that v(r) = -Z/r at r=(np-)*dr
  while ( fabs(v(Z,a,b,c,(np-1)*dr)+Z/((np-1)*dr)) > 1.e-10 )
  {
    np += 100;
    assert(np<=10001);
  }

  // compute norm of exact solution outside of r=1/Z
  // Normalized exact solution is Z^(3/2)/sqrt(pi) exp(-Z*r)
  double npn = 1.0/(Z*dr);
  int npout = 20 * np;
  double r = dr*npn;
  double t = r * (sqrt(Z*Z*Z)/sqrt(M_PI)) * exp(-Z*r);
  double exact_norm2 = 0.5*t*t*dr;
  for ( int i = npn+1; i < npout; i++ )
  {
    r = dr * i;
    t = r * (sqrt(Z*Z*Z)/sqrt(M_PI)) * exp(-Z*r);
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
    for ( int i = 0; i < npout; i++ )
    {
      r = dr * i;
      t = r * phi(Z,a,b,c,r);
      sum += t*t*dr;
    }
    fac = 1.0 / sqrt(4.0*M_PI*sum);

    // compute norm2 of normalized phi outside of r=1/Z
    // trapezoidal rule: first point with weight 0.5
    r = dr*npn;
    t = r * fac * phi(Z,a,b,c,r);
    pseudo_norm2 = 0.5*t*t*dr;
    for ( int i = npn+1; i < npout; i++ )
    {
      r = dr * i;
      t = r * fac * phi(Z,a,b,c,r);
      pseudo_norm2 += t*t*dr;
    }
    pseudo_norm2 *= 4 * M_PI;
    //cerr << "pseudo_norm2-exact_norm2="
    //     << pseudo_norm2-exact_norm2 << endl;
    if ( fabs(pseudo_norm2-exact_norm2) > epsilon )
      b = finder.next(b,pseudo_norm2-exact_norm2);
      c = czab(Z,a,b);
  }
  cerr << "Pseudo norm2 in [1/Z,infinity]: " << pseudo_norm2
       << " b=" << b << endl;
  return fac;
}

////////////////////////////////////////////////////////////////////////////////
// norm2 of pseudofunction outside of 1/Z
double psnorm2(int Z, double a, double b)
{
  double c = czab(Z,a,b);
  int np = 501;
  const double dr = 0.002/Z;
  int npout = 20 * np;
  // compute normalization factor of phi(r)
  double r,t;
  double sum = 0.0;
  for ( int i = 0; i < npout; i++ )
  {
    r = dr * i;
    t = r * phi(Z,a,b,c,r);
    sum += t*t*dr;
  }
  const double fac = 1.0 / sqrt(4.0*M_PI*sum);

  // compute norm2 of normalized phi outside of r=1/Z
  double npn = 1.0/(Z*dr);
  r = npn*dr;
  t = r * fac * phi(Z,a,b,c,r);
  double pseudo_norm2 = 0.5*t*t*dr;
  for ( int i = npn+1; i < npout; i++ )
  {
    r = dr * i;
    t = r * fac * phi(Z,a,b,c,r);
    pseudo_norm2 += t*t*dr;
  }
  pseudo_norm2 *= 4 * M_PI;
  return pseudo_norm2;
}
