//
// vaefour.cpp
// Fourier transform of an AE potential
//
#include<iostream>
#include<vector>
#include<cmath>
using namespace std;
#include "sinft.h"
#include "vae.h"

// f(r) = exp(-a^2 * r^2)
const double a = 1.5;

double f(double r)
{
  const double a = 1.5;
  return exp(-a*a*r*r);
}

double simpsn ( int n, double *t );

int main(int argc, char **argv)
{
  int np = 500;
  double h = 0.02;
  double a = atof(argv[1]);
  double b,c;
  double b = -1.0/(a*sqrt(M_PI));
  vsetb(a,b);
  c = cab(a,b);

  vector<double> f(np),fint(np);

  // q=0 integral
  for ( int i = 0; i < np; i++ )
  {
    double r = h * i;
    fint[i] = 4.0 * M_PI * r * r * phi1(a,b,c,r);
  }
  double fq0 = h * simpsn(np,&fint[0]);

  // f(q)
  for ( int i = 0; i < np; i++ )
  {
    double r = h * i;
    f[i] = 4.0 * M_PI * r * phi1(a,b,c,r) * h;
  }
  for ( int i = 0; i < np; i++ )
    cout << i*h << " " << f[i] << endl;
  cout << endl << endl;

  sinft(np,&f[0]);

  cout << 0 << " " << fq0 << endl;
  for ( int i = 1; i < np; i++ )
  {
    double q = i * M_PI / (h * np);
    cout << q << " " << f[i]/q << endl;
  }

#if 0
  // reference
  // v(q) = v(q=0)*exp(-q^2/(4 a^2))
  // v(q=0) = 4 pi int r^2 exp(-a^2 r^2) dr = pi^(3/2) / a^3
  cout << endl << endl;
  const double vq0ref = pow(M_PI,1.5)/(a*a*a);
  for ( int i = 0; i < np; i++ )
  {
    double q = i * M_PI / (h * np);
    cout << q << " " << vq0ref * exp(-q*q/(4*a*a)) << endl;
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
double simpsn ( int n, double *t )
{
  const double c0 =  17.0/48.0, c1 = 59.0/48.0,
               c2 =  43.0/48.0, c3 = 49.0/48.0;
  double sum = c0 * ( t[0] + t[n-1] ) +
               c1 * ( t[1] + t[n-2] ) +
               c2 * ( t[2] + t[n-3] ) +
               c3 * ( t[3] + t[n-4] );

  for ( int i = 4; i < n-4; i++ )
  {
    sum += t[i];
  }
  return sum;
}
