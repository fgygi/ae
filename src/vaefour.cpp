//
// vaefour.cpp
// Fourier transform of an AEPW potential
// Output:
//  plot of Vae(r)
//  plot of q^2 * Vae(q)
//
#include<iostream>
#include<vector>
#include<cmath>
using namespace std;
#include "sinft.h"
#include "vae.h"

// verf(Z,a,r) = -Z * erf(a*r)/r
double verf(int Z, double a, double r)
{
  if ( r == 0.0 )
    // M_2_SQRTPI = 2 / sqrt(pi)
    return -a * Z * M_2_SQRTPI;
  return -Z*erf(a*r)/r;
}

double simpsn ( int n, double *t );

int main(int argc, char **argv)
{
  if ( argc < 2 )
  {
    cerr << " Use: vaefour a [b [c] ]" << endl;
    return 1;
  }
  int np = 1000;
  const int Z = 1;
  double a = atof(argv[1]);
  double h = 0.02/a;
  double b = -1.0/(a*sqrt(M_PI));
  if ( argc > 2 )
    b = atof(argv[2]);
  else
    vsetb(a,b);
  double c = cab(a,b);
  if ( argc > 3 )
    c = atof(argv[3]);

  cerr << "Z=" << Z << " a=" << a << " b=" << b << " c=" << c << endl;

  // plot Vae(r)
  for ( int i = 0; i < np; i++ )
  {
    double r = h * i;
    cout << r << " " << v(Z,a,b,c,r) << endl;
  }
  cout << endl << endl;

  vector<double> f(np),fint(np);
  // compute Fourier ( vae(r)-verf(r) ) + Fourier ( verf(r) )
  for ( int i = 0; i < np; i++ )
  {
    double r = h * i;
    f[i] = 4.0 * M_PI * r * h * ( v(Z,a,b,c,r) - verf(Z,a,r) );
  }
  // q=0 integral
  for ( int i = 0; i < np; i++ )
  {
    double r = h * i;
    fint[i] = f[i] * r;
  }
  double fq0 = simpsn(np,&fint[0]);

  sinft(np,&f[0]);

  // plot q^2 * ( V(q) + Verf(q) )
  cout << 0 << " " << fq0 << endl;
  for ( int i = 0; i < np; i++ )
  {
    double q = i * M_PI / (h * np);
    //cout << q << " " << f[i]/q << endl;
    double q2_v_q = q * f[i];
    double q2_verf_q = 4 * M_PI * exp(-q*q/(4*a*a));
    cout << q << " " << q2_v_q + q2_verf_q << endl;
  }

#if 0
  // reference gaussian
  // v(r) = exp(-a^2 r^2)
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
