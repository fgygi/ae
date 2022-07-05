//
// gfour.cpp
// Fourier transform of a 3d gaussian function
//
#include<iostream>
#include<vector>
#include<cmath>
using namespace std;
#include "sinft.h"

// f(r) = exp(-a^2 * r^2)
const double a = 1.5;

double f(double r)
{
  const double a = 1.5;
  return exp(-a*a*r*r);
}

double simpsn ( int n, double *t );

int main()
{
  int np = 500;
  double h = 0.02;

  vector<double> v(np),fint(np);

  // q=0 integral
  for ( int i = 0; i < np; i++ )
  {
    double r = h * i;
    fint[i] = 4.0 * M_PI * r * r * f(r);
  }
  double vq0 = h * simpsn(np,&fint[0]);

  // v(q)
  for ( int i = 0; i < np; i++ )
  {
    double r = h * i;
    v[i] = 4.0 * M_PI * r * f(r) * h;
  }
  for ( int i = 0; i < np; i++ )
    cout << i*h << " " << v[i] << endl;
  cout << endl << endl;

  sinft(np,&v[0]);

  cout << 0 << " " << vq0 << endl;
  for ( int i = 1; i < np; i++ )
  {
    double q = i * M_PI / (h * np);
    cout << q << " " << v[i]/q << endl;
  }

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
