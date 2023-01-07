////////////////////////////////////////////////////////////////////////////////
//
//  aecheck_erfgau.cpp
//  Solve the non-interacting problem in an erfgau potential
//  [1] C. E. Gonzalez-Espinoza et al, Theor. Chem. Acc 135, 256 (2016)
//
//  Note: the parameter a is mu in Ref [1]
//  Input: Z, a, l, rmax, hZ
//  Z: atomic number
//  a parameters of the regularized potential
//  l: angular momentum
//  hZ: produce of mesh spacing times Z (typical values: 0.0002 - 0.002)
//
////////////////////////////////////////////////////////////////////////////////

#include "nradsolve.h"
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <valarray>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

double verfgau(int Z, double mu, double r)
{
  const double kappa = 1.5;
  const double spi = sqrt(M_PI);
  const double gamma = 27.0 / (8.0 * spi);
  const double alpha0 = (8.0 + sqrt(64.0-18.0*M_PI))/(12.0 * spi);
  const double c0 = 27.0 * alpha0 / (4.0 * spi);
  const double alpha = kappa * mu + alpha0;
  const double c = gamma * mu + c0;
  assert(Z==1.0);
  if ( r == 0.0 )
    return - (c + 2.0 * mu * Z / sqrt(M_PI));
  else
    return - (c * exp(-alpha*alpha*r*r) + erf(mu*r)/r);
}

int main(int argc, char **argv)
{
  // use: aecheck_erf Z a l rmax hZ
  if ( !(argc == 6 ) )
  {
    cerr << " use: aecheck_erf Z a l rmax hZ" << endl;
    return 1;
  }
  int iarg = 1;
  int Z = atoi(argv[iarg++]);
  double a = atof(argv[iarg++]);
  int l = atoi(argv[iarg++]);
  double rmax = atof(argv[iarg++]);
  const double hZ = atof(argv[iarg++]);

  const double h = hZ/Z;
  const int np = (int)(rmax/h);
  rmax = h*(np+1);
  cout << " Z = " << Z << " a = " << a
       << " l = " << l << endl;
  cout << " rmax = " << rmax << " hZ = " << h*Z << endl;
  cout << " np = " << np << " h = " << h << endl;

  valarray<double> r(np);
  valarray<double> vext(np);
  for ( int i = 0; i < np; i++ )
  {
    r[i] = (i+1)*h;
    vext[i] = verfgau(Z,a,r[i]);
  }

  //////////////////////////////////////////////////////////////////////////////
  // solve radial Schroedinger equation for potential v

  int neig = 7-l;
  valarray<double> eig(neig);
  valarray<double>u(neig*np);

  nradsolve(rmax,l,neig,vext,u,eig);

  cout << setprecision(8);
  cout.setf(ios::fixed,ios::floatfield);
  cout.setf(ios::right,ios::adjustfield);
  for ( int i = 0; i < neig; i++ )
  {
    int n = i+l+1;
    cout << " E_" << n << " = " << setw(15) << eig[i]
         << "   " << setw(15) << -(0.5*Z*Z)/(n*n)
         << "   " << setw(15) << eig[i]+(0.5*Z*Z)/(n*n) << endl;
  }

  //////////////////////////////////////////////////////////////////////////////
  // plot output
#if 0
  ofstream outfile("aecheck.dat");
  for ( int ie = 0; ie < neig; ie++ )
  {
    int n = ie+1;
    outfile << "# u(r) i=" << n << " E = " << eig[ie] << endl;
    outfile << "0.0  0.0" << endl;
    for ( int i = 0; i < np; i++ )
      outfile << r[i] << " " << u[i+ie*np] << endl;
    outfile << endl << endl;
  }
  outfile.close();
#endif

  return 0;
}
