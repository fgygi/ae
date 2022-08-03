////////////////////////////////////////////////////////////////////////////////
//
//  aecheck.cpp
//  Solve the non-interacting problem in a regularized potential
//
//  Use: aecheck [-zc] Z, a, l, rmax, hZ
//  If -zc is used, zero-curvature of V(r) is enforced at r=0
//  Z: atomic number
//  a parameter of the regularized potential
//  The b parameter is calculated by enforcing norm conservation of the 1s state
//  l: angular momentum
//  hZ: product of mesh spacing times Z (typical values: 0.0002 - 0.002)
//
////////////////////////////////////////////////////////////////////////////////

#include "vae.h"
#include "nradsolve.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <valarray>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

const char *const version = "v4.0";

int main(int argc, char **argv)
{
  // use: aecheck [-zc] Z a l rmax hZ
  if ( !(argc == 6 || argc == 7) )
  {
    cerr << "aecheck " << version << endl;
    cerr << "Use: aecheck [-zc] Z a l rmax hZ" << endl;
    return 1;
  }
  bool zc = false;
  double c = 0.0;
  int iarg = 1;
  if ( !strcmp(argv[iarg],"-zc") )
  {
    assert(argc==7);
    // impose zero-curvature condition at r=0
    zc = true;
    iarg++;
  }
  int Z = atoi(argv[iarg++]);
  double a = atof(argv[iarg++]);
  // initial value of b
  double b = -1.0/(a*sqrt(M_PI));

  if ( zc )
  {
    vsetbc(a,b,c);
  }
  else
  {
    vsetb(a,b);
  }

  int l = atoi(argv[iarg++]);
  double rmax = atof(argv[iarg++]);
  const double hZ = atof(argv[iarg++]);

  const double h = hZ/Z;
  const int np = (int)(rmax/h);
  rmax = h*(np+1);

  cout << " aecheck " << version << endl;
  cout << " Z = " << Z << " a = " << a << " b = " << b << " c = " << c
       << " l = " << l << endl;
  cout << " rmax = " << rmax << " hZ = " << h*Z << endl;
  cout << " np = " << np << " h = " << h << endl;

  valarray<double> r(np);
  valarray<double> vext(np);
  for ( int i = 0; i < np; i++ )
  {
    r[i] = (i+1)*h;
    vext[i] = v(Z,a,b,c,r[i]);
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
