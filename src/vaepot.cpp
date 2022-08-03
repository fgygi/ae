////////////////////////////////////////////////////////////////////////////////
//
// Generate analytic regularized Coulomb pseudopotential
//
// use: vaepot [-zc] Z a
//
// The parameter b is determined by enforcing norm-conservation.
// The default value of parameter c is 0
// If the -zc option is used, the parameter c is determined by enforcing
// zero curvature of V(r) at r=0.
//
////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<cassert>
#include "RootFinder.h"
#include "vae.h"
using namespace std;

const char *const version = "v4.0";

int main(int argc, char **argv)
{
  if ( !(argc == 3 || argc == 4) )
  {
    cerr << "vaepot " << version << endl;
    cerr << "Use: vaepot [-zc] Z a" << endl;
    return 1;
  }
  bool zc = false;
  double c = 0.0;
  int iarg = 1;
  if ( !strcmp(argv[iarg],"-zc") )
  {
    assert(argc==4);
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

  cerr << "vaepot " << version << endl;
  cerr << "Z=" << Z << " a=" << a << " b=" << b << " c=" << c << endl;
  const double dr = 0.002/Z;
  int np = 501;
  // adjust np so that v(r) = -Z/r at r=(np-)*dr
  while ( fabs(v(Z,a,b,c,(np-1)*dr)+Z/((np-1)*dr)) > 1.e-10 )
  {
    np += 100;
    assert(np<=10001);
  }
  cerr << "rmax=" << (np-1)*dr << " np=" << np << endl;

  // output XML potential file
  cout.setf(ios::scientific,ios::floatfield);
  cout << setprecision(10);

  cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  cout << "<fpmd:species xmlns:fpmd="
          "\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0\"" << endl;
  cout << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << endl;
  cout << "xsi:schemaLocation=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0 species.xsd\">" << endl;
  cout << "<description>" << endl;
  cout << "  Analytic regularized all-electron potential " << version << endl;
  cout << "  Z=" << Z << " a=" << a << " b=" << b << " c=" << c
       << " dr=" << dr << endl;
  cout << "</description>" << endl;
  cout << "<symbol>X</symbol>" << endl;
  cout << "<atomic_number>" << Z << "</atomic_number>" << endl;
  cout << "<mass>1.00</mass>" << endl;
  cout << "<norm_conserving_pseudopotential>" << endl;
  cout << "<valence_charge>" << Z << "</valence_charge>" << endl;
  cout << "<lmax>0</lmax>" << endl;
  cout << "<llocal>0</llocal>" << endl;
  cout << "<nquad>0</nquad>" << endl;
  cout << "<rquad>0</rquad>" << endl;
  cout << "<mesh_spacing> " << dr << " </mesh_spacing>" << endl;
  cout << "<projector l=\"0\" size=\"" << np << "\">" << endl;
  cout << "<radial_potential>" << endl;
  for ( int i = 0; i < np; i++ )
  {
    double r = dr * i;
    cout << v(Z,a,b,c,r) << endl;
  }
  cout << "</radial_potential>" << endl;

  cout << "<radial_function>" << endl;
  for ( int i = 0; i < np; i++ )
  {
    double r = dr * i;
    cout << phi(Z,a,b,c,r) << endl;
  }
  cout << "</radial_function>" << endl;
  cout << "</projector>" << endl;
  cout << "</norm_conserving_pseudopotential>" << endl;
  cout << "</fpmd:species>" << endl;
}
