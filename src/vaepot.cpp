////////////////////////////////////////////////////////////////////////////////
//
// generate analytic hydrogen pseudopotential
//
// use: v Z a [b]
//
// If not provided, the parameter b is determined by enforcing
// norm-conservation at r=1/Z
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

int main(int argc, char **argv)
{
  if ( !(argc == 3 || argc == 4) )
  {
    cerr << "Use: v Z a [b]" << endl;
    return 1;
  }
  int Z = atoi(argv[1]);
  double a = atof(argv[2]);
  double b = 0.0;
  if ( argc == 4 )
    b = atof(argv[3]);

  double fac = vsetb(Z,a,b);

  // output XML potential file
  const double dr = 0.002/Z;
  int np = 501;

  cout.setf(ios::scientific,ios::floatfield);
  cout << setprecision(10);

  cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  cout << "<fpmd:species xmlns:fpmd="
          "\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0\"" << endl;
  cout << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << endl;
  cout << "xsi:schemaLocation=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0 species.xsd\">" << endl;
  cout << "<description>" << endl;
  cout << "  Analytic regularized all-electron potential" << endl;
  cout << "  Z=" << Z << " a=" << a << " b=" << b << " dr=" << dr << endl;
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
    cout << v(Z,a,b,r) << endl;
  }
  cout << "</radial_potential>" << endl;

  cout << "<radial_function>" << endl;
  for ( int i = 0; i < np; i++ )
  {
    double r = dr * i;
    cout << fac*phi(Z,a,b,r) << endl;
  }
  cout << "</radial_function>" << endl;
  cout << "</projector>" << endl;
  cout << "</norm_conserving_pseudopotential>" << endl;
  cout << "</fpmd:species>" << endl;
}
