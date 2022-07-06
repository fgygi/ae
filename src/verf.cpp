////////////////////////////////////////////////////////////////////////////////
//
// erf(ar)/r hydrogen pseudopotential
//
// use: verf Z a
//
////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<cassert>
#include "RootFinder.h"
using namespace std;

int Z;
double a;

double v(double r)
{
  if ( r == 0.0 )
    // M_2_SQRTPI = 2 / sqrt(pi)
    return -a * Z * M_2_SQRTPI;
  return -Z*erf(a*r)/r;
}

int main(int argc, char **argv)
{
  if ( !(argc == 3) )
  {
    cerr << "Use: verf Z a" << endl;
    return 1;
  }
  Z = atoi(argv[1]);
  a = atof(argv[2]);
  const double dr = 0.002/Z;
  int np = 501;
  // adjust np so that v(r) = -Z/r at r=(np-)*dr
  while ( fabs(v((np-1)*dr)+Z/((np-1)*dr)) > 1.e-10 )
  {
    np += 100;
    assert(np<=10001);
  }

  // output XML potential file
  cout.setf(ios::scientific,ios::floatfield);
  cout << setprecision(10);

  cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  cout << "<fpmd:species xmlns:fpmd="
          "\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0\"" << endl;
  cout << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << endl;
  cout << "xsi:schemaLocation=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0 species.xsd\">" << endl;
  cout << "<description>" << endl;
  cout << "  Verf all-electron potential" << endl;
  cout << "  Z=" << Z << " a=" << a << " dr=" << dr << endl;
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
    cout << v(r) << endl;
  }
  cout << "</radial_potential>" << endl;
  cout << "</projector>" << endl;
  cout << "</norm_conserving_pseudopotential>" << endl;
  cout << "</fpmd:species>" << endl;
}
