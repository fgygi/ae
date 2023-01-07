////////////////////////////////////////////////////////////////////////////////
//
// erfgau hydrogen pseudopotential
//
// use: verfgau a
//
////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<cassert>
#include "RootFinder.h"
using namespace std;

const int Z = 1;

double v(double mu, double r)
{
  const double kappa = 1.5;
  const double spi = sqrt(M_PI);
  const double gamma = 27.0 / (8.0 * spi);
  const double alpha0 = (8.0 + sqrt(64.0-18.0*M_PI))/(12.0 * spi);
  const double c0 = 27.0 * alpha0 / (4.0 * spi);
  const double alpha = kappa * mu + alpha0;
  const double c = gamma * mu + c0;
  if ( r == 0.0 )
    return - (c + 2.0 * mu * Z / sqrt(M_PI));
  else
    return - (c * exp(-alpha*alpha*r*r) + erf(mu*r)/r);
}

int main(int argc, char **argv)
{
  if ( !(argc == 2) )
  {
    cerr << "Use: verfgau a" << endl;
    return 1;
  }
  const double mu = atof(argv[1]);
  const double dr = 0.002/Z;
  int np = 501;
  // adjust np so that v(r) = -Z/r at r=(np-)*dr
  while ( fabs(v(mu,(np-1)*dr)+Z/((np-1)*dr)) > 1.e-10 )
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
  cout << "  Verfgau all-electron potential" << endl;
  cout << "  Z= " << Z << " mu=" << mu << " dr=" << dr << endl;
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
    cout << v(mu,r) << endl;
  }
  cout << "</radial_potential>" << endl;
  cout << "</projector>" << endl;
  cout << "</norm_conserving_pseudopotential>" << endl;
  cout << "</fpmd:species>" << endl;
}
