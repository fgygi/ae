////////////////////////////////////////////////////////////////////////////////
//
//  nradsolve.cpp
//  solve the radial Schroedinger equation on a linear radial mesh
//  compute n eigenvectors and eigenvalues
//
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <valarray>
using namespace std;
extern "C"
{
void dstevx_(char*,char*,int*,double*,double*,double*,double*,int*,int*,double*,
          int*,double*,double*,int*,double*,int*,int*,int*);
}

void nradsolve(double rmax, int l, int neig, valarray<double>& v,
  valarray<double>& u, valarray<double>& eig)
{
  // Solve the radial Schroedinger equation on a linear mesh
  // The linear mesh is defined by rmax and n=v.size()
  // h = rmax/(n+1)
  // r[0] = h
  // r[i] = (i+1) * h
  // r[n-1] = rmax - h
  // Wavefunctions u = r * phi are zero at the origin and at rmax
  // Wavefunctions are only represented where they are finite, i.e. from
  // r[0] to r[n-1] inclusive
  // u: valarray of size neig*n
  // eig: valarray of size neig

  // Returns the n lowest eigenfunctions and lowest eigenvalues

  int n = v.size();
  assert(u.size()==neig*n);
  assert(eig.size()==neig);
  const double h = rmax/(n+1);
  const double h2inv = 1.0 / ( h * h );
  // solve symmetric tridiagonal eigenvalue problem
  // diagonal d, subdiagonal e
  valarray<double> d(n),e(n);
  // Note: -0.5*Laplacian + V in next lines
  for ( int i = 0; i < n; i++ )
  {
    const double r = (i+1)*h;
    d[i] = h2inv + v[i] + 0.5*l*(l+1)/(r*r);
  }
  for ( int i = 0; i < n; i++ )
    e[i] = -0.5 * h2inv;

  // LAPACK call: find lowest eigenvalue and eigenvector
  valarray<double> work(5*n);
  valarray<int> iwork(5*n),ifail(n);
  char jobz = 'V';
  char range = 'I';
  double vl,vu;
  int il=1,iu=neig,ldz=n,m;
  double abstol = 1.e-8;
  int info;
  dstevx_(&jobz,&range,&n,&d[0],&e[0],&vl,&vu,&il,&iu,&abstol,&m,&eig[0],
          &u[0],&ldz,&work[0],&iwork[0],&ifail[0],&info);

  // normalize the solutions
  // the diagonalization routine normalizes the solution so that sum(u_i^2)=1
  // Normalization over the interval [0,rmax] implies h*sum(u_i^2) = 1
  // with h = rmax/(np+1);

  // correction: multiply by sqrt((n+1)/rmax) == 1/sqrt(h)
  // all eigenvectors are normalized at the same time
  const double sqrthinv = 1.0 / sqrt(h);
  u *= sqrthinv;;

  assert(info==0);
  return;
}
