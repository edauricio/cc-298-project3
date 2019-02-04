#include <vector>
#include "LinearSysSolvers.h"
#include "Operators.h"

/* TRIDIAGONAL LINEAR SYSTEM SOLVER FUNCTION
 * Description:
 * Solves linear system of equations of the form Ax = f of size n, for a tridiagonal coefficient matrix, that is, 
 * a matrix whose coefficients are zero for all |i - j| > 1.
 * Can handle systems of any size.
 * 
 * Usage:
 * Arguments are matrix A and vector f. Return solution vector x.
 * Of course, A must be a square matrix, and f dimensions must be consistent with A (i.e. n).
 * If matrix A is not tridiagonal, it solves the system using Gauss(...) method instead.
 * */
template <typename T>
std::vector<T> Trid(const std::vector<std::vector<T> > &A, const std::vector<T> &f) {
  const size_t n = f.size();
  if (A[0].size() != 3) {
    std::cout << "The system is not tridiagonal or coefficient matrix is not in compact form. Exiting. . ." << std::endl;
    exit(-1);
  } else if (A.size() != f.size()) {
    std::cout << "Size of matrix and vector are not compatible. Exiting. . ." << std::endl;
    exit(-1);
  }

  std::vector<T> alfa(n), gama(n), g(n), x(n);

  //Constructing diagonals alfa and gama, and solving for the intermediate solution Lg = f
  alfa[0] = A[0][1];
  gama[0] = A[0][2]/alfa[0];
  g[0] = f[0]/alfa[0];
  for (size_t i = 1; i != n-1; i++) {
    alfa[i] = A[i][1] - A[i][0]*gama[i-1];
    gama[i] = A[i][2]/alfa[i];
    g[i] = (f[i] - A[i][0]*g[i-1])/alfa[i];
  }
  alfa[n-1] = A[n-1][1] - A[n-1][0]*gama[n-2];
  g[n-1] = (f[n-1] - A[n-1][0]*g[n-2])/alfa[n-1];
  
  //Solving for Ux = g
  x[n-1] = g[n-1];
  for (size_t i = n-1; i-- > 0 ;) {
    x[i] = g[i] - gama[i]*x[i+1];
  }

  return x;   
}

/* PENTADIAGONAL LINEAR SYSTEM SOLVER FUNCTION
 * Description:
 * Solves linear system of equations of the form Ax = f, of size n, where A is a pentadiagonal
 * banded matrix.
 * Can handle systems of any size.

 * Usage:
 * Argument matrix A must be in 'compact form', where only the nonzero values are stored.
 * Thus, A is a n-by-5 matrix, on which the 5 columns are the subdiagonals, diagonal and superdiagonals
 * nonzero values.
 * No restrictions on vector f.
 * Returns solution vector x.
 *
 * Reference: 
 * Askar, S. S., Karawia, A. A., "On Solving Pentadiagonal Linear Systems via Transformations"
 * Mathematical problems in engineering, volume 2015, Article ID 232456.
 * http://dx.doi.org/10.1155/2015/232456
 */
template <typename T>
std::vector<T> Penta(const std::vector<std::vector<T> > &A, const std::vector<T> &f) {
	if (A[0].size() != 5) {
		cout << "The system is not pentadiagonal or coefficient matrix is in invalid form. Exiting. . ." << endl;
		exit(-1);
	} else if (A.size() != f.size()) {
		cout << "Dimensions of matrix A and vector f don't match. Exiting. . ." << endl;
		exit(-1);
	}
	const size_t n = f.size();
	std::vector<T> alfa(n), beta(n), gama(n), mi(n), z(n);
	std::vector<T> x(n);

	// P/ i = 1;
	mi[0] = A[0][2];
	alfa[0] = A[0][3]/mi[0];
	beta[0] = A[0][4]/mi[0];
	z[0] = f[0]/mi[0];

	// P/ i = 2;
	gama[1] = A[1][1];
	mi[1] = A[1][2] - alfa[0]*gama[1];
	alfa[1] = (A[1][3] - beta[0]*gama[1])/mi[1];
	beta[1] = A[1][4]/mi[1];
	z[1] = (f[1] - z[0]*gama[1])/mi[1];

	// P/ i = 3 ate (n-2)
	for (size_t i = 2; i != (n-2); ++i) {
		gama[i] = A[i][1] - alfa[i-2]*A[i][0];
		mi[i] = A[i][2] - beta[i-2]*A[i][0] - alfa[i-1]*gama[i];
		alfa[i] = (A[i][3] - beta[i-1]*gama[i])/mi[i];
		beta[i] = A[i][4]/mi[i];
		z[i] = (f[i] - z[i-2]*A[i][0] - z[i-1]*gama[i])/mi[i];
	}
	gama[n-2] = A[n-2][1] - alfa[n-4]*A[n-2][0];
	mi[n-2] = A[n-2][2] - beta[n-4]*A[n-2][0] - alfa[n-3]*gama[n-2];
	alfa[n-2] = (A[n-2][3] - beta[n-3]*gama[n-2])/mi[n-2];
	z[n-2] = (f[n-2] - z[n-4]*A[n-2][0] - z[n-3]*gama[n-2])/mi[n-2];

	gama[n-1] = A[n-1][1] - alfa[n-3]*A[n-1][0];
	mi[n-1] = A[n-1][2] - beta[n-3]*A[n-1][0] - alfa[n-2]*gama[n-1];
	z[n-1] = (f[n-1] - z[n-3]*A[n-1][0] - z[n-2]*gama[n-1])/mi[n-1];

	// Backward substitution
	x[n-1] = z[n-1];
	x[n-2] = z[n-2] - alfa[n-2]*x[n-1];
	for (size_t i = n-2; i-- > 0;) {
		x[i] = z[i] - alfa[i]*x[i+1] - beta[i]*x[i+2];
	}

	return x;
}

/* GAUSS ELIMINATION LINEAR SYSTEM SOLVER FUNCTION
 * Description:
 * Solves linear system of equations of the form Ax = f of size n.
 * Can handle systems of any size.
 * 
 * Usage:
 * Arguments are matrix A and vector f. Return solution vector x.
 * Of course, A must be a square matrix, and f dimensions must be consistent with A (i.e. n).
 * */
template <typename T>
std::vector<T> Gauss(const std::vector<std::vector<T> > &A, const std::vector<T> &f) {
  if (A[0].size() != A.size()) {
    std::cout << "Matrix is not square. Exiting. . ." << std::endl;
    exit(-1);
  } else if (A[0].size() != f.size()) {
    std::cout << "Size of matrix and vector are not compatible. Exiting. . ." << std::endl;
    exit(-1);
  }
  const size_t n = f.size();

  std::vector<std::vector<std::vector<T> > > U(n);
  std::vector<std::vector<T> > g(n);
  std::vector<T> x(n);
  
  // Initialization
  for (size_t k = 0; k != n; k++) {
    U[k].resize(n);
    g[k].resize(n);
    for (size_t p = 0; p != n; ++p)
    	U[k][p].resize(n);
  }
  
  g[0] = f;
  U[0] = A;

  for (size_t k = 1; k != n; k++) {
      for (size_t i = 0; i != n; i++) {
          for (size_t j = 0; j != n; j++) {
              if (i >= k) {
              	g[k][i] = g[k-1][i] - (U[k-1][i][k-1]/U[k-1][k-1][k-1])*g[k-1][k-1];
              	if (j <= (k-1)) {
              		U[k][i][j] = 0;
              	} else if (j >= k) {
              		U[k][i][j] = U[k-1][i][j] - (U[k-1][i][k-1]/U[k-1][k-1][k-1])*U[k-1][k-1][j];
              	}
              } else if (i <= (k-1)) {
              	g[k][i] = g[k-1][i];
              	U[k][i][j] = U[k-1][i][j];
              }
          }
      }
  }

  double sum = 0;
  x[n-1] = g[n-1][n-1]/U[n-1][n-1][n-1];
  for (size_t i = n-1; i-- > 0 ;) {
      sum = 0;
      for (size_t j = i+1; j < n; j++) {
          sum = sum + U[n-1][i][j]*x[j];
      }
      x[i] = (g[n-1][i] - sum)/U[n-1][i][i];
  }

  return x;
}

/* BLOCK-TRIDIAGONAL SYSTEM SOLVER
 * Description:
 * Solves block-matrices systems of the form Mx = f of size n, where M is a tridiagonal block-matrix of size n,
 * on which each element is a square matrix of size m, and f is a vector of size n, on which each element of f 
 * is a vector of size m.
 * Can handle systems of any size.
 * 
 * Usage:
 * Arguments are the n elements of main diagonal, A, the n elements of the diagonal below main diagonal, B,
 * the n elements of the diagonal above main diagonal, C, and size of system n.
 * Returns solution vector x, where each element of x is a vector of size m.
 * */
template <typename T>
std::vector<std::vector<T> > Block(const std::vector<std::vector<std::vector<T> > > &A, const std::vector<std::vector<std::vector<T> > > &B, const std::vector<std::vector<std::vector<T> > > &C, const std::vector<std::vector<T> > &f) {
  const size_t n = A.size();
  const size_t m = A[0].size();

  std::vector<std::vector<std::vector<T> > > AL(n), gama(n);
  std::vector<std::vector<T> > x(n), y(n);
  std::vector<T> auxc(m);
  
  //Initializing matrices Aem and gama, and soluion vector x
  for (auto k = 0; k != n; ++k) {
    AL[k].resize(m);
    gama[k].resize(m);
    x[k].resize(m);
    y[k].resize(m);
    for (auto j = 0; j != m; ++j) {
    	AL[k][j].resize(m);
    	gama[k][j].resize(m);
    }
  }

  AL[0] = A[0];
  for (size_t j = 0; j != m; ++j) {
    auxc = Gauss(A[0], C[0]<<j);
    for (size_t p = 0; p != m; ++p) {
      gama[0][p][j] = auxc[p];
    }
  }

  for (size_t k = 1; k != n-1; k++) {
    AL[k] = A[k] - B[k]*gama[k-1];
    for (size_t j = 0; j < m; j++) {
      auxc = Gauss(AL[k], C[k]<<j);
      for (size_t p = 0; p < m; p++) {
        gama[k][p][j] = auxc[p];
      }
    }
  }
  AL[n-1] = A[n-1] - B[n-1]*gama[n-2];

  y[0] = Gauss(AL[0], f[0]);
  for (int k = 1; k < n; k++) {
    y[k] = Gauss(AL[k], (f[k] - B[k]*y[k-1]));
  }

  x[n-1] = y[n-1];
  for (int k = n-2; k >= 0; k--) {
    x[k] = y[k] - gama[k]*x[k+1];
  }

  return x;
}