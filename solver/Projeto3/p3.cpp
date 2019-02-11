#include <iostream>
#include <cmath>
#include <vector>
#include "p3Header.h"
#include "Mesh.h"
#include "Operators.h"

extern ExactSol exact;

int main() {

  int convCrash = 0;

	// Declaring flow/problem variables
	Matrix<double> dt(mesh.xsize());
	for (size_t i = 0; i != mesh.xsize(); ++i) {
		dt[i].resize(mesh.ysize());
	}

	Matrix<Vector<double> > Q(mesh.xsize()), E_p(mesh.xsize()), E_m(mesh.xsize()),
													F_p(mesh.xsize()), F_m(mesh.xsize()), RES(mesh.xsize());
	for (size_t i = 0; i != mesh.xsize(); ++i) {
		Q[i].resize(mesh.ysize());
		E_p[i].resize(mesh.ysize());
		E_m[i].resize(mesh.ysize());
		F_p[i].resize(mesh.ysize());
		F_m[i].resize(mesh.ysize());
		RES[i].resize(mesh.ysize());
		for (size_t j = 0; j != mesh.ysize(); ++j) {
			Q[i][j].resize(4);
			E_p[i][j].resize(4);
			E_m[i][j].resize(4);
			F_p[i][j].resize(4);
			F_m[i][j].resize(4);
			RES[i][j].resize(4);
		}
	}

	Matrix<Matrix<double> > A_p(mesh.xsize()), A_m(mesh.xsize()), B_p(mesh.xsize()), B_m(mesh.xsize());
	for (size_t i = 0; i != mesh.xsize(); ++i) {
		A_p[i].resize(mesh.ysize());
		A_m[i].resize(mesh.ysize());
		B_p[i].resize(mesh.ysize());
		B_m[i].resize(mesh.ysize());
		for (size_t j = 0; j != mesh.ysize(); ++j) {
			A_p[i][j].resize(4);
			A_m[i][j].resize(4);
			B_p[i][j].resize(4);
			B_m[i][j].resize(4);
			for (size_t k = 0; k != 4; ++k) {
				A_p[i][j][k].resize(4);
				A_m[i][j][k].resize(4);
				B_p[i][j][k].resize(4);
				B_m[i][j][k].resize(4);
			}
		}
	}

  Matrix<Vector<double> > dQstar(mesh.xsize()-2);
  for (auto i = 0; i != mesh.xsize()-2; ++i) {
    dQstar[i].resize(mesh.ysize()-2);
    for (auto j = 0; j != mesh.ysize()-2; ++j) {
      dQstar[i][j].resize(4);
    }
  }

  Matrix<Vector<double> > dQ(mesh.ysize()-2);
  for (auto j = 0; j != mesh.ysize()-2; ++j) {
    dQ[j].resize(mesh.xsize()-2);
    for (auto i = 0; i != mesh.xsize()-2; ++i) {
      dQ[j][i].resize(4);
    }
  }

	// Setting initial flow conditions
	setInitialCond(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m);
	// Setting boundary conditions from initial flow conditions
	updateBoundaryCond(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m);

	writeResult(Q, RES, "results/malha0.vtk");

	// Solution
  switch(exp_imp) {
    case 0:
    	for (int n = 1; n != NMAX+1; ++n) {

    		// Calculate local dt and residuals
    		for (size_t i = 1; i != mesh.xsize()-1; ++i) {
    			for (size_t j = 1; j != mesh.ysize()-1; ++j) {
    				dt[i][j] = (CFL*dx) / (Q[i][j][1]/Q[i][j][0] + sqrt(gama*calcP(Q,i,j)/Q[i][j][0]));

    				calcRES(Q, RES, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m, dt, i, j);
    			}
    		}

    		// Update vector Q
    		for (size_t i = 1; i != mesh.xsize()-1; ++i) {
    			for (size_t j = 1; j != mesh.ysize()-1; ++j) {
    				Q[i][j] = Q[i][j] + RES[i][j];

    				// Update flux vectors
    				calcFluxJ(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m, i, j);
    			}
    		}

    		// Update boundary conditions based on new Q
    		updateBoundaryCond(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m);

        // Check convergence or crash
        if (n != 1) {
          convCrash = checkConvCrash(RES);
        }

    		// Write iteration result to file
    		if (!convCrash && (!(n%write_int) || (n == 1))) {
    			writeResult(Q, RES, "results/malha" + std::to_string(n) + ".vtk");
          writeResidual(RES, n);
    	  } else if (convCrash) {
          if (convCrash > 0) {
            std::cout << "++++++++++ CONVERGIU" << "\n"
                      << "Iteracao: " << n << std::endl;
            writeResult(Q, RES, "results/malha" + std::to_string(n) + ".vtk");
            writeResidual(RES, n);
            int sys = system("gnuplot -p plots/plotres.plt");
            exit(-1);
          } else {
            std::cout << "---------- EXPLODIU" << "\n"
                      << "Iteracao: " << n << std::endl;
            writeResult(Q, RES, "results/malha" + std::to_string(n) + ".vtk");
            writeResidual(RES, n);
            int sys = system("gnuplot -p plots/plotres.plt");
            exit(-1);
          }
        }
      }
    break;

    case 1:
      for (int n = 1; n != NMAX+1; ++n) {

        // Step 1: Residual calc and solution of lower triangular matrix
        for (size_t j = 1; j != mesh.ysize()-1; ++j) {
          for (size_t i = 1; i != mesh.xsize()-1; ++i) {
            dt[i][j] = (CFL*dx) / (Q[i][j][1]/Q[i][j][0] + sqrt(gama*calcP(Q,i,j)/Q[i][j][0]));

            calcRES(Q, RES, E_p, E_m, F_p, F_m, A_m, A_p, B_m, B_p, dt, i, j);
          }
          calcSWLU(dQstar, dQ, RES, A_p, B_p, dt, j, 1);
        }

        // Step 2: Solution of upper triangular matrix (dQ calc)
        for (size_t j = mesh.ysize()-1; j-- > 1 ;) {
          calcSWLU(dQstar, dQ, RES, A_m, B_m, dt, j, 2);
        }

        // Step 3: Update vector Q
        for (size_t i = 1; i != mesh.xsize()-1; ++i) {
          for (size_t j = 1; j != mesh.ysize()-1; ++j) {
            Q[i][j] = Q[i][j] + dQ[j-1][i-1];

            // Update flux vectors
            calcFluxJ(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m, i, j);
          }
        }

        // Update boundary conditions based on new Q
        updateBoundaryCond(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m);

        // Check convergence or crash
        if (n != 1) {
          convCrash = checkConvCrash(RES);
        }

        // Write iteration result to file
        if (!convCrash && (!(n%write_int) || (n == 1))) {
          writeResult(Q, RES, "results/malha" + std::to_string(n) + ".vtk");
          writeResidual(RES, n);
        } else if (convCrash) {
          if (convCrash > 0) {
            std::cout << "++++++++++ CONVERGIU" << "\n"
                      << "Iteracao: " << n << std::endl;
            writeResult(Q, RES, "results/malha" + std::to_string(n) + ".vtk");
            writeResidual(RES, n);
            int sys = system("gnuplot -p plots/plotres.plt");
            exit(-1);
          } else {
            std::cout << "---------- EXPLODIU" << "\n"
                      << "Iteracao: " << n << std::endl;
            writeResult(Q, RES, "results/malha" + std::to_string(n) + ".vtk");
            writeResidual(RES, n);
            int sys = system("gnuplot -p plots/plotres.plt");
            exit(-1);
          }
        }


      }
    break;

    default:
      std::cout << "Invalid time-advancing scheme. Please double check." << std::endl;
      exit(-1);
  }
  std::cout << "+-+-+-+-+- FIM DAS ITERACOES" << "\n"
            << "Nao houve convergencia." << std::endl;
  int sys = system("gnuplot -p plots/plotres.plt");

	return 0;
}