#include <iostream>
#include <cmath>
#include <vector>
#include "p3Header.h"
#include "Mesh.h"
#include "Operators.h"

extern ExactSol exact;

int main() {

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

	std::cout << calcP(Q,1,1) << std::endl;

	// Setting initial flow conditions
	setInitialCond(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m);
	// Setting boundary conditions from initial flow conditions
	updateBoundaryCond(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m);

	writeResult(Q, RES, "results/malha0.vtk");

	// Solution
	for (int n = 1; n != NMAX+1; ++n) {

		// Calculate local dt and residuals
		for (size_t i = 1; i != mesh.xsize()-1; ++i) {
			for (size_t j = 1; j != mesh.ysize()-1; ++j) {
				//dt[i][j] = (CFL*dx) / (Q[i][j][1]/Q[i][j][0] + sqrt(gama*calcP(Q,i,j)/Q[i][j][0]));

				calcRES(RES, Q, E_p, E_m, F_p, F_m, dt, i, j);
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

		// Write iteration result to file
		if (!(n%write_int))
			writeResult(Q, RES, "results/malha" + std::to_string(n) + ".vtk");
	}

	return 0;
}