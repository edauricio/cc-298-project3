#include <fstream>
#include <cmath>
#include <iostream>
#include "p3Header.h"
#include "Operators.h"
#include "json.hpp"
//#define MESHDEBUG
#define NUMDEBUG

using json = nlohmann::json;

ExactSol exact;
double norm_cont, norm_xm, norm_ym, norm_en;

double InitFlow(std::string varName) {
	json Init;
	std::ifstream infile("config");
	Init = json::parse(infile);

	return Init["Initial Conditions"][varName];
}

double InitNum(std::string varName) {
	json Init;
	std::ifstream infile("config");
	Init = json::parse(infile);

	return Init["Numerics"][varName];
}

std::string readMesh(std::string mFile) {
	json MeshData;
	std::ifstream infile("config");
	MeshData = json::parse(infile);

	return MeshData["MeshFile"];
}

std::string InitNumStr(std::string varName) {
  json Init;
  std::ifstream infile("config");
  Init = json::parse(infile);

  return Init["Numerics"][varName];
}

void setInitialCond(Matrix<Vector<double> > &Q, Matrix<Vector<double> > &E_p, Matrix<Vector<double> > &E_m, Matrix<Vector<double> > &F_p, Matrix<Vector<double> > &F_m,
										Matrix<Matrix<double> > &A_p, Matrix<Matrix<double> > &A_m, Matrix<Matrix<double> > &B_p, Matrix<Matrix<double> > &B_m) {
	for (size_t i = 1; i != mesh.xsize()-1; ++i) {
		for (size_t j = 1; j != mesh.ysize()-1; ++j) {
			Q[i][j][0] = rho1/rho_ref;
			Q[i][j][1] = (rho1*u1)/(rho_ref*uref);
			Q[i][j][2] = (rho1*v1)/(rho_ref*uref);
			Q[i][j][3] = (p1/(gama-1) + 0.5*rho1*(pow(u1,2) + pow(v1,2)))/(rho_ref*pow(uref,2));

			calcFluxJ(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m, i, j);
		}
	}
}

void updateBoundaryCond(Matrix<Vector<double> > &Q, Matrix<Vector<double> > &E_p, Matrix<Vector<double> > &E_m, Matrix<Vector<double> > &F_p, Matrix<Vector<double> > &F_m,
										Matrix<Matrix<double> > &A_p, Matrix<Matrix<double> > &A_m, Matrix<Matrix<double> > &B_p, Matrix<Matrix<double> > &B_m) {
	static int cnt = 0;

	// Inflow: Dirichlet. Only set for initial flow; won't change later.
	if (!cnt) {
		for (size_t j = 0; j != mesh.ysize()-1; ++j) {
			Q[0][j][0] = rho1/rho_ref;
			Q[0][j][1] = (rho1*u1)/(rho_ref*uref);
			Q[0][j][2] = (rho1*v1)/(rho_ref*uref);
			Q[0][j][3] = (p1/(gama-1) + 0.5*rho1*(pow(u1,2) + pow(v1,2)))/(rho_ref*pow(uref,2));

			calcFluxJ(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m, 0, j);
		}
	}

	// Top boundary: Dirichlet. Same strategy as inflow.
	if (!cnt) {
		double a2, V2, u2, v2;
		a2 = sqrt(gama*exact.p2()/exact.rho2());
		V2 = exact.M2()*a2;
		u2 = V2*cos(V1dir-exact.theta());
		v2 = V2*sin(V1dir-exact.theta());
		for (size_t i = 0; i != mesh.xsize()-1; ++i) {
			Q[i][mesh.ysize()-1][0] = exact.rho2()/rho_ref;
			Q[i][mesh.ysize()-1][1] = (exact.rho2()*u2)/(rho_ref*uref);
			Q[i][mesh.ysize()-1][2] = (exact.rho2()*v2)/(rho_ref*uref);
			Q[i][mesh.ysize()-1][3] = (exact.p2()/(gama-1) + 0.5*exact.rho2()*(pow(u2,2) + pow(v2,2)))/(rho_ref*pow(uref,2));

			calcFluxJ(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m, i, mesh.ysize()-1);
		}
		++cnt;
	}

	// Wall: Euler BC; dp/dn = dT/dn = 0.
	for (size_t i = 1; i != mesh.xsize()-1; ++i) {
		Q[i][0][0] = calcP(Q,i,1)/(Rnon*calcT(Q,i,1));
		Q[i][0][1] = Q[i][1][1];
		Q[i][0][2] = 0.0;
		Q[i][0][3] = calcP(Q,i,1)/(gama-1) + 0.5*Q[i][0][0]*(pow(Q[i][0][1]/Q[i][0][0],2) + pow(Q[i][0][2]/Q[i][0][0],2));

		calcFluxJ(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m, i, 0);
	}

	// Exit: supersonic
	for (size_t j = 0; j != mesh.ysize(); ++j) {
		Q[mesh.xsize()-1][j][0] = Q[mesh.xsize()-2][j][0];
		Q[mesh.xsize()-1][j][1] = Q[mesh.xsize()-2][j][1];
		Q[mesh.xsize()-1][j][2] = Q[mesh.xsize()-2][j][2];
		Q[mesh.xsize()-1][j][3] = Q[mesh.xsize()-2][j][3];

    calcFluxJ(Q, E_p, E_m, F_p, F_m, A_p, A_m, B_p, B_m, mesh.xsize()-1, j);
	}
}

void calcFluxJ(Matrix<Vector<double> > &Q, Matrix<Vector<double> > &E_p, Matrix<Vector<double> > &E_m, Matrix<Vector<double> > &F_p, Matrix<Vector<double> > &F_m,
										Matrix<Matrix<double> > &A_p, Matrix<Matrix<double> > &A_m, Matrix<Matrix<double> > &B_p, Matrix<Matrix<double> > &B_m, const size_t i, const size_t j) {
	static Matrix<double> P(4), P_i(4), M(4), M_i(4), Diag(4);
	static int cnt = 0;
	double k1, k2, k1t, k2t, c;

	if (!cnt++)
		for (size_t i = 0; i != 4; ++i) {
			P[i].resize(4);
			P_i[i].resize(4);
			M[i].resize(4);
			M_i[i].resize(4);
			Diag[i].resize(4);
	}

	k1 = 1.0;
	k2 = 0;
	k1t = 1.0;
	k2t = 0;
	c = sqrt(gama*calcP(Q,i,j)/Q[i][j][0]);

	M[0][0] = 1.0;
	M[1][0] = Q[i][j][1]/Q[i][j][0];
	M[1][1] = Q[i][j][0];
	M[2][0] = Q[i][j][2]/Q[i][j][0];
	M[2][2] = Q[i][j][0];
	M[3][0] = 0.5*(pow(Q[i][j][1]/Q[i][j][0],2) + pow(Q[i][j][2]/Q[i][j][0], 2));
	M[3][1] = Q[i][j][1];
	M[3][2] = Q[i][j][2];
	M[3][3] = 1./(gama-1);

	M_i[0][0] = 1.0;
	M_i[1][0] = -(Q[i][j][1]/Q[i][j][0])/Q[i][j][0];
	M_i[1][1] = 1./Q[i][j][0];
	M_i[2][0] = -(Q[i][j][2]/Q[i][j][0])/Q[i][j][0];
	M_i[2][2] = 1./Q[i][j][0];
	M_i[3][0] = 0.5*(gama-1)*(pow(Q[i][j][1]/Q[i][j][0],2) + pow(Q[i][j][2]/Q[i][j][0],2));
	M_i[3][1] = -(gama-1)*(Q[i][j][1]/Q[i][j][0]);
	M_i[3][2] = -(gama-1)*(Q[i][j][2]/Q[i][j][0]);
	M_i[3][3] = gama-1;

  P[0][0] = k1t;
  P[0][1] = 0.0;
  P[0][2] = Q[i][j][0]/(sqrt(2)*c);
  P[0][3] = Q[i][j][0]/(sqrt(2)*c);
  P[1][0] = 0.0;
  P[1][1] = k2t;
  P[1][2] = k1t/sqrt(2);
  P[1][3] = -k1t/sqrt(2);
  P[2][0] = 0.0;
  P[2][1] = -k1t;
  P[2][2] = k2t/sqrt(2);
  P[2][3] = -k2t/sqrt(2);
  P[3][0] = 0.0;
  P[3][1] = 0.0;
  P[3][2] = Q[i][j][0]*c/sqrt(2);
  P[3][3] = Q[i][j][0]*c/sqrt(2);

  P_i[0][0] = k1t;
  P_i[0][1] = 0.0;
  P_i[0][2] = 0.0;
  P_i[0][3] = -k1t/pow(c,2);
  P_i[1][0] = 0.0;
  P_i[1][1] = k2t;
  P_i[1][2] = -k1t;
  P_i[1][3] = 0.0;
  P_i[2][0] = 0.0;
  P_i[2][1] = k1t/sqrt(2);
  P_i[2][2] = k2t/sqrt(2);
  P_i[2][3] = 1./(sqrt(2)*Q[i][j][0]*c);
  P_i[3][0] = 0.0;
  P_i[3][1] = -k1t/sqrt(2);
  P_i[3][2] = -k2t/sqrt(2);
  P_i[3][3] = 1./(sqrt(2)*Q[i][j][0]*c);

	Diag[0][0] = calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], 0, 1);
	Diag[1][1] = calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], 0, 1);
	Diag[2][2] = calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], c, 1);
	Diag[3][3] = calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], -1.*c, 1);

	A_p[i][j] = M*(P*Diag*P_i)*M_i;

	Diag[0][0] = calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], 0, -1);
	Diag[1][1] = calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], 0, -1);
	Diag[2][2] = calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], c, -1);
	Diag[3][3] = calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], -1.*c, -1);

  A_m[i][j] = M*(P*Diag*P_i)*M_i;

  E_p[i][j][0] = (Q[i][j][0]/(2*gama))*(2*(gama-1)*calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], 0, 1) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], c, 1) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], -1*c, 1));
  E_p[i][j][1] = (Q[i][j][0]/(2*gama))*(2*(gama-1)*calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], 0, 1)*(Q[i][j][1]/Q[i][j][0]) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], c, 1)*((Q[i][j][1]/Q[i][j][0]) + c*k1t) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], -1*c, 1)*((Q[i][j][1]/Q[i][j][0]) - c*k1t));
  E_p[i][j][2] = (Q[i][j][0]/(2*gama))*(2*(gama-1)*calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], 0, 1)*(Q[i][j][2]/Q[i][j][0]) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], c, 1)*((Q[i][j][2]/Q[i][j][0]) + c*k2t) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], -1*c, 1)*((Q[i][j][2]/Q[i][j][0]) - c*k2t));
  E_p[i][j][3] = (Q[i][j][0]/(2*gama))*((gama-1)*calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], 0, 1)*(pow(Q[i][j][1]/Q[i][j][0],2) + pow(Q[i][j][2]/Q[i][j][0],2)) + 0.5*calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], c, 1)*(pow((Q[i][j][1]/Q[i][j][0] + c*k1t),2) + pow((Q[i][j][2]/Q[i][j][0] + c*k2t),2)) + 0.5*calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], -1*c, 1)*(pow((Q[i][j][1]/Q[i][j][0] - c*k1t),2) + pow((Q[i][j][2]/Q[i][j][0] - c*k2t),2)) + 
                  ((3-gama)*(calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], c, 1) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], -1*c, 1))*pow(c,2))/(2*(gama-1)));


	E_m[i][j][0] = (Q[i][j][0]/(2*gama))*(2*(gama-1)*calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], 0, -1) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], c, -1) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], -1*c, -1));
  E_m[i][j][1] = (Q[i][j][0]/(2*gama))*(2*(gama-1)*calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], 0, -1)*(Q[i][j][1]/Q[i][j][0]) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], c, -1)*((Q[i][j][1]/Q[i][j][0]) + c*k1t) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], -1*c, -1)*((Q[i][j][1]/Q[i][j][0]) - c*k1t));
  E_m[i][j][2] = (Q[i][j][0]/(2*gama))*(2*(gama-1)*calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], 0, -1)*(Q[i][j][2]/Q[i][j][0]) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], c, -1)*((Q[i][j][2]/Q[i][j][0]) + c*k2t) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], -1*c, -1)*((Q[i][j][2]/Q[i][j][0]) - c*k2t));
  E_m[i][j][3] = (Q[i][j][0]/(2*gama))*((gama-1)*calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], 0, -1)*(pow(Q[i][j][1]/Q[i][j][0],2) + pow(Q[i][j][2]/Q[i][j][0],2)) + 0.5*calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], c, -1)*(pow((Q[i][j][1]/Q[i][j][0] + c*k1t),2) + pow((Q[i][j][2]/Q[i][j][0] + c*k2t),2)) + 0.5*calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], -1*c, -1)*(pow((Q[i][j][1]/Q[i][j][0] - c*k1t),2) + pow((Q[i][j][2]/Q[i][j][0] - c*k2t),2)) + 
                  ((3-gama)*(calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], c, -1) + calcEigValPlusMinus(Q[i][j][1]/Q[i][j][0], -1*c, -1))*pow(c,2))/(2*(gama-1)));

	k1 = 0;
	k2 = 1.0;
	k1t = 0;
	k2t = 1.0;

	P[0][0] = k1t;
  P[0][1] = 0.0;
  P[0][2] = Q[i][j][0]/(sqrt(2)*c);
  P[0][3] = Q[i][j][0]/(sqrt(2)*c);
  P[1][0] = 0.0;
  P[1][1] = k2t;
  P[1][2] = k1t/sqrt(2);
  P[1][3] = -k1t/sqrt(2);
  P[2][0] = 0.0;
  P[2][1] = -k1t;
  P[2][2] = k2t/sqrt(2);
  P[2][3] = -k2t/sqrt(2);
  P[3][0] = 0.0;
  P[3][1] = 0.0;
  P[3][2] = Q[i][j][0]*c/sqrt(2);
  P[3][3] = Q[i][j][0]*c/sqrt(2);

  P_i[0][0] = k1t;
  P_i[0][1] = 0.0;
  P_i[0][2] = 0.0;
  P_i[0][3] = -k1t/pow(c,2);
  P_i[1][0] = 0.0;
  P_i[1][1] = k2t;
  P_i[1][2] = -k1t;
  P_i[1][3] = 0.0;
  P_i[2][0] = 0.0;
  P_i[2][1] = k1t/sqrt(2);
  P_i[2][2] = k2t/sqrt(2);
  P_i[2][3] = 1./(sqrt(2)*Q[i][j][0]*c);
  P_i[3][0] = 0.0;
  P_i[3][1] = -k1t/sqrt(2);
  P_i[3][2] = -k2t/sqrt(2);
  P_i[3][3] = 1./(sqrt(2)*Q[i][j][0]*c);

	Diag[0][0] = calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], 0, 1);
	Diag[1][1] = calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], 0, 1);
	Diag[2][2] = calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], c, 1);
	Diag[3][3] = calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], -1.*c, 1);

	B_p[i][j] = M*(P*Diag*P_i)*M_i;

	Diag[0][0] = calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], 0, -1);
	Diag[1][1] = calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], 0, -1);
	Diag[2][2] = calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], c, -1);
	Diag[3][3] = calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], -1.*c, -1);

	B_m[i][j] = M*(P*Diag*P_i)*M_i;

	F_p[i][j][0] = (Q[i][j][0]/(2*gama))*(2*(gama-1)*calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], 0, 1) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], c, 1) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], -1*c, 1));
  F_p[i][j][1] = (Q[i][j][0]/(2*gama))*(2*(gama-1)*calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], 0, 1)*(Q[i][j][1]/Q[i][j][0]) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], c, 1)*((Q[i][j][1]/Q[i][j][0]) + c*k1t) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], -1*c, 1)*((Q[i][j][1]/Q[i][j][0]) - c*k1t));
  F_p[i][j][2] = (Q[i][j][0]/(2*gama))*(2*(gama-1)*calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], 0, 1)*(Q[i][j][2]/Q[i][j][0]) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], c, 1)*((Q[i][j][2]/Q[i][j][0]) + c*k2t) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], -1*c, 1)*((Q[i][j][2]/Q[i][j][0]) - c*k2t));
  F_p[i][j][3] = (Q[i][j][0]/(2*gama))*((gama-1)*calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], 0, 1)*(pow(Q[i][j][1]/Q[i][j][0],2) + pow(Q[i][j][2]/Q[i][j][0],2)) + 0.5*calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], c, 1)*(pow((Q[i][j][1]/Q[i][j][0] + c*k1t),2) + pow((Q[i][j][2]/Q[i][j][0] + c*k2t),2)) + 0.5*calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], -1*c, 1)*(pow((Q[i][j][1]/Q[i][j][0] - c*k1t),2) + pow((Q[i][j][2]/Q[i][j][0] - c*k2t),2)) + 
                  ((3-gama)*(calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], c, 1) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], -1*c, 1))*pow(c,2))/(2*(gama-1)));


  F_m[i][j][0] = (Q[i][j][0]/(2*gama))*(2*(gama-1)*calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], 0, -1) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], c, -1) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], -1*c, -1));
  F_m[i][j][1] = (Q[i][j][0]/(2*gama))*(2*(gama-1)*calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], 0, -1)*(Q[i][j][1]/Q[i][j][0]) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], c, -1)*((Q[i][j][1]/Q[i][j][0]) + c*k1t) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], -1*c, -1)*((Q[i][j][1]/Q[i][j][0]) - c*k1t));
  F_m[i][j][2] = (Q[i][j][0]/(2*gama))*(2*(gama-1)*calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], 0, -1)*(Q[i][j][2]/Q[i][j][0]) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], c, -1)*((Q[i][j][2]/Q[i][j][0]) + c*k2t) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], -1*c, -1)*((Q[i][j][2]/Q[i][j][0]) - c*k2t));
  F_m[i][j][3] = (Q[i][j][0]/(2*gama))*((gama-1)*calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], 0, -1)*(pow(Q[i][j][1]/Q[i][j][0],2) + pow(Q[i][j][2]/Q[i][j][0],2)) + 0.5*calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], c, -1)*(pow((Q[i][j][1]/Q[i][j][0] + c*k1t),2) + pow((Q[i][j][2]/Q[i][j][0] + c*k2t),2)) + 0.5*calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], -1*c, -1)*(pow((Q[i][j][1]/Q[i][j][0] - c*k1t),2) + pow((Q[i][j][2]/Q[i][j][0] - c*k2t),2)) + 
                  ((3-gama)*(calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], c, -1) + calcEigValPlusMinus(Q[i][j][2]/Q[i][j][0], -1*c, -1))*pow(c,2))/(2*(gama-1)));

}

void calcRES(Matrix<Vector<double> > &RES, Matrix<Vector<double> > &E_p, Matrix<Vector<double> > &E_m, Matrix<Vector<double> > &F_p, Matrix<Vector<double> > &F_m, Matrix<double> &dt, const size_t i, const size_t j) {
  switch(ss_order) {
    case 1:
	   RES[i][j] = -dt[i][j]*((1./dx)*(E_p[i][j] - E_p[i-1][j]) + (1./dx)*(E_m[i+1][j] - E_m[i][j]) + (1./dy)*(F_p[i][j] - F_p[i][j-1]) + (1./dy)*(F_m[i][j+1] - F_m[i][j]));
    break;

    case 2:
      if (i == 1) {
        if (j == 1) {
          RES[i][j] = -dt[i][j]*((1./dx)*(E_p[i][j] - E_p[i-1][j]) + (1./(2*dx))*(-3*E_m[i][j] + 4*E_m[i+1][j] - E_m[i+2][j]) + (1./dy)*(F_p[i][j] - F_p[i][j-1]) + (1./(2*dy))*(-3*F_m[i][j] + 4*F_m[i][j+1] - F_m[i][j+2]));
        } else if (j == mesh.ysize()-2) {
          RES[i][j] = -dt[i][j]*((1./dx)*(E_p[i][j] - E_p[i-1][j]) + (1./(2*dx))*(-3*E_m[i][j] + 4*E_m[i+1][j] - E_m[i+2][j]) + (1./(2*dy))*(3*F_p[i][j] - 4*F_p[i][j-1] + F_p[i][j-2]) + (1./dy)*(F_m[i][j+1] - F_m[i][j]));
        } else {
          RES[i][j] = -dt[i][j]*((1./dx)*(E_p[i][j] - E_p[i-1][j]) + (1./(2*dx))*(-3*E_m[i][j] + 4*E_m[i+1][j] - E_m[i+2][j]) + (1./(2*dy))*(3*F_p[i][j] - 4*F_p[i][j-1] + F_p[i][j-2]) + (1./(2*dy))*(-3*F_m[i][j] + 4*F_m[i][j+1] - F_m[i][j+2]));
        }
      } else if (i == mesh.xsize()-2) {
        if (j == 1) {
          RES[i][j] = -dt[i][j]*((1./(2*dx))*(3*E_p[i][j] - 4*E_p[i-1][j] + E_p[i-2][j]) + (1./dx)*(E_m[i+1][j] - E_m[i][j]) + (1./dy)*(F_p[i][j] - F_p[i][j-1]) + (1./(2*dy))*(-3*F_m[i][j] + 4*F_m[i][j+1] - F_m[i][j+2]));
        } else if (j == mesh.ysize()-2) {
          RES[i][j] = -dt[i][j]*((1./(2*dx))*(3*E_p[i][j] - 4*E_p[i-1][j] + E_p[i-2][j]) + (1./dx)*(E_m[i+1][j] - E_m[i][j]) + (1./(2*dy))*(3*F_p[i][j] - 4*F_p[i][j-1] + F_p[i][j-2]) + (1./dy)*(F_m[i][j+1] - F_m[i][j]));
        } else {
          RES[i][j] = -dt[i][j]*((1./(2*dx))*(3*E_p[i][j] - 4*E_p[i-1][j] + E_p[i-2][j]) + (1./dx)*(E_m[i+1][j] - E_m[i][j]) + (1./(2*dy))*(3*F_p[i][j] - 4*F_p[i][j-1] + F_p[i][j-2]) + (1./(2*dy))*(-3*F_m[i][j] + 4*F_m[i][j+1] - F_m[i][j+2]));
        }
      } else {
        if (j == 1) {
          RES[i][j] = -dt[i][j]*((1./(2*dx))*(3*E_p[i][j] - 4*E_p[i-1][j] + E_p[i-2][j]) + (1./(2*dx))*(-3*E_m[i][j] + 4*E_m[i+1][j] - E_m[i+2][j]) + (1./dy)*(F_p[i][j] - F_p[i][j-1]) + (1./(2*dy))*(-3*F_m[i][j] + 4*F_m[i][j+1] - F_m[i][j+2]));
        } else if (j == mesh.ysize()-2) {
          RES[i][j] = -dt[i][j]*((1./(2*dx))*(3*E_p[i][j] - 4*E_p[i-1][j] + E_p[i-2][j]) + (1./(2*dx))*(-3*E_m[i][j] + 4*E_m[i+1][j] - E_m[i+2][j]) + (1./(2*dy))*(3*F_p[i][j] - 4*F_p[i][j-1] + F_p[i][j-2]) + (1./dy)*(F_m[i][j+1] - F_m[i][j]));
        } else {
          RES[i][j] = -dt[i][j]*((1./(2*dx))*(3*E_p[i][j] - 4*E_p[i-1][j] + E_p[i-2][j]) + (1./(2*dx))*(-3*E_m[i][j] + 4*E_m[i+1][j] - E_m[i+2][j]) + (1./(2*dy))*(3*F_p[i][j] - 4*F_p[i][j-1] + F_p[i][j-2]) + (1./(2*dy))*(-3*F_m[i][j] + 4*F_m[i][j+1] - F_m[i][j+2]));
        }
      }
    break;

    default:
      std::cout << "Invalid steady-state order. Please double check." << "\n";
      exit(-1);
    break;
  }
}

double ExactSol::findTheta(double M, double beta) {
	static double theta = 0.1;
	double num = (gama-1)*pow(M,2)*pow(sin(beta),2) + 2;
	double den = (gama+1)*pow(M,2)*pow(sin(beta),2);
	double lhs = tan(beta-theta);
	double rhs = num*tan(beta)/den;
	if (std::abs(lhs-rhs) > 0.00001) {
    if ((lhs-rhs) > 0) {
      theta = theta + (theta/100);
    } else {
      theta = theta - (theta/100);
    }
    findTheta(M, beta);
  }
	return theta;	
}

double ExactSol::findBeta(double M, double theta) {
	static double beta = 0.1;
	double num = (gama-1)*pow(M,2)*pow(sin(beta),2) + 2;
	double den = (gama+1)*pow(M,2)*pow(sin(beta),2);
	double lhs = tan(beta-theta)/tan(beta);
	double rhs = num/den;
	if (std::abs(lhs - rhs) > 0.00001) {
    if ((lhs-rhs) > 0) {
      beta = beta - (beta/100);
    } else {
      beta = beta + (beta/100);
    }
    findBeta(M, theta);
  }
	return beta;	
}

double ExactSol::calc_postM(double M, double beta, double theta) {
	double num = 1 + ((gama-1)/2)*pow(M,2)*pow(sin(beta),2);
	double den = (gama*pow(M,2)*pow(sin(beta),2) - (gama-1)/2)*pow(sin(beta-theta),2);
	return sqrt(num/den);
}

double ExactSol::calc_rhoRatio(double M, double beta) {
	double num = (gama+1)*pow(M,2)*pow(sin(beta),2);
	double den = (gama-1)*pow(M,2)*pow(sin(beta),2) + 2;
	return num/den;
}

double ExactSol::calc_pRatio(double M, double beta) {
	double mult = (2*gama)/(gama+1);
	double par = pow(M,2)*pow(sin(beta),2) - 1;
	return mult*par + 1;
}

double ExactSol::calc_TRatio(double M, double beta) {
	double mult1 = (2*(gama-1))/pow(gama+1,2);
	double mult2_num = pow(M,2)*pow(sin(beta),2) - 1;
	double mult2_den = pow(M,2)*pow(sin(beta),2);
	double mult2 = mult2_num/mult2_den;
	double par = gama*pow(M,2)*pow(sin(beta),2) + 1;
	return 1.0 + mult1 * mult2 * par;
}

void ExactSol::calc_all() {
	Theta = findTheta(M1, beta1);
	Mach2 = calc_postM(M1, beta1, Theta);
	P2 = p1*calc_pRatio(M1, beta1);
	Rho2 = rho1*calc_rhoRatio(M1, beta1);
	Temp2 = T1*calc_TRatio(M1, beta1);
	Beta2 = findBeta(Mach2, Theta);
	BetaPrime = Beta2 - Theta;
	Mach3 = calc_postM(Mach2, Beta2, Theta);
	Rho3 = Rho2*calc_rhoRatio(Mach2, Beta2);
	P3 = P2*calc_pRatio(Mach2, Beta2);
	Temp3 = Temp2*calc_TRatio(Mach2, Beta2);
}

void writeResult(Matrix<Vector<double> > &Q, Matrix<Vector<double> > &RES, std::string fName) {
	ofstream wFile;
  wFile.open(fName);
  wFile.precision(6);
  wFile << scientific;
  wFile << "# vtk DataFile Version 3.0" << "\n"
        << "Resultado Projeto 3" << "\n"
        << "ASCII" << "\n"
        << "DATASET STRUCTURED_GRID" << "\n"
        << "DIMENSIONS " << mesh.xsize() << " " << mesh.ysize() << " " << mesh.zsize() << " " << "\n"
        << "POINTS " << mesh.xsize()*mesh.ysize()*mesh.zsize() << " double" << "\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
    for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
      for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
        wFile << mesh.x(i,j,k) << " " << mesh.y(i,j,k) << " " << mesh.z(i,j,k) << "\n";
      }
    }
  }
  wFile << "POINT_DATA " << mesh.xsize()*mesh.ysize()*mesh.zsize() << "\n";
#ifdef MESHDEBUG  
  wFile << "SCALARS Jacob double" << "\n"
  << "LOOKUP_TABLE default" << "\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << mesh.J(i,j,k) << "\n";
  		}
  	}
  }
  wFile << "SCALARS xi_x double" << "\n"
  << "LOOKUP_TABLE default" << "\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << mesh.xi_x(i,j,k) << "\n";
  		}
  	}
  }
  wFile << "SCALARS xi_y double" << "\n"
  << "LOOKUP_TABLE default" << "\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << mesh.xi_y(i,j,k) << "\n";
  		}
  	}
  }
  wFile << "SCALARS eta_x double" << "\n"
  << "LOOKUP_TABLE default" << "\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << mesh.eta_x(i,j,k) << "\n";
  		}
  	}
  }
  wFile << "SCALARS eta_y double" << "\n"
  << "LOOKUP_TABLE default" << "\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << mesh.eta_y(i,j,k) << "\n";
  		}
  	}
  }
#endif
#ifdef NUMDEBUG
  wFile << "SCALARS RES_Cont double" << "\n"
  << "LOOKUP_TABLE default" << "\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << RES[i][j][0] << "\n";
  		}
  	}
  }
  wFile << "SCALARS RES_X-Momentum double" << "\n"
  << "LOOKUP_TABLE default" << "\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << RES[i][j][1] << "\n";
  		}
  	}
  }
  wFile << "SCALARS RES_Y-Momentum double" << "\n"
  << "LOOKUP_TABLE default" << "\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << RES[i][j][2] << "\n";
  		}
  	}
  }
  wFile << "SCALARS RES_Energy double" << "\n"
  << "LOOKUP_TABLE default" << "\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << RES[i][j][3] << "\n";
  		}
  	}
  }
#endif

  wFile << "SCALARS Rho double\n"
  			<< "LOOKUP_TABLE default\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << Q[i][j][0] << "\n";
  		}
  	}
  }
  wFile << "SCALARS Energy double\n"
  			<< "LOOKUP_TABLE default\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << Q[i][j][3] << "\n";
  		}
  	}
  }
  wFile << "SCALARS Pressure double\n"
  			<< "LOOKUP_TABLE default\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << calcP(Q,i,j) << "\n";
  		}
  	}
  }
  wFile << "SCALARS Temperature double\n"
  			<< "LOOKUP_TABLE default\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << calcT(Q,i,j) << "\n";
  		}
  	}
  }
  wFile << "VECTORS Velocity double\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << Q[i][j][1]/Q[i][j][0] << " " << Q[i][j][2]/Q[i][j][0] << " 0" << "\n";
  		}
  	}
  }
  wFile << "VECTORS Mach double\n";
  for (Mesh::mesh_index k = 0; k != mesh.zsize(); ++k) {
  	for (Mesh::mesh_index j = 0; j != mesh.ysize(); ++j) {
  		for (Mesh::mesh_index i = 0; i != mesh.xsize(); ++i) {
  			wFile << (Q[i][j][1]/Q[i][j][0])/sqrt(gama*calcP(Q,i,j)/Q[i][j][0]) << " " << (Q[i][j][2]/Q[i][j][0])/sqrt(gama*calcP(Q,i,j)/Q[i][j][0]) << " 0" << "\n";
  		}
  	}
  }
  wFile << std::flush;
  wFile.close();
}

double calcLinfRes(Matrix<Vector<double> > &RES, const size_t comp) {
  double max_res = RES[1][1][comp];
  for (size_t i = 1; i != mesh.xsize()-1; ++i) {
    for (size_t j = 1; j != mesh.ysize()-1; ++j) {
      if (std::abs(RES[i][j][comp]) > max_res)
        max_res = std::abs(RES[i][j][comp]);
    }
  }
  return max_res;
}

int checkConvCrash(Matrix<Vector<double> > &RES) {
  if (std::isnan(calcLinfRes(RES,0)) || std::isnan(calcLinfRes(RES,1)) || std::isnan(calcLinfRes(RES,2)) || std::isnan(calcLinfRes(RES,3)) ||
      std::isinf(calcLinfRes(RES,0)) || std::isinf(calcLinfRes(RES,1)) || std::isinf(calcLinfRes(RES,2)) || std::isinf(calcLinfRes(RES,3)))
    return -1;
  else if (max(max(max(calcLinfRes(RES,0)/norm_cont, calcLinfRes(RES,1)/norm_xm), calcLinfRes(RES,2)/norm_ym), calcLinfRes(RES,3)/norm_en) < pow(10, -Conv)) {
    return 1;
  }

  return 0;
}

void writeResidual(Matrix<Vector<double> > &RES, const int n) {
  static int cnt = 0;
  std::ofstream wFile;
  if (!cnt++) {
    norm_cont = calcLinfRes(RES,0);
    norm_xm = calcLinfRes(RES,1);
    norm_ym = calcLinfRes(RES,2);
    norm_en = calcLinfRes(RES,3);
    wFile.open("plots/RES.data");
    if (!wFile.is_open()) {
      std::cout << "Error opening residual file" << std::endl;
      exit(-1);
    }
    wFile.precision(6);
    wFile << "#n\tContinuity\tX-Momentum\tY-Momentum\tEnergy\n"
          << n << "\t\t" << log10(calcLinfRes(RES, 0)/norm_cont) << "\t\t" << log10(calcLinfRes(RES, 1)/norm_xm) << "\t\t" << log10(calcLinfRes(RES, 2)/norm_ym) << "\t\t" << log10(calcLinfRes(RES, 3)/norm_en) << "\n";
  } else {
    wFile.open("plots/RES.data", std::ofstream::app);
    if (!wFile.is_open()) {
      std::cout << "Error opening residual file" << std::endl;
      exit(-1);
    }
    wFile << n << "\t\t" << log10(calcLinfRes(RES, 0)/norm_cont) << "\t\t" << log10(calcLinfRes(RES, 1)/norm_xm) << "\t\t" << log10(calcLinfRes(RES, 2)/norm_ym) << "\t\t" << log10(calcLinfRes(RES, 3)/norm_en) << "\n";
  }
  wFile.close();
}