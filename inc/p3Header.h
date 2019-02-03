#ifndef P3HEADER_H
#define P3HEADER_H

#include <string>
#include <vector>
#include <cmath>
#include <initializer_list>
#include "Mesh.h"

template <typename T>
using Vector = std::vector<T>;

template <typename T>
using Matrix = std::vector<std::vector<T> >;

double InitFlow(std::string);
double InitNum(std::string);
std::string InitNumStr(std::string);
std::string readMesh(std::string);
void setInitialCond(Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&,
										Matrix<Matrix<double> >&,Matrix<Matrix<double> >&, Matrix<Matrix<double> >&, Matrix<Matrix<double> >&);
void updateBoundaryCond(Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&,
										Matrix<Matrix<double> >&,Matrix<Matrix<double> >&, Matrix<Matrix<double> >&, Matrix<Matrix<double> >&);
void calcFluxJ(Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&,
										Matrix<Matrix<double> >&,Matrix<Matrix<double> >&, Matrix<Matrix<double> >&, Matrix<Matrix<double> >&, const size_t, const size_t);
void calcRES(Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<double>&, const size_t, const size_t);

void writeResult(Matrix<Vector<double> >&, Matrix<Vector<double> >&, std::string);
double calcLinfRes(Matrix<Vector<double> >&, const size_t);
int checkConvCrash(Matrix<Vector<double> >&);
void writeResidual(Matrix<Vector<double> >&, const int);

const double PI = 3.141592;
const double d2r = PI/180, r2d = 180./PI;

const double gama = InitFlow("Gamma"),
						 M1 = InitFlow("Mach1"),
						 beta1 = InitFlow("Beta1")*d2r,
						 p1 = InitFlow("P1"),
						 T1 = InitFlow("T1"),
						 Rg = InitFlow("Rgas"),
						 V1dir = InitFlow("InflowDir"),
						 rho1 = p1/(Rg*(T1+273)),
						 a1 = sqrt(gama*p1/rho1),
						 V1 = M1*a1,
						 u1 = V1*cos(V1dir*d2r),
						 v1 = V1*sin(V1dir*d2r),
						 p0 = p1*pow((1 + ((gama-1)/2)*pow(M1,2)), (gama/(gama-1))),
						 pref = p0,
						 T0 = (T1+273)*(1 + ((gama-1)/2)*pow(M1,2)),
						 Tref = T0,
						 a0 = sqrt((1 + ((gama-1)/2)*pow(M1,2))*pow(a1,2)),
						 astar = sqrt((2 / (gama+1))*pow(a0,2)),
						 uref = astar,
						 rho_ref = p0/pow(uref,2),
						 Rref = pow(uref,2)/T0,
						 Rnon = Rg/Rref,

						 CFL = InitNum("CFL"),
             Conv = InitNum("ConvergenceCriteria");

const int NMAX = InitNum("IterationMax"),
					write_int = InitNum("WriteInterval"),
          ss_order = (InitNumStr("SteadyStateSpatialOrder") == "first") ? 1 : ((InitNumStr("SteadyStateSpatialOrder") == "second") ? 2 : 0),
          exp_imp = (InitNumStr("Time") == "explicit") ? 0 : ((InitNumStr("Time") == "implicit") ? 1 : 0);

inline double calcP(const Matrix<Vector<double> > &Q, const size_t i, const size_t j) {
	return (gama-1)*(Q[i][j][3] - 0.5*Q[i][j][0]*(pow(Q[i][j][1]/Q[i][j][0],2) + pow(Q[i][j][2]/Q[i][j][0],2)));
}
inline double calcT(const Matrix<Vector<double> > &Q, const size_t i, const size_t j) {
  return calcP(Q,i,j) / (Q[i][j][0]*Rnon);
}
inline double calcEigValPlusMinus(const double Vel, const double c, const int s_s) {
  return 0.5*((Vel + c) + s_s*std::abs(Vel + c));
}

const Mesh mesh(readMesh("MeshFile"));
const double dx = mesh.x(1,0) - mesh.x(0,0),
						 dy = mesh.y(0,1) - mesh.y(0,0);

class ExactSol {
	public:
		ExactSol() { calc_all(); };
		double theta(int conv = 0) { return conv ? Theta*r2d : Theta; };
		double beta2(int conv = 0) { return conv ? Beta2*r2d : Beta2; };
		double betaLinha(int conv = 0) { return conv ? BetaPrime*r2d : BetaPrime; };
		double M2() { return Mach2; };
		double M3() { return Mach3; };
		double rho2() { return Rho2; };
		double rho3() { return Rho3; };
		double T2() { return Temp2; };
		double T3() { return Temp3; };
		double p2() { return P2; };
		double p3() { return P3; };
	private:
		double Theta, Beta2, BetaPrime;
		double Mach2, Mach3;
		double Rho2, Temp2, P2;
		double Rho3, Temp3, P3;
		void calc_all() ;
		double findTheta(double, double);
		double findBeta(double, double);
		double calc_postM(double, double, double);
		double calc_rhoRatio(double, double);
		double calc_pRatio(double, double);
		double calc_TRatio(double, double);
};

#endif