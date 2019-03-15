#ifndef P3HEADER_H
#define P3HEADER_H

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <initializer_list>
#include "Mesh.h"
#include "LinearSysSolvers.h"

template <typename P>
using Vector = std::vector<P>;

template <typename P>
using Matrix = std::vector<std::vector<P> >;

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
std::vector<double> calcNumFlux(Matrix<Vector<double> >&, const size_t, const size_t, const int, const int);
void calcAD(Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<double>&, const size_t, const size_t, const int);
void calcRES(Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Matrix<double> >&, Matrix<Matrix<double> >&, Matrix<Matrix<double> >&, Matrix<Matrix<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<double>&, const size_t, const size_t);

void writeResult(Matrix<Vector<double> >&, Matrix<Vector<double> >&, std::string);
double calcLinfRes(Matrix<Vector<double> >&, const size_t);
int checkConvCrash(Matrix<Vector<double> >&);
void writeResidual(Matrix<Vector<double> >&, const int);
void calcSWLU(Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Vector<double> >&, Matrix<Matrix<double> >&, Matrix<Matrix<double> >&, Matrix<double>&, const size_t&, const int&);

const double PI = 3.141592;
const double d2r = PI/180, r2d = 180./PI;

const double gama = InitFlow("Gamma"),
						 M1 = InitFlow("Mach1"),
						 beta1 = InitFlow("Beta1")*d2r,
						 p1 = InitFlow("P1"),
						 T1 = InitFlow("T1")+273,
						 Rg = InitFlow("Rgas"),
						 V1dir = InitFlow("InflowDir"),
						 rho1 = p1/(Rg*T1),
						 a1 = sqrt(gama*p1/rho1),
						 V1 = M1*a1,
						 u1 = V1*cos(V1dir*d2r),
						 v1 = V1*sin(V1dir*d2r),
						 p0 = p1*pow((1 + ((gama-1)/2)*pow(M1,2)), (gama/(gama-1))),
						 pref = p0,
						 T0 = T1*(1 + ((gama-1)/2)*pow(M1,2)),
						 Tref = T0,
						 a0 = sqrt((1 + ((gama-1)/2)*pow(M1,2))*pow(a1,2)),
						 astar = sqrt((2 / (gama+1))*pow(a0,2)),
						 uref = astar,
						 rho_ref = p0/pow(uref,2),
						 Rref = pow(uref,2)/Tref,
						 Rnon = Rg/Rref,

						 CFL = InitNum("CFL"),
             Conv = InitNum("ConvergenceCriteria"),
             alphaL = InitNum("AlphaLiou"),
             betaL = InitNum("BetaLiou"),
             eps_e = InitNum("ConstantEpsilon"),
             kappa2 = InitNum("Kappa2"),
             kappa4 = InitNum("Kappa4");

const int NMAX = InitNum("IterationMax"),
					write_int = InitNum("WriteInterval"),
          method = (InitNumStr("Scheme") == "StegerWarming") ? 1 : ((InitNumStr("Scheme") == "vanLeer") ? 2 : ((InitNumStr("Scheme") == "Liou") ? 3 : ((InitNumStr("Scheme") == "BeamWarming") ? 4 : 0))),
          ss_order = (InitNumStr("SteadyStateSpatialOrder") == "first") ? 1 : ((InitNumStr("SteadyStateSpatialOrder") == "second") ? 2 : 0),
          exp_imp = (InitNumStr("Time") == "explicit") ? 0 : ((InitNumStr("Time") == "implicit") ? 1 : 0),
          imp_treat = (InitNumStr("TransientSpatialOrder") == "first") ? 1 : ((InitNumStr("TransientSpatialOrder") == "second") ? 2 : 0),
          ad_type = (InitNumStr("ArtificialDissipationType") == "Linear") ? 1 : ((InitNumStr("ArtificialDissipationType") == "Isotropic") ? 2 : ((InitNumStr("ArtificialDissipationType") == "Anisotropic") ? 3 : 0));
template <typename T> inline int sgnf(T val) {
  return (val > 0) - (val < 0);
}
inline double calcP(const Matrix<Vector<double> > &Q, const size_t i, const size_t j) {
	return (gama-1)*(Q[i][j][3] - 0.5*Q[i][j][0]*(pow(Q[i][j][1]/Q[i][j][0],2) + pow(Q[i][j][2]/Q[i][j][0],2)));
}
inline double calcT(const Matrix<Vector<double> > &Q, const size_t i, const size_t j) {
  return calcP(Q,i,j) / (Q[i][j][0]*Rnon);
}
inline double calcEigValPlusMinus(const double Vel, const double c, const int s_s) {
  return 0.5*((Vel + c) + s_s*std::abs(Vel + c));
}
inline double calcaStarL(Matrix<Vector<double> > &Q, const size_t i, const size_t j) {
  return sqrt((2*(Q[i][j][3] + calcP(Q,i,j))*(gama-1))/(Q[i][j][0]*(gama+1)));
}
inline double calcaTildeE(Matrix<Vector<double> > &Q, const size_t i, const size_t j) {
  return pow(calcaStarL(Q,i,j), 2)/std::max(calcaStarL(Q,i,j), std::abs(Q[i][j][1]/Q[i][j][0]));
}
inline double calcaTildeF(Matrix<Vector<double> > &Q, const size_t i, const size_t j) {
  return pow(calcaStarL(Q,i,j), 2)/std::max(calcaStarL(Q,i,j), std::abs(Q[i][j][2]/Q[i][j][0]));
}
inline double calcMStyle(const double Ma, const int sign) {
  return (std::abs(Ma) >= 1.0) ? (0.5*(Ma + sign*std::abs(Ma))) : (sign*0.25*pow(Ma + sign, 2) + sign*betaL*pow(pow(Ma,2) - 1, 2));
}
inline double calcmFaceE(Matrix<Vector<double> > &Q, const double aface, const size_t i, const size_t j) {
  return calcMStyle((Q[i][j][1]/Q[i][j][0])/aface, 1) + calcMStyle((Q[i+1][j][1]/Q[i+1][j][0])/aface, -1);
}
inline double calcmFaceF(Matrix<Vector<double> > &Q, const double aface, const size_t i, const size_t j) {
  return calcMStyle((Q[i][j][2]/Q[i][j][0])/aface, 1) + calcMStyle((Q[i][j+1][2]/Q[i][j+1][0])/aface, -1);
}
inline double calcPStyle(const double Ma, const int sign) {
  return (std::abs(Ma) >= 1.0) ? (0.5*(1 + sign*sgnf(Ma))) : (0.25*pow(Ma + sign, 2)*(2 - sign*Ma) + sign*alphaL*Ma*pow(pow(Ma,2) - 1, 2));
}
inline double calcSpeedSound(Matrix<Vector<double> > &Q, const size_t i, const size_t j) {
	return sqrt(gama*calcP(Q,i,j)/Q[i][j][0]);
}
inline double calcSigmaIso(Matrix<Vector<double> > &Q, const size_t i, const size_t j) {
  return std::abs(Q[i][j][1]/Q[i][j][0]) + std::abs(Q[i][j][2]/Q[i][j][0]) + 2*calcSpeedSound(Q, i, j);
}
inline double calcSigmaAniso(Matrix<Vector<double> > &Q, const size_t i, const size_t j, const size_t eq) {
  return std::abs(Q[i][j][eq]/Q[i][j][0]) + calcSpeedSound(Q, i, j);
}
inline double calcSIGMA(Matrix<Vector<double> > &Q, const size_t i, const size_t j, const size_t disp_i, const size_t disp_j) {
  return std::abs(calcP(Q, i+disp_i, j+disp_j) - 2*calcP(Q, i, j) + calcP(Q, i-disp_i, j-disp_j))/std::abs(calcP(Q, i+disp_i, j+disp_j) + 2*calcP(Q, i, j) + calcP(Q, i-disp_i, j-disp_j));
}
inline double calcEps2(Matrix<Vector<double> > &Q, Matrix<double> &dt, const size_t i, const size_t j, const size_t disp_i, const size_t disp_j) {
  return kappa2*dt[i][j]*std::max({calcSIGMA(Q, i+disp_i, j+disp_j, disp_i, disp_j), calcSIGMA(Q, i, j, disp_i, disp_j), calcSIGMA(Q, i-disp_i, j-disp_j, disp_i, disp_j)});
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