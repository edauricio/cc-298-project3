#ifndef P3HEADER_H
#define P3HEADER_H

#include <fstream>
#include "json.hpp"

using json = nlohmann::json;

double Init(std::string);

const double PI = 3.141592;
const double d2r = PI/180, r2d = 180./PI;

const double gama = Init("Gamma"),
						 M1 = Init("Mach1"),
						 beta1 = Init("Beta1")*d2r,
						 rho1 = Init("Rho1"),
						 p1 = Init("P1"),
						 T1 = Init("T1");

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