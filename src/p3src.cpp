#include <cmath>
#include <fstream>
#include <string>
#include "p3Header.h"
#include "json.hpp"

using json = nlohmann::json;

double Init(std::string varName) {
	json Init;
	std::ifstream infile("flow");
	Init = json::parse(infile);

	return Init["Initial Conditions"][varName];
}

// double Init(double *Mach, double *Gama, double *Beta, double *P0, double *T0, double *Rho0,
// 						double *P, double *Rho, double *T) {
// 	json Init;
// 	std::ifstream infile("../solver/flow");
// 	Init = json::parse(infile);
// }

double ExactSol::findTheta(double M, double beta) {
	static double theta = 0.1;
	double num = (gama-1)*pow(M,2)*pow(sin(beta),2) + 2;
	double den = (gama+1)*pow(M,2)*pow(sin(beta),2);
	double lhs = tan(beta-theta);
	double rhs = num*tan(beta)/den;
	if (std::abs(lhs-rhs) < 0.00001) {
		return theta;
	} else {
		if ((lhs-rhs) > 0) {
			theta = theta + (theta/100);
		} else {
			theta = theta - (theta/100);
		}
		findTheta(M, beta);
	}
}

double ExactSol::findBeta(double M, double theta) {
	static double beta = 0.1;
	double num = (gama-1)*pow(M,2)*pow(sin(beta),2) + 2;
	double den = (gama+1)*pow(M,2)*pow(sin(beta),2);
	double lhs = tan(beta-theta)/tan(beta);
	double rhs = num/den;
	if (std::abs(lhs - rhs) < 0.00001) {
		return beta;
	} else {
		if ((lhs-rhs) > 0) {
			beta = beta - (beta/100);
		} else {
			beta = beta + (beta/100);
		}
		findBeta(M, theta);
	}
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