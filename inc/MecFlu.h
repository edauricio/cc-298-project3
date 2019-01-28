#ifndef MECFLU_H
#define MECFLU_H

#include <vector>
#include <string>

class Fluid {
	public:
		typedef std::vector<double>::size_type fluid_index;

		explicit Fluid(fluid_index x_length, fluid_index y_length = 1, fluid_index z_length = 1, double gam = 1.4, double cp = 1.0, double r = 1.0) : Gamma(gam), Cp(cp), Rgas(r), IMAX(x_length), JMAX(y_length), KMAX(z_length) { prop_size(); };
		double &gamma() { return Gamma; };
		double &cp() { return Cp; };
		double &R() { return Rgas; };
		double &rho() { return RhoC; };
		double &rho(fluid_index i = 0, fluid_index j = 0, fluid_index k = 0) { return Rho[i][j][k]; };
		double &mu() { return MuC; };
		double &mu(fluid_index i, fluid_index j = 0, fluid_index k = 0) { return Mu[i][j][k]; };
		double ni(fluid_index i, fluid_index j = 0, fluid_index k = 0) { return Mu[i][j][k]/Rho[i][j][k]; };
    fluid_index xsize() { return IMAX; };
    fluid_index ysize() { return JMAX; };
    fluid_index zsize() { return KMAX; };
	private:
		double Gamma, Cp, Rgas, RhoC, MuC;
		const fluid_index IMAX, JMAX, KMAX;
		std::vector<std::vector<std::vector<double> > > Rho, Mu;
		void prop_size();
		friend class ScalarFlowProp;
};

class ScalarFlowProp {
	public:
		typedef std::vector<double>::size_type prop_index;
		ScalarFlowProp(std::string pName, prop_index x_length, prop_index y_length = 1, prop_index z_length = 1) : IMAX(x_length), JMAX(y_length), KMAX(z_length), PropName(pName) { prop_size(); };
		void setName(std::string pName) { PropName = pName; };
    std::string name() { return PropName; };
		double &at(prop_index i, prop_index j = 0, prop_index k = 0) { return Prop[i][j][k]; };
    prop_index xsize() { return IMAX; };
    prop_index ysize() { return JMAX; };
    prop_index zsize() { return KMAX; };
	private:
		prop_index IMAX, JMAX, KMAX;
		std::string PropName;
		std::vector<std::vector<std::vector<double> > > Prop;
		void prop_size();
};


class VectorFlowProp {
	public:
		typedef std::vector<double>::size_type prop_index;
		VectorFlowProp(std::string pName, prop_index x_length, prop_index y_length = 1, prop_index z_length = 1) : IMAX(x_length), JMAX(y_length), KMAX(z_length), PropName(pName) { prop_size(); };
		void setName(std::string pName) { PropName = pName; };
    std::string name() { return PropName; };
    double &X(prop_index i, prop_index j = 0, prop_index k = 0) { return PropX[i][j][k]; };
    double &Y(prop_index i, prop_index j = 0, prop_index k = 0) { return PropY[i][j][k]; };
    double &Z(prop_index i, prop_index j = 0, prop_index k = 0) { return PropZ[i][j][k]; };
    prop_index xsize() { return IMAX; };
    prop_index ysize() { return JMAX; };
    prop_index zsize() { return KMAX; };
	private:
		prop_index IMAX, JMAX, KMAX;
		std::string PropName;
		std::vector<std::vector<std::vector<double> > > PropX, PropY, PropZ;
		void prop_size();
};

#endif