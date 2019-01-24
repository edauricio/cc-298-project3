#ifndef BULK_H
#define BULK_H

#include <vector>
#include <string>

class Fluid {
	public:
		typedef std::vector<double>::size_type fluid_index;

		Fluid(fluid_index x_length, fluid_index y_length = 1, fluid_index z_length = 1, double gam = 1.4, double cp = 1.0, double r = 1.0) : Gamma(gam), Cp(cp), Rgas(r), IMAX(x_length), JMAX(y_length), KMAX(z_length) { prop_size(); };
		double &gamma() { return Gamma; };
		double &cp() { return Cp; };
		double &R() { return Rgas; };
		double &rho() { return RhoC; };
		double &rho(fluid_index i = 0, fluid_index j = 0, fluid_index k = 0) { return Rho[i][j][k]; };
		double &mu() { return MuC; };
		double &mu(fluid_index i, fluid_index j = 0, fluid_index k = 0) { return Mu[i][j][k]; };
		double ni(fluid_index i, fluid_index j = 0, fluid_index k = 0) { return Mu[i][j][k]/Rho[i][j][k]; };
	private:
		double Gamma, Cp, Rgas, RhoC, MuC;
		const fluid_index IMAX, JMAX, KMAX;
		std::vector<std::vector<std::vector<double> > > Rho, Mu;
		void prop_size();
		friend class ScalarFlowProp;
		
};

// Variable property sizing (after constructor with IMAX & JMAX provided)
void Fluid::prop_size() {
	Rho.resize(IMAX);
	Mu.resize(IMAX);
	for (fluid_index i = 0; i != IMAX; ++i) {
		Rho[i].resize(JMAX);
		Mu[i].resize(JMAX);
		for (fluid_index j = 0; j != JMAX; ++j) {
			Rho[i][j].resize(KMAX);
			Mu[i][j].resize(KMAX);
		}
	}
}

class ScalarFlowProp {
	public:
		typedef std::vector<double>::size_type prop_index;
		ScalarFlowProp(prop_index x_length, prop_index y_length = 1, prop_index z_length = 1) : IMAX(x_length), JMAX(y_length), KMAX(z_length) { prop_size(); };
		void setName(std::string pName) { PropName = pName; };
		//std::vector<std::vector<double> > &operator[] (prop_index ind) { return Prop[ind]; };
		double &val(prop_index i, prop_index j = 0, prop_index k = 0) { return Prop[i][j][k]; };
	private:
		prop_index IMAX, JMAX, KMAX;
		std::string PropName;
		std::vector<std::vector<std::vector<double> > > Prop;
		void prop_size();
};

void ScalarFlowProp::prop_size() {
	Prop.resize(IMAX);
	for (prop_index i = 0; i != IMAX; ++i) {
		Prop[i].resize(JMAX);
		for (prop_index j = 0; j != JMAX; ++j) {
			Prop[i][j].resize(KMAX);
		}
	}
}

class VectorFlowProp {
	public:
		typedef std::vector<double>::size_type prop_index;
		VectorFlowProp(prop_index x_length, prop_index y_length = 1, prop_index z_length = 1) : IMAX(x_length), JMAX(y_length), KMAX(z_length) { prop_size(); };
		void setName(std::string pName) { PropName = pName; };
	private:
		prop_index IMAX, JMAX, KMAX;
		std::string PropName;
		std::vector<std::vector<std::vector<double> > > PropX, PropY, PropZ;
		void prop_size();
};

void VectorFlowProp::prop_size() {
	PropX.resize(IMAX);
	PropY.resize(IMAX);
	PropZ.resize(IMAX);
}

#endif