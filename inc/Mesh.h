#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cctype>

class Mesh {
  public:
  	typedef std::vector<double>::size_type mesh_index;
  	Mesh(std::string fName) { read(fName); metrics(); }
  	Mesh(std::string fName, std::string syscoord) : coord(syscoord) { for (auto &c : coord) c = toupper(c); read(fName); metrics(); };
    explicit Mesh(mesh_index x_length = 1, mesh_index y_length = 1, mesh_index z_length = 1) : IMAX(x_length), JMAX(y_length), KMAX(z_length) { sizing(); } ;
    const Mesh &write(std::string, std::string = "NO", std::string = "Mesh") const;
    mesh_index xsize() const { return IMAX; };
    mesh_index ysize() const { return JMAX; };
    mesh_index zsize() const { return KMAX; };
    double x(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return X[i][j][k]; };
    double y(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return Y[i][j][k]; };
    double z(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return Z[i][j][k]; };
    double xi_x(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return XI_X[i][j][k]; };
    double xi_y(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return XI_Y[i][j][k]; };
    double xi_z(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return XI_Z[i][j][k]; };
    double eta_x(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return ETA_X[i][j][k]; };
    double eta_y(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return ETA_Y[i][j][k]; };
    double eta_z(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return ETA_Z[i][j][k]; };
    double zeta_x(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return ZETA_X[i][j][k]; };
    double zeta_y(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return ZETA_Y[i][j][k]; };
    double zeta_z(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return ZETA_Z[i][j][k]; };
    double J(mesh_index i, mesh_index j = 0, mesh_index k = 0) const { return Jacob[i][j][k]; };
    std::string coord_sys() const { return coord; };
    void calc_metrics() { metrics(); };

    double &set_x(mesh_index i, mesh_index j = 0, mesh_index k = 0) { return X[i][j][k]; };
    double &set_y(mesh_index i, mesh_index j = 0, mesh_index k = 0) { return Y[i][j][k]; };
    double &set_z(mesh_index i, mesh_index j = 0, mesh_index k = 0) { return Z[i][j][k]; };
  private:
  	mesh_index IMAX, JMAX, KMAX;
    std::vector<std::vector<std::vector<double> > > X, Y, Z, XI_X, XI_Y, XI_Z, ETA_X, ETA_Y, ETA_Z, ZETA_X, ZETA_Y, ZETA_Z, Jacob;
    std::string coord = "CARTESIAN";

    void sizing();
    void metrics();
    Mesh &read(std::string);
};

class Boundary {
  public:
    Boundary(std::string pName, std::string loc = "", std::string type = "", double val = 0, double a = 1, double b = 0) : PropName(pName), Location(loc), Type(type), Value(val), Alfa(a), Beta(b) {};
    std::string &type() { return Type; };
    std::string &propName() { return PropName; };
    std::string &location() { return Location; };
    double &value() { return Value; };
    double &alfa() { return Alfa; };
    double &beta() { return Beta; };
  private:
    std::string PropName, Location, Type;
    double Value, Alfa, Beta;
};

#endif