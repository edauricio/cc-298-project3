#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>

using namespace std;

class Mesh {
    vector<double>::size_type IMAX, JMAX, KMAX;
    vector<vector<vector<double> > > X, Y, Z, XI_X, XI_Y, XI_Z, ETA_X, ETA_Y, ETA_Z, Jacob;

    void sizing();
  public:
    explicit Mesh(vector<double>::size_type x_length = 1, vector<double>::size_type y_length = 1, vector<double>::size_type z_length = 1) : IMAX(x_length), JMAX(y_length), KMAX(z_length) { sizing(); } ;
    Mesh &read(string);
    const Mesh &write(string, string = "Mesh") const;
    const Mesh &write(string, const vector<double>::size_type, const vector<double>::size_type, const vector<double>::size_type, string = "Mesh") const;
    vector<double>::size_type xsize() const { return IMAX; };
    vector<double>::size_type ysize() const { return JMAX; };
    vector<double>::size_type zsize() const { return KMAX; };
    double x(vector<double>::size_type i, vector<double>::size_type j = 0, vector<double>::size_type k = 0) const { return X[i][j][k]; };
    double y(vector<double>::size_type i, vector<double>::size_type j = 0, vector<double>::size_type k = 0) const { return Y[i][j][k]; };
    double z(vector<double>::size_type i, vector<double>::size_type j = 0, vector<double>::size_type k = 0) const { return Z[i][j][k]; };
    double xi_x(vector<double>::size_type i, vector<double>::size_type j = 0, vector<double>::size_type k = 0) const { return XI_X[i][j][k]; };
    double xi_y(vector<double>::size_type i, vector<double>::size_type j = 0, vector<double>::size_type k = 0) const { return XI_Y[i][j][k]; };
    double xi_z(vector<double>::size_type i, vector<double>::size_type j = 0, vector<double>::size_type k = 0) const { return XI_Z[i][j][k]; };
    double eta_x(vector<double>::size_type i, vector<double>::size_type j = 0, vector<double>::size_type k = 0) const { return ETA_X[i][j][k]; };
    double eta_y(vector<double>::size_type i, vector<double>::size_type j = 0, vector<double>::size_type k = 0) const { return ETA_Y[i][j][k]; };
    double eta_z(vector<double>::size_type i, vector<double>::size_type j = 0, vector<double>::size_type k = 0) const { return ETA_Z[i][j][k]; };
    double J(vector<double>::size_type i, vector<double>::size_type j = 0, vector<double>::size_type k = 0) const { return Jacob[i][j][k]; };
};

void Mesh::sizing() {
  X.resize(IMAX);
  Y.resize(IMAX);
  Z.resize(IMAX);
  XI_X.resize(IMAX);
  XI_Y.resize(IMAX);
  ETA_X.resize(IMAX);
  ETA_Y.resize(IMAX);
  Jacob.resize(IMAX);
  for (auto i = 0; i != IMAX; ++i) {
    X[i].resize(JMAX);
    Y[i].resize(JMAX);
    Z[i].resize(JMAX);
    XI_X[i].resize(JMAX);
    XI_Y[i].resize(JMAX);
    ETA_X[i].resize(JMAX);
    ETA_Y[i].resize(JMAX);
    Jacob[i].resize(JMAX);
    for (auto j = 0; j != JMAX; ++j) {
      X[i][j].resize(KMAX);
      Y[i][j].resize(KMAX);
      Z[i][j].resize(KMAX);
      XI_X[i][j].resize(KMAX);
      XI_Y[i][j].resize(KMAX);
      ETA_X[i][j].resize(KMAX);
      ETA_Y[i][j].resize(KMAX);
      Jacob[i][j].resize(KMAX);
    }
  }
}

Mesh &Mesh::read(string fName) {
  ifstream rFile;
  rFile.open(fName);
  rFile.precision(6);
  rFile >> scientific;
  int aux, header = 1;
  string lineBuff;
  while(header && getline(rFile, lineBuff)) {
    if (lineBuff.substr(0,10) == "DIMENSIONS") {
      IMAX = atoi(lineBuff.substr(11, lineBuff.find_first_of(' ', 11)-11).c_str());
      aux = lineBuff.find_first_of(' ', 11);
      JMAX = atoi(lineBuff.substr(aux+1, lineBuff.find_first_of(' ', aux+1)-(aux+1)).c_str());
      aux = lineBuff.find_first_of(' ', aux+1);
      KMAX = atoi(lineBuff.substr(aux+1, lineBuff.find_first_of(' ', aux+1)-(aux+1)).c_str());
    } else if (lineBuff.substr(0, 6) == "POINTS") header = 0;
  }
  sizing();
  for (vector<double>::size_type k = 0; k != KMAX; ++k) {
    for (vector<double>::size_type j = 0; j != JMAX; ++j) {
      for (vector<double>::size_type i = 0; i != IMAX; ++i) {
        getline(rFile, lineBuff);
        X[i][j][k] = atof(lineBuff.substr(0, lineBuff.find_first_of(' ', 0)).c_str());
        aux = lineBuff.find_first_of(' ', 0);
        Y[i][j][k] = atof(lineBuff.substr(aux+1, lineBuff.find_first_of(' ', aux+1)-(aux+1)).c_str());
        aux = lineBuff.find_first_of(' ', aux+1);
        Z[i][j][k] = atof(lineBuff.substr(aux+1).c_str());
      }
    }
  }
  rFile.close();
  return *this;
}

const Mesh &Mesh::write(string fName, string mName) const {
  ofstream wFile;
  wFile.open(fName);
  wFile.precision(6);
  wFile << scientific;
  wFile << "# vtk DataFile Version 3.0" << endl
        << mName << endl
        << "ASCII" << endl
        << "DATASET STRUCTURED_GRID" << endl
        << "DIMENSIONS " << IMAX << " " << JMAX << " " << KMAX << " " << endl
        << "POINTS " << IMAX*JMAX*KMAX << "double" << endl;
  for (vector<double>::size_type k = 0; k != KMAX; ++k) {
    for (vector<double>::size_type j = 0; j != JMAX; ++j) {
      for (vector<double>::size_type i = 0; i != IMAX; ++i) {
        wFile << X[i][j][k] << " " << Y[i][j][k] << " " << Z[i][j][k] << endl;
      }
    }
  }
  wFile.close();
  return *this;
}

const Mesh &Mesh::write(string fName, const vector<double>::size_type i_end, const vector<double>::size_type j_end, const vector<double>::size_type k_end, string mName) const {
  ofstream wFile;
  wFile.open(fName);
  wFile.precision(6);
  wFile << scientific;
  wFile << "# vtk DataFile Version 3.0" << endl
        << mName << endl
        << "ASCII" << endl
        << "DATASET STRUCTURED_GRID" << endl
        << "DIMENSIONS " << IMAX << " " << JMAX << " " << KMAX << " " << endl
        << "POINTS " << IMAX*JMAX*KMAX << "double" << endl;
  for (vector<double>::size_type k = 0; k != k_end; ++k) {
    for (vector<double>::size_type j = 0; j != j_end; ++j) {
      for (vector<double>::size_type i = 0; i != i_end; ++i) {
        wFile << X[i][j][k] << " " << Y[i][j][k] << " " << Z[i][j][k] << endl;
      }
    }
  }
  wFile.close();
  return *this;
}

#endif