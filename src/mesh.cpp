#include <string>
#include <vector>
#include <iostream>
#include "Mesh.h"

using namespace std;

void Mesh::sizing() {
  X.resize(IMAX);
  Y.resize(IMAX);
  Z.resize(IMAX);
  XI_X.resize(IMAX);
  XI_Y.resize(IMAX);
  XI_Z.resize(IMAX);
  ETA_X.resize(IMAX);
  ETA_Y.resize(IMAX);
  ETA_Z.resize(IMAX);
  ZETA_X.resize(IMAX);
  ZETA_Y.resize(IMAX);
  ZETA_Z.resize(IMAX);
  Jacob.resize(IMAX);
  for (auto i = 0; i != IMAX; ++i) {
    X[i].resize(JMAX);
    Y[i].resize(JMAX);
    Z[i].resize(JMAX);
    XI_X[i].resize(JMAX);
    XI_Y[i].resize(JMAX);
    XI_Z[i].resize(JMAX);
    ETA_X[i].resize(JMAX);
    ETA_Y[i].resize(JMAX);
    ETA_Z[i].resize(JMAX);
    ZETA_X[i].resize(JMAX);
    ZETA_Y[i].resize(JMAX);
    ZETA_Z[i].resize(JMAX);
    Jacob[i].resize(JMAX);
    for (auto j = 0; j != JMAX; ++j) {
      X[i][j].resize(KMAX);
      Y[i][j].resize(KMAX);
      Z[i][j].resize(KMAX);
      XI_X[i][j].resize(KMAX);
      XI_Y[i][j].resize(KMAX);
      XI_Z[i][j].resize(KMAX);
      ETA_X[i][j].resize(KMAX);
      ETA_Y[i][j].resize(KMAX);
      ETA_Z[i][j].resize(KMAX);
      ZETA_X[i][j].resize(KMAX);
      ZETA_Y[i][j].resize(KMAX);
      ZETA_Z[i][j].resize(KMAX);
      Jacob[i][j].resize(KMAX);
    }
  }
}

void Mesh::metrics() {
  if (coord != "GENERALIZED" && coord != "CARTESIAN") {
    cout << "Invalid coordinate system provided to Constructor. Should be either Cartesian or Generalized." << endl;
    exit(-1);
  } else {
    if (KMAX == 1) {
      // Interior points
      for (mesh_index i = 1; i != IMAX-1; ++i) {
        for (mesh_index j = 1; j != JMAX-1; ++j) {
          Jacob[i][j][0] = 1./(0.5*(X[i+1][j][0] - X[i-1][j][0])*0.5*(Y[i][j+1][0] - Y[i][j-1][0]) - 0.5*(X[i][j+1][0] - X[i][j-1][0])*0.5*(Y[i+1][j][0] - Y[i-1][j][0]));
          XI_X[i][j][0] = Jacob[i][j][0]*0.5*(Y[i][j+1][0] - Y[i][j-1][0]);
          XI_Y[i][j][0] = -Jacob[i][j][0]*0.5*(X[i][j+1][0] - X[i][j-1][0]);
          ETA_X[i][j][0] = -Jacob[i][j][0]*0.5*(Y[i+1][j][0] - Y[i-1][j][0]);
          ETA_Y[i][j][0] = Jacob[i][j][0]*0.5*(X[i+1][j][0] - X[i-1][j][0]);
        }
      }

      // Left boundary
      for (mesh_index j = 1; j != JMAX-1; ++j) {
        Jacob[0][j][0] = 1./((X[1][j][0] - X[0][j][0])*0.5*(Y[0][j+1][0] - Y[0][j-1][0]) - 0.5*(X[0][j+1][0] - X[0][j-1][0])*(Y[1][j][0] - Y[0][j][0]));
        XI_X[0][j][0] = Jacob[0][j][0]*0.5*(Y[0][j+1][0] - Y[0][j-1][0]);
        XI_Y[0][j][0] = -Jacob[0][j][0]*0.5*(X[0][j+1][0] - X[0][j-1][0]);
        ETA_X[0][j][0] = -Jacob[0][j][0]*(Y[1][j][0] - Y[0][j][0]);
        ETA_Y[0][j][0] = Jacob[0][j][0]*(X[1][j][0] - X[0][j][0]);
      }

      // Right boundary
      for (mesh_index j = 1; j != JMAX-1; ++j) {
        Jacob[IMAX-1][j][0] = 1./((X[IMAX-1][j][0] - X[IMAX-2][j][0])*0.5*(Y[IMAX-1][j+1][0] - Y[IMAX-1][j-1][0]) - 0.5*(X[IMAX-1][j+1][0] - X[IMAX-1][j-1][0])*(Y[IMAX-1][j][0] - Y[IMAX-2][j][0]));
        XI_X[IMAX-1][j][0] = Jacob[IMAX-1][j][0]*0.5*(Y[IMAX-1][j+1][0] - Y[IMAX-1][j-1][0]);
        XI_Y[IMAX-1][j][0] = -Jacob[IMAX-1][j][0]*0.5*(X[IMAX-1][j+1][0] - X[IMAX-1][j-1][0]);
        ETA_X[IMAX-1][j][0] = -Jacob[IMAX-1][j][0]*(Y[IMAX-1][j][0] - Y[IMAX-2][j][0]);
        ETA_Y[IMAX-1][j][0] = Jacob[IMAX-1][j][0]*(X[IMAX-1][j][0] - X[IMAX-2][j][0]);
      }

      // Bottom boundary
      for (mesh_index i = 1; i != IMAX-1; ++i) {
        Jacob[i][0][0] = 1./(0.5*(X[i+1][0][0] - X[i-1][0][0])*(Y[i][1][0] - Y[i][0][0]) - (X[i][1][0] - X[i][0][0])*0.5*(Y[i+1][0][0] - Y[i-1][0][0]));
        XI_X[i][0][0] = Jacob[i][0][0]*(Y[i][1][0] - Y[i][0][0]);
        XI_Y[i][0][0] = -Jacob[i][0][0]*(X[i][1][0] - X[i][0][0]);
        ETA_X[i][0][0] = -Jacob[i][0][0]*0.5*(Y[i+1][0][0] - Y[i-1][0][0]);
        ETA_Y[i][0][0] = Jacob[i][0][0]*0.5*(X[i+1][0][0] - X[i-1][0][0]);
      }

      // Upper boundary
      for (mesh_index i = 1; i != IMAX-1; ++i) {
        Jacob[i][JMAX-1][0] = 1./(0.5*(X[i+1][JMAX-1][0] - X[i-1][JMAX-1][0])*(Y[i][JMAX-1][0] - Y[i][JMAX-2][0]) - (X[i][JMAX-1][0] - X[i][JMAX-2][0])*0.5*(Y[i+1][JMAX-1][0] - Y[i-1][JMAX-1][0]));
        XI_X[i][JMAX-1][0] = Jacob[i][JMAX-1][0]*(Y[i][JMAX-1][0] - Y[i][JMAX-2][0]);
        XI_Y[i][JMAX-1][0] = -Jacob[i][JMAX-1][0]*(X[i][JMAX-1][0] - X[i][JMAX-2][0]);
        ETA_X[i][JMAX-1][0] = -Jacob[i][JMAX-1][0]*0.5*(Y[i+1][JMAX-1][0] - Y[i-1][JMAX-1][0]);
        ETA_Y[i][JMAX-1][0] = Jacob[i][JMAX-1][0]*0.5*(X[i+1][JMAX-1][0] - X[i-1][JMAX-1][0]);
      }

      // Corners
      Jacob[0][0][0] = 1./((X[1][0][0] - X[0][0][0])*(Y[0][1][0] - Y[0][0][0]) - (X[0][1][0] - X[0][0][0])*(Y[1][0][0] - Y[0][0][0]));
      XI_X[0][0][0] = Jacob[0][0][0]*(Y[0][1][0] - Y[0][0][0]);
      XI_Y[0][0][0] = -Jacob[0][0][0]*(X[0][1][0] - X[0][0][0]);
      ETA_X[0][0][0] = -Jacob[0][0][0]*(Y[1][0][0] - Y[0][0][0]);
      ETA_Y[0][0][0] = Jacob[0][0][0]*(X[1][0][0] - X[0][0][0]);

      Jacob[IMAX-1][0][0] = 1./((X[IMAX-1][0][0] - X[IMAX-2][0][0])*(Y[IMAX-1][1][0] - Y[IMAX-1][0][0]) - (X[IMAX-1][1][0] - X[IMAX-1][0][0])*(Y[IMAX-1][0][0] - Y[IMAX-2][0][0]));
      XI_X[IMAX-1][0][0] = Jacob[IMAX-1][0][0]*(Y[IMAX-1][1][0] - Y[IMAX-1][0][0]);
      XI_Y[IMAX-1][0][0] = -Jacob[IMAX-1][0][0]*(X[IMAX-1][1][0] - X[IMAX-1][0][0]);
      ETA_X[IMAX-1][0][0] = -Jacob[IMAX-1][0][0]*(Y[IMAX-1][0][0] - Y[IMAX-2][0][0]);
      ETA_Y[IMAX-1][0][0] = Jacob[IMAX-1][0][0]*(X[IMAX-1][0][0] - X[IMAX-2][0][0]);

      Jacob[0][JMAX-1][0] = 1./((X[1][JMAX-1][0] - X[0][JMAX-1][0])*(Y[0][JMAX-1][0] - Y[0][JMAX-2][0]) - (X[0][JMAX-1][0] - X[0][JMAX-2][0])*(Y[1][JMAX-1][0] - Y[0][JMAX-1][0]));
      XI_X[0][JMAX-1][0] = Jacob[0][JMAX-1][0]*(Y[0][JMAX-1][0] - Y[0][JMAX-2][0]);
      XI_Y[0][JMAX-1][0] = -Jacob[0][JMAX-1][0]*(X[0][JMAX-1][0] - X[0][JMAX-2][0]);
      ETA_X[0][JMAX-1][0] = -Jacob[0][JMAX-1][0]*(Y[1][JMAX-1][0] - Y[0][JMAX-1][0]);
      ETA_Y[0][JMAX-1][0] = Jacob[0][JMAX-1][0]*(X[1][JMAX-1][0] - X[0][JMAX-1][0]);

      Jacob[IMAX-1][JMAX-1][0] = 1./((X[IMAX-1][JMAX-1][0] - X[IMAX-2][JMAX-1][0])*(Y[IMAX-1][JMAX-1][0] - Y[IMAX-1][JMAX-2][0]) - (X[IMAX-1][JMAX-1][0] - X[IMAX-1][JMAX-2][0])*(Y[IMAX-1][JMAX-1][0] - Y[IMAX-2][JMAX-1][0]));
      XI_X[IMAX-1][JMAX-1][0] = Jacob[IMAX-1][JMAX-1][0]*(Y[IMAX-1][JMAX-1][0] - Y[IMAX-1][JMAX-2][0]);
      XI_Y[IMAX-1][JMAX-1][0] = -Jacob[IMAX-1][JMAX-1][0]*(X[IMAX-1][JMAX-1][0] - X[IMAX-1][JMAX-2][0]);
      ETA_X[IMAX-1][JMAX-1][0] = -Jacob[IMAX-1][JMAX-1][0]*(Y[IMAX-1][JMAX-1][0] - Y[IMAX-2][JMAX-1][0]);
      ETA_Y[IMAX-1][JMAX-1][0] = Jacob[IMAX-1][JMAX-1][0]*(X[IMAX-1][JMAX-1][0] - X[IMAX-2][JMAX-1][0]);
    } else {
      cout << "3D mesh not supported yet." << endl;
      exit(-1);
    }
  }
}

Mesh &Mesh::read(string fName) {
  ifstream rFile;
  rFile.open(fName);
  if (rFile.is_open()) {
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
    for (mesh_index k = 0; k != KMAX; ++k) {
      for (mesh_index j = 0; j != JMAX; ++j) {
        for (mesh_index i = 0; i != IMAX; ++i) {
          getline(rFile, lineBuff);
          X[i][j][k] = atof(lineBuff.substr(0, lineBuff.find_first_of(' ', 0)).c_str());
          aux = lineBuff.find_first_of(' ', 0);
          Y[i][j][k] = atof(lineBuff.substr(aux+1, lineBuff.find_first_of(' ', aux+1)-(aux+1)).c_str());
          aux = lineBuff.find_first_of(' ', aux+1);
          Z[i][j][k] = atof(lineBuff.substr(aux+1).c_str());
        }
      }
    }
  } else {
    cout << "Error opening file to read mesh. Please check filename." << endl;
    exit(-1);
  }
  rFile.close();
  return *this;
}

const Mesh &Mesh::write(string fName, string pMetrics, string mName) const {
  for (auto &c : pMetrics)
    c = toupper(c);
  ofstream wFile;
  wFile.open(fName);
  wFile.precision(6);
  wFile << scientific;
  wFile << "# vtk DataFile Version 3.0" << endl
        << mName << endl
        << "ASCII" << endl
        << "DATASET STRUCTURED_GRID" << endl
        << "DIMENSIONS " << IMAX << " " << JMAX << " " << KMAX << " " << endl
        << "POINTS " << IMAX*JMAX*KMAX << " double" << endl;
  for (mesh_index k = 0; k != KMAX; ++k) {
    for (mesh_index j = 0; j != JMAX; ++j) {
      for (mesh_index i = 0; i != IMAX; ++i) {
        wFile << X[i][j][k] << " " << Y[i][j][k] << " " << Z[i][j][k] << endl;
      }
    }
  }
  if (pMetrics == "YES") {
    wFile << "POINT_DATA " << IMAX*JMAX*KMAX << endl
          << "SCALARS Jacob double" << endl
          << "LOOKUP_TABLE default" << endl;
    for (mesh_index k = 0; k != KMAX; ++k) {
      for (mesh_index j = 0; j != JMAX; ++j) {
        for (mesh_index i = 0; i != IMAX; ++i) {
          wFile << Jacob[i][j][k] << endl;
        }
      }
    }
    wFile << "SCALARS xi_x double" << endl
          << "LOOKUP_TABLE default" << endl;
    for (mesh_index k = 0; k != KMAX; ++k) {
      for (mesh_index j = 0; j != JMAX; ++j) {
        for (mesh_index i = 0; i != IMAX; ++i) {
          wFile << XI_X[i][j][k] << endl;
        }
      }
    }
    wFile << "SCALARS xi_y double" << endl
          << "LOOKUP_TABLE default" << endl;
    for (mesh_index k = 0; k != KMAX; ++k) {
      for (mesh_index j = 0; j != JMAX; ++j) {
        for (mesh_index i = 0; i != IMAX; ++i) {
          wFile << XI_Y[i][j][k] << endl;
        }
      }
    }
    wFile << "SCALARS eta_x double" << endl
          << "LOOKUP_TABLE default" << endl;
    for (mesh_index k = 0; k != KMAX; ++k) {
      for (mesh_index j = 0; j != JMAX; ++j) {
        for (mesh_index i = 0; i != IMAX; ++i) {
          wFile << ETA_X[i][j][k] << endl;
        }
      }
    }
    wFile << "SCALARS eta_y double" << endl
          << "LOOKUP_TABLE default" << endl;
    for (mesh_index k = 0; k != KMAX; ++k) {
      for (mesh_index j = 0; j != JMAX; ++j) {
        for (mesh_index i = 0; i != IMAX; ++i) {
          wFile << ETA_Y[i][j][k] << endl;
        }
      }
    }
  }
  wFile.close();
  return *this;
}