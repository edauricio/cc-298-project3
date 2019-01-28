#include "MecFlu.h"

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

void ScalarFlowProp::prop_size() {
  Prop.resize(IMAX);
  for (prop_index i = 0; i != IMAX; ++i) {
    Prop[i].resize(JMAX);
    for (prop_index j = 0; j != JMAX; ++j) {
      Prop[i][j].resize(KMAX);
    }
  }
}

void VectorFlowProp::prop_size() {
  PropX.resize(IMAX);
  PropY.resize(IMAX);
  PropZ.resize(IMAX);
  for (prop_index i = 0; i != IMAX; ++i) {
    PropX[i].resize(JMAX);
    PropY[i].resize(JMAX);
    PropZ[i].resize(JMAX);
    for (prop_index j = 0; j != JMAX; ++j) {
      PropX[i][j].resize(KMAX);
      PropY[i][j].resize(KMAX);
      PropZ[i][j].resize(KMAX);      
    }
  }
}