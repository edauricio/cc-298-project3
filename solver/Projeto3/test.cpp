#include <iostream>
#include "Mesh.h"
#include "MecFlu.h"
#include "FlowSolver.h"

int main() {
  const Mesh Malha("malha.vtk");
  ScalarFlowProp T("Temp", 10,10), p("Pressure", 10, 10);
  SolveEqn("Dxx + Dyy = 0", {&T, &T}, Malha);
  std::cout << T.at(5,9) << std::endl;
  std::cout << p.at(1,1) << std::endl;
}