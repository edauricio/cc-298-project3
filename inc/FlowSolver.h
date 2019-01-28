#ifndef FLOWSOLVER_H
#define FLOWSOLVER_H

#include <string>
#include <vector>
#include <initializer_list>
#include "Mesh.h"
#include "MecFlu.h"

bool issign(char);
double DDt(ScalarFlowProp* Prop);
void SolveEqn(const std::string&, std::initializer_list<ScalarFlowProp*>, const Mesh&);
ScalarFlowProp &setInitialFlow(ScalarFlowProp&);
ScalarFlowProp &setBoundaryCond(ScalarFlowProp&, std::vector<Boundary> &, int = 1);
std::vector<Boundary> &defineBoundaries(std::vector<Boundary>&, std::string);

#endif