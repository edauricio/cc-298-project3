#include <initializer_list>
#include <iostream>
#include <cctype>
#include <fstream>
#include <string>
#include "Mesh.h"
#include "MecFlu.h"
#include "FlowSolver.h"
#include "Numerics.h"
#include "json.hpp"

using json = nlohmann::json;

bool issign(char c) {
  return c == '+' || c == '-' || c == '/' || c == '*';
}

double DDt(ScalarFlowProp* Prop) {
  Prop->at(1,1) = 0.3;
  return 1.1;
}


void SolveEqn(const std::string &Eqn, std::initializer_list<ScalarFlowProp*> ScalarProp, const Mesh& mesh) {
  std::vector<std::string> LHS, RHS;
  std::vector<std::string>::size_type eqSign;
  std::string FullLHS, FullRHS;

  eqSign = Eqn.find_first_of('=');
  FullLHS = Eqn.substr(0, eqSign);
  FullRHS = Eqn.substr(eqSign+1);

  for (std::string::size_type i = 0; i != FullLHS.size(); ++i) {
    if (isspace(FullLHS[i])) {
      FullLHS = FullLHS.substr(0, i) + FullLHS.substr(i+1);
      --i;
    }
  }
  for (std::string::size_type i = 0; i != FullRHS.size(); ++i) {
    if (isspace(FullRHS[i])) {
      FullRHS = FullRHS.substr(0, i) + FullRHS.substr(i+1);
      --i;
    }
  }

  // Vector to hold direction of derivative (space/time) and order of all variables
  std::vector<std::vector<int> > derivativeOrder(ScalarProp.size());
  for (auto &v : derivativeOrder)
    v.resize(4);

  // Parsing LHS operators (signs)
  std::vector<std::string::size_type> Lsigns;
  for (std::string::size_type i = 0; i != FullLHS.size(); ++i) {
    if (issign(FullLHS[i])) Lsigns.push_back(i);
  }
  Lsigns.push_back(FullLHS.size());

  // Parsing LHS derivative order
  size_t firstC = 0;
  std::vector<int>::size_type propN = 0;
  int termsL = 0, termsR = 0, termsT = 0;
  for (std::vector<std::string::size_type>::size_type k = 0; k != Lsigns.size(); ++k) {
    std::string term = FullLHS.substr(firstC, Lsigns[k] - firstC);
    termsL++;
    termsT++;
    if (termsL > ScalarProp.size()) { 
      std::cout << "In " << __func__ << "():" << std::endl
                << "Mismatch between number of terms in the equation and arguments (variables) passed." << std :: endl; 
      exit(-1); 
    }
    int first = 1;
    for (auto e : term) {
      if (first && e != 'D') {
        std::cout << "In " << __func__ << "():" << std::endl
                  << "Equation written in wrong format. Please refer to the documentation." << std::endl;
        exit(-1);
      } else if (first) {
        first = 0;
      } else if (e == 't') {
        derivativeOrder[propN][0] += 1;
      } else if (e == 'x') {
        derivativeOrder[propN][1] += 1;
      } else if (e == 'y') {
        derivativeOrder[propN][2] += 1;
      } else if (e == 'z') {
        derivativeOrder[propN][3] += 1;
      }
    }
    propN++;
    firstC = Lsigns[k]+1;
  }

  // Parsing RHS operators (signs)
  std::vector<std::string::size_type> Rsigns;
  for (std::string::size_type i = 0; i != FullRHS.size(); ++i) {
    if (issign(FullRHS[i])) Rsigns.push_back(i);
  }
  Rsigns.push_back(FullRHS.size());

  // Parsing LHS derivative order
  firstC = 0;
  for (std::vector<std::string::size_type>::size_type k = 0; k != Rsigns.size(); ++k) {
    std::string term = FullRHS.substr(firstC, Rsigns[k] - firstC);
    if (term == "0") break;
    termsR++;
    termsT++;
    if (termsT > ScalarProp.size()) { 
      std::cout << "Numbers of terms and derivatives provided to SolveEqn() don't match." << std :: endl; 
      exit(-1); 
    }
    int first = 1;
    for (auto e : term) {
      if (first && e != 'D') {
        std::cout << "Equation written in wrong format. Please refer to the documentation." << std::endl;
        exit(-1);
      } else if (first) {
        first = 0;
      } else if (e == 't') {
        derivativeOrder[propN][0] += 1;
      } else if (e == 'x') {
        derivativeOrder[propN][1] += 1;
      } else if (e == 'y') {
        derivativeOrder[propN][2] += 1;
      } else if (e == 'z') {
        derivativeOrder[propN][3] += 1;
      }
    }
    propN++;
    firstC = Rsigns[k]+1;
  }

  // Checking if there are more parameters than terms in the equation
  if (ScalarProp.size() > termsT) {
    std::cout << "In " << __func__ << "():" << std::endl
              << "More arguments (variables) passed than terms in the equation." << std::endl
              << "Ignoring spare arguments." << std::endl;
  }

  // Parsing method for solving the equation
  Method scheme("numerics");
  for (std::initializer_list<ScalarFlowProp*>::size_type i = 0; i != ScalarProp.size(); ++i) {
    if (derivativeOrder[i][0]) {
      json ImpData;
      std::ifstream inF(scheme.file());
      ImpData = json::parse(inF);
      scheme.expimp() = ImpData["Method"]["TimeType"];
      inF.close();
    }
  }

  // Creating Boundary Objects
  std::vector<vector<Boundary> > Boundaries(ScalarProp.size());
  std::vector<vector<Boundary> >::size_type B_index = 0;
  for (std::initializer_list<ScalarFlowProp*>::iterator it = ScalarProp.begin(); it != ScalarProp.end(); ++it) {
    defineBoundaries(Boundaries[B_index++], (*it)->name());
  }

  // Setting Initial Conditions
  for (std::initializer_list<ScalarFlowProp*>::iterator it = ScalarProp.begin(); it != ScalarProp.end(); ++it) {
    setInitialFlow(**it);
  }

  // Setting Boundary Values
  B_index = 0;
  for (std::initializer_list<ScalarFlowProp*>::iterator it = ScalarProp.begin(); it != ScalarProp.end(); ++it) {
    setBoundaryCond(**it, Boundaries[B_index++], 0);
  }


  // Starting solver
}

std::vector<Boundary> &defineBoundaries(std::vector<Boundary> &B, std::string pName) {
  std::vector<std::string> Location = {"Top", "Bottom", "Left", "Right"};
  Boundary obj(pName);
  json Data;
  std::ifstream inFile("flow");
  Data = json::parse(inFile);

  for (std::vector<std::string>::size_type i = 0; i != Location.size(); ++i) {
    obj.location() = Location[i];
    obj.type() = Data["BoundaryCond"][obj.location()][obj.propName()]["type"];
    obj.value() = Data["BoundaryCond"][obj.location()][obj.propName()]["value"];
    if (obj.type() == "Robin") {
      obj.beta() = Data["BoundaryCond"][obj.location()][obj.propName()]["beta"];
      obj.alfa() = Data["BoundaryCond"][obj.location()][obj.propName()]["alfa"];
    }
    B.push_back(obj);
  }

  inFile.close();

  return B;
}

ScalarFlowProp &setInitialFlow(ScalarFlowProp& Prop) {
  json Data;
  std::ifstream inFile("flow");
  Data = json::parse(inFile);

  if (Prop.zsize() == 1) {
  for (ScalarFlowProp::prop_index i = 1; i != Prop.xsize()-1; ++i) {
    for (ScalarFlowProp::prop_index j = 1; j != Prop.ysize()-1; ++j) {
        Prop.at(i,j) = Data["InitCond"][Prop.name()];
      }
    }
  } else {
    std::cout << "3D not supported yet." << std::endl;
    exit(-1);
  }

  inFile.close();

  return Prop;
}

ScalarFlowProp &setBoundaryCond(ScalarFlowProp &Prop, std::vector<Boundary> &BC, int cnt) {
  if (Prop.zsize() == 1) {
    for (std::vector<Boundary>::iterator vit = BC.begin(); vit != BC.end(); ++vit) {
      
      if ((vit->type() == "Dirichlet") && cnt == 0) {
        if (vit->location() == "Top") {
          for (ScalarFlowProp::prop_index i = 0; i != Prop.xsize(); ++i) {
            Prop.at(i,Prop.ysize()-1) = vit->value();
          }
        }
        if (vit->location() == "Bottom") {
          for (ScalarFlowProp::prop_index i = 0; i != Prop.xsize(); ++i) {
            Prop.at(i,0) = vit->value();
          }
        }
        if (vit->location() == "Left") {
          for (ScalarFlowProp::prop_index j = 0; j != Prop.ysize(); ++j) {
            Prop.at(0,j) = vit->value();
          }
        }
        if (vit->location() == "Right") {
          for (ScalarFlowProp::prop_index j = 0; j != Prop.ysize(); ++j) {
            Prop.at(Prop.xsize()-1,j) = vit->value();
          }
        }
      }

      else if (vit->type() == "Neumann") {

      }
    }
  } else {
    std::cout << "3D not supported yet." << std::endl;
    exit(-1);
  }

  return Prop;
}