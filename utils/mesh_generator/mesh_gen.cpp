#include <iostream>
#include "Mesh.h"

using namespace std;

int main(int argc, char* argv[]) {
  Mesh mesh("malha0.vtk");
  //mesh.read("malha0.vtk");
  cout.precision(6);
  cout << fixed << mesh.x(30,5) << endl;
  mesh.set_x(30,5) = 101.24;
  cout << mesh.x(30,5) << endl;
  mesh.write("teste_malha.vtk");

  return 0;
}