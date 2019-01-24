#include <iostream>
#include "Mesh.h"

using namespace std;

int main(int argc, char* argv[]) {
  Mesh mesh;
  mesh.read("malha0.vtk");
  cout << mesh.x(339, 1) << endl;
  mesh.write("teste_malha.vtk");

  return 0;
}