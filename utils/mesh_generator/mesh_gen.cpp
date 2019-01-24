#include <iostream>
#include "Mesh.h"
#include "Aux.h"

using namespace std;

int main(int argc, char* argv[]) {
	Mesh::mesh_index IMAX, JMAX;

	// Checking arguments
	if (argc > 1) {
		if (!(argc%2) || (argc < 5)) { help(); exit(-1); }
		for (auto i = 1; i != argc; ++i) {
			if (string(argv[i]) == "-I") IMAX = atoi(argv[i+1]);
			if (string(argv[i]) == "-J") JMAX = atoi(argv[i+1]);
		}
	} else {
		help();
		exit(-1);
	}

	// Start
	Mesh mesh(IMAX, JMAX);
	double L = 3.5, H = 1.0;
	double dx = L/(IMAX-1), dy = H/(JMAX-1);

	for (Mesh::mesh_index i = 0; i != IMAX; ++i) {
		for (Mesh::mesh_index j = 0; j != JMAX; ++j) {
			mesh.set_x(i,j) = i*dx;
			mesh.set_y(i,j) = j*dy;
		}
	}

	mesh.calc_metrics();
	mesh.write("malha.vtk");
  return 0;
}