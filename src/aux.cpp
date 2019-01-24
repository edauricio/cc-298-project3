#include <iostream>
#include "Aux.h"

void help() {
	std::cout << "Usage: ./meshgen <-I #> <-J #>"	 << std::endl << std::endl
						<< "-I,\t\tNumber of points in X" << std::endl
						<< "-J,\t\tNumber of points in Y" << std::endl;
}