#include <iostream>
#include "Bulk.h"

int main() {
	ScalarFlowProp Q_(10, 10);
	std::cout << Q_.val(0,1) << std::endl;
	Q_.val(0,1) = 0.12;
	std::cout << Q_.val(0,1) << std::endl;
}