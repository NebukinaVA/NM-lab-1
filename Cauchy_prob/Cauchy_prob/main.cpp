#include "test.h"
#include "task1.h"
#include <iostream>

/*
        x0 = _x0;
		v0 = _v0;
		h = _h;
		n = _n;
		eps = _eps;
		xmax = _xmax;
		prec = _prec;
*/

int main()
{
	test example(0.0, 1.0, 0.001, 1000, 1e-6, 2, 1e-3);
	example.calculate_w_error();
	std::cout << example;
	int i;
	std::cin >> i;
	return 0;
}