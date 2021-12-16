#include "test.h"
#include "task1.h"
#include "task2.h"
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
//	test example(0.0, 1.0, 0.001, 1000, 1e-5, 2, 1e-3);
//	task1 example(0.0, 1.0, 0.001, 1000, 1e-5, 2, 1e-3);
	task2 example(0.0, 1.0, 1.0, 2.0, 1.0, 1.0, 0.001, 1000, 1e-5, 10.0, 1e-3);
	example.calculate_w_error();
	std::cout << example;
	int i;
	std::cin >> i;
	return 0;
}