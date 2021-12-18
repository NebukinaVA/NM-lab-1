#include "test.h"
#include "task1.h"
#include "task2.h"
#include <iostream>

int main()
{
//	int i;
	setlocale(LC_ALL, "Russian");
//	test example(0.0, 1.0, 0.001, 1000, 1e-5, 2, 1e-3);
//	task1 example(0.0, 1.0, 0.001, 200, 1e-5, 2, 1e-3);
//	task2 example(0.0, 1.0, 1.0, 2.0, 1.0, 1.0, 0.001, 1000, 1e-5, 10.0, 1e-3);
//	example.calculate_w_error();
//	std::cout << example;
//	std::cin >> i;

	// ------------------------

	
	int task, n;
	double x0, u0, v0, h, eps, xmax, prec, a, b, c;
	eps = 1e-5;
	bool flag = true;
	while (flag) {
		std::cout << "�������� ��� ������:\n�������� - 0\n�������� 1 - 1\n�������� 2 - 2\n";
		std::cin >> task;
		if (task == 0 || task == 1)
		{
			std::cout << "������� ��������� �������:\n";
			std::cout << "x0 = ";
			std::cin >> x0;
			std::cout << "u0 = ";
			std::cin >> u0;
		}
		else
		{
			std::cout << "������� ��������� �������:\n";
			std::cout << "x0 = ";
			std::cin >> x0;
			std::cout << "u0 = ";
			std::cin >> u0;
			std::cout << "u'0 = ";
			std::cin >> v0;
			std::cout << "������� ��������� �������:\n";
			std::cout << "a = ";
			std::cin >> a;
			std::cout << "b = ";
			std::cin >> b;
			std::cout << "c = ";
			std::cin >> c;
		}
		std::cout << "������� ��������� ���:\nh = ";
		std::cin >> h;
		std::cout << "������� ������������ ����� �����:\nN = ";
		std::cin >> n;
		std::cout << "������� ������ ������� �� �:\nxmax = ";
		std::cin >> xmax;
		std::cout << "������� �������� ������ �� ������ �������:\nprecision = ";
		std::cin >> prec;
		bool error;
		std::cout << "����������� ������ � ������� ��������� �����������?\n�� - 1\n��� - 0\n";
		std::cin >> error;
		if (error)
		{
			std::cout << "������� �������� ��������� �����������:\nepsilon = ";
			std::cin >> eps;
		}
		if (task == 0)
		{
			test equation(x0, u0, h, n, eps, xmax, prec);
			if (error)
				equation.calculate_w_error();
			else equation.calculate();
			std::cout << equation;
			equation.help();
		}
		else if (task == 1)
		{
			task1 equation(x0, u0, h, n, eps, xmax, prec);
			if (error)
				equation.calculate_w_error();
			else equation.calculate();
			std::cout << equation;
			equation.help();
		}
		else if (task == 2)
		{
			task2 equation(x0, u0, v0, a, b, c, h, n, eps, xmax, prec);
			if (error)
				equation.calculate_w_error();
			else equation.calculate();
			std::cout << equation;
			equation.help();
		}
		else
		{
			std::cout << "�������� ������";
			flag = false;
		}
		std::cout << "\n������ ���������� ������?\n�� - 1\n��� - 0\n";
		std::cin >> flag;
	}
	
	return 0;
}