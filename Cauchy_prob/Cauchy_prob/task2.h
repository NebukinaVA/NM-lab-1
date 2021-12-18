#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>

class task2
{
private:
	double u0, v0, x0;
	double a, b, c;
	double h, eps, xmax, prec;
	int n, N;                   // number of steps, total steps
	std::vector<double> arg; //x
	std::vector<double> ures; //v
	std::vector<double> vres; //u
	double Xn, Un, Vn;
	int inc = 0;
	int dec = 0; // total inc dec
	int Smin = 0, Smax = 0; // number of string with Smin, Smax
	int hmax, hmin; //number of string with hmax, hmin
	std::vector<std::pair<double, double>> reswcap; //(u, v) with cap
	std::vector<double> steps; //h
	std::vector<double> ss;    //S*
	std::vector<int> hinc;  // total step increases
	std::vector<int> hdec;  // total step decreases
	double func1(double x, double u, double v) // u' = v
	{
		return v;
	}
	double func2(double x, double u, double v) // v'= u'' = -a * v*abs(v) - b * v - c * u
	{
		return (-a * v*abs(v) - b * v - c * u);

	}
public:
	task2(double _x0, double _u0, double _v0, double _a, double _b, double _c, double _h, int _n, double _eps, double _xmax, double _prec)
	{
		x0 = _x0;
		u0 = _u0;
		v0 = _v0;
		a = _a;
		b = _b;
		c = _c;
		h = _h;
		n = _n;
		eps = _eps;
		xmax = _xmax;
		prec = _prec;
	}
	void reset(double _x0, double _u0, double _v0, double _a, double _b, double _c, double _h, int _n, double _eps, double _xmax, double _prec)
	{
		arg.clear();
		ures.clear();
		vres.clear();
		reswcap.clear();
		steps.clear();
		ss.clear();
		hinc.clear();
		hdec.clear();
		N = 0;
		inc = 0;
		dec = 0;
		Smin = 0;
		Smax = 0;
		x0 = _x0;
		u0 = _u0;
		v0 = _v0;
		a = _a;
		b = _b;
		c = _c;
		h = _h;
		n = _n;
		eps = _eps;
		xmax = _xmax;
		prec = _prec;
	}
	std::pair<double, double> RK4(double xn, double un, double vn, double h)
	{
		double ku1, ku2, ku3, ku4;  // for un
		double kv1, kv2, kv3, kv4;  // for vn

		ku1 = func1(xn, un, vn);
		kv1 = func2(xn, un, vn);

		ku2 = func1(xn + h / 2.0, un + (h / 2.0)*ku1, vn + (h / 2.0)*kv1);
		kv2 = func2(xn + h / 2.0, un + (h / 2.0)*ku1, vn + (h / 2.0)*kv1);

		ku3 = func1(xn + h / 2.0, un + (h / 2.0)*ku2, vn + (h / 2.0)*kv2);
		kv3 = func2(xn + h / 2.0, un + (h / 2.0)*ku2, vn + (h / 2.0)*kv2);

		ku4 = func1(xn + h, un + h * ku3, vn + h * kv3);
		kv4 = func2(xn + h, un + h * ku3, vn + h * kv3);

		un = un + h * (ku1 + 2.0 * ku2 + 2.0 * ku3 + ku4) / 6.0;
		vn = vn + h * (kv1 + 2.0 * kv2 + 2.0 * kv3 + kv4) / 6.0;

		return std::make_pair(un, vn);
	}
	std::pair<std::vector<double>, std::vector<double>> calculate()
	{
		arg.push_back(x0);
		auto result = std::make_pair(u0, v0);
		ures.push_back(u0);
		vres.push_back(v0);
		hinc.push_back(0);
		ss.push_back(0.0);
		hdec.push_back(0);
		steps.push_back(0.0);
		reswcap.push_back(std::make_pair(0.0, 0.0));
		double xn = x0;
		int i = 0;
		while (i < n)
		{
			if ((xn > (xmax - prec)) && (xn < xmax))
			{
				break;
			}
			else {
				if ((xn + h) > xmax)
				{
					while (((xn + h) > xmax) && (xn < (xmax - prec)))
					{
						h /= 2.0;
					}
					result = RK4(xn, result.first, result.second, h);
					xn += h;
				}
				else
				{
					result = RK4(xn, result.first, result.second, h);
					xn += h;
				}
				arg.insert(arg.begin() + i + 1, xn);
				ures.insert(ures.begin() + i + 1, result.first);
				vres.insert(vres.begin() + i + 1, result.second);
				steps.insert(steps.begin() + i + 1, h);
				hinc.insert(hinc.begin() + i + 1, 0);
				hdec.insert(hdec.begin() + i + 1, 0);
				reswcap.insert(reswcap.begin() + i + 1, std::make_pair(0.0, 0.0));
				ss.insert(ss.begin() + i + 1, 0.0);
				i++;
			}
		}
		N = i;
		Xn = arg[i];
		Un = ures[i];
		Vn = vres[i];
		if (i == n)
			hmin = 1;
		else hmin = i;
		hmax = 1;
		std::pair<std::vector<double>, std::vector<double>> res;
		res.first = ures;
		res.second = vres;
		return res;
	}
	std::pair<std::vector<double>, std::vector<double>> calculate_w_error()
	{
		auto result = std::make_pair(u0, v0);
		ss.push_back(0.0);
		arg.push_back(x0);
		ures.push_back(u0);
		vres.push_back(v0);
		steps.push_back(0.0);
		reswcap.push_back(std::make_pair(0.0, 0.0));
		hinc.push_back(0);
		hdec.push_back(0);
		hinc.push_back(0);
		hdec.push_back(0);
		double xn = x0;
		double xhalf = x0;
		auto half = result;
		auto cap = result;
		auto next = result;
		double S, S1, S2;
		int i = 0;
		while (i < n)
		{
			if ((xn > (xmax - prec)) && (xn <= xmax))
			{
				break;
			}
			else
			{
				next = RK4(xn, result.first, result.second, h);
				half = RK4(xn, result.first, result.second, h / 2.0);
				xhalf = xn + h / 2.0;
				cap = RK4(xhalf, half.first, half.second, h / 2.0);
				S1 = abs(cap.first - next.first) / 15.0;
				S2 = abs(cap.second - next.second) / 15.0;
				S = std::max(S1, S2);
				if ((S >= (eps / 32.0)) && (S <= eps))
				{
					if ((xn + h) > xmax)
					{
						while ((!((xn > (xmax - prec)) && (xn <= xmax))) && ((xn + h)>xmax))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						result = RK4(xn, result.first, result.second, h);
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, ++dec);
					}
					else
					{
						result = next;
						xn += h;
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, dec);
					}

				}
				else if (S < (eps / 32.0))
				{
					if ((xn + h) > xmax)
					{
						while ((!((xn > (xmax - prec)) && (xn <= xmax))) && ((xn + h) > xmax))
						{
							//((xn > (xmax - prec)) && (xn < xmax))
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						result = RK4(xn, result.first, result.second, h);
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, ++dec);
					}
					else {
						result = next;
						xn += h;
						steps.insert(steps.begin() + i + 1, h);
						h *= 2.0;
						hinc.insert(hinc.begin() + i + 2, ++inc);
						hdec.insert(hdec.begin() + i + 2, dec);
					}
				}
				else if (S > eps)
				{
					if ((xn + h) > xmax)
					{
						while ((!((xn > (xmax - prec)) && (xn <= xmax))) && ((xn + h) > xmax))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						result = RK4(xn, result.first, result.second, h);
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, ++dec);
					}
					else {
						h = h / 2.0;
						result = RK4(xn, result.first, result.second, h);
						xn += h;
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, ++dec);
					}
				}
				i++;
				S *= 16.0;
				if (i == 1)
				{
					Smin = i;
					Smax = i;
					hmin = i;
					hmax = i;
				}
				else if (S < ss[Smin])
					Smin = i;
				else if (S > ss[Smax])
					Smax = i;
				if (h < steps[hmin])
					hmin = i;
				if (h > steps[hmax])
					hmax = i;
				reswcap.insert(reswcap.begin() + i, cap);
				ss.insert(ss.begin() + i, S);
				arg.insert(arg.begin() + i, xn);
				ures.insert(ures.begin() + i, result.first);
				vres.insert(vres.begin() + i, result.second);
			}
		}
		N = i;
		Xn = arg[i];
		Un = ures[i];
		Vn = vres[i];
		inc = hinc[i];
		dec = hdec[i];
		std::pair<std::vector<double>, std::vector<double>> res;
		res.first = ures;
		res.second = vres;
		return res;
	}
	friend std::ostream & operator<<(std::ostream &out, task2 &vc)
	{
		if (vc.ures.empty() || vc.ures.empty())
			out << "There are no calculated results yet.";
		else {
						out << "Основная задача 2" << std::endl;
			out << std::setw(5) << "n" << std::setw(12) << "h n-1" << std::setw(15) << "x" << std::setw(15) << "un" << std::setw(15) << "vn" << std::setw(15) << "u^" << std::setw(15) << "v^" << std::setw(15) << "S*" << std::setw(5) << "inc" << std::setw(4) << "dec" << std::endl;
			out << "--------------------------------------------------------------------------------------------------------------------" << std::endl;

			for (int i = 0; i < vc.ures.size(); i++)
			{
				out << std::setw(5) << i << std::setw(12) << vc.steps[i] << std::setw(15) << vc.arg[i] << std::setw(15) << vc.ures[i] << std::setw(15) << vc.vres[i] << std::setw(15) << vc.reswcap[i].first << std::setw(15) << vc.reswcap[i].second << std::setw(15) << vc.ss[i] << std::setw(5) << vc.hinc[i] << std::setw(4) << vc.hdec[i] << std::endl;
			}
		}
		return out;
	}

	void help()
	{
		std::cout << std::setw(50) << "Справка" << std::endl;
		std::cout << "Метод Рунге-Кутта порядка p = 4" << std::endl;
		std::cout << "Количество шагов n = " << N << std::endl;
		std::cout << "xmax - xn = " << (xmax - Xn) << std::endl;
		std::cout << "max |S*| = " << ss[Smax] << "  при х = " << arg[Smax] << std::endl;
		std::cout << "min |S*| = " << ss[Smin] << "  при х = " << arg[Smin] << std::endl;
		std::cout << "Общее число увеличений шага = " << inc << std::endl;
		std::cout << "Общее число уменьшений шага = " << dec << std::endl;
		std::cout << "Максимальный шаг h = " << steps[hmax] << "  при х = " << arg[hmax] << std::endl;
		std::cout << "Минимальный шаг h = " << steps[hmin] << "  при х = " << arg[hmin] << std::endl;
	}

	/*
	double reth(int count)
	{
		return steps[count];
	}

	double retx(int count)
	{
		return arg[count];
	}

	double retun(int count)
	{
		return ures[count];
	}

	double retvn(int count)
	{
		return vres[count];
	}

	double retun_with_cap(int count)
	{
		return reswcap[count].first;
	}

	double retvn_with_cap(int count)
	{
		return reswcap[count].second;
	}

	double rets(int count)
	{
		return ss[count];
	}

	int rethinc(int count)
	{
		return hinc[count];
	}

	int rethdec(int count)
	{
		return hdec[count];
	}

	QVector<double> retGraphIn()
	{

		return QVector<double>::fromStdVector(res);
	}

	QVector<double> retGraphI()
	{

		return QVector<double>::fromStdVector(exres);
	}

	QVector<double> retGraphX()
	{
		return QVector<double>::fromStdVector(arg);
	}

	double minX(void)
	{
		return arg[0];
	}

	double maxX(void)
	{
		return arg[arg.size() - 1];
	}

	double minI(void)
	{
		double min_val = exres[0];

		for (unsigned int i = 1; i < exres.size(); i++)
		{
			if (res[i] < min_val)
				min_val = res[i];
			if (exres[i] < min_val)
				min_val = exres[i];
		}
		return min_val;
	}

	double maxI(void)
	{
		double max_val = arg[0];

		for (unsigned int i = 1; i < arg.size(); i++)
		{
			if (res[i] > max_val)
				max_val = res[i];
			if (exres[i] > max_val)
				max_val = exres[i];
		}
		return max_val;
	}

	int retN(void) //индекс последнего шага
	{
		return n;
	}

	double retxN(void) //значение Xn (реальное значение)
	{
		return Xn;
	}

	double retxMax(void) //значение Xn - заданная правая граница
	{
		return xmax;
	}

	double retvN(void) //значение In
	{
		return Vn;
	}

	double retSmax(void) //max|S|
	{
		return ss[Smax] / 8;
	}

	double retSmin(void) //min|S|
	{
		return ss[Smin] / 8;
	}

	double retX_Smax(void) //max|S| при
	{
		return arg[Smax];
	}

	double retX_Smin(void) //min|S| при
	{
		return arg[Smin];
	}

	int retHinc(void) //всего увеличений шага
	{
		return inc;
	}

	int rethdec(void) //всего уменьшений шага
	{
		return dec;
	}

	double retHmax(void) //max hn
	{
		return steps[hmax];
	}

	double retHmin(void) //min hn
	{
		return steps[hmin];
	}

	double retX_hmax(void) //max hn при
	{
		return arg[hmax];
	}

	double retX_hmin(void) //min hn при
	{
		return arg[hmin];
	}
	*/
};