#pragma once
#include <iostream>
#include <cmath>
#include <vector>

// calculation stops when n steps are done

class task1 // variable current
{
private:
	double v0, x0;
	double h, eps, xmax, prec;
	int n; // number of steps
	int N = 0; // total steps
	double Xn, Vn;
	int inc = 0;
	int dec = 0; // total inc dec
	int Smin = 0, Smax = 0; // number of string with Smin, Smax
	int hmax, hmin; //number of string with hmax, hmin
	std::vector<double> arg; //x
	std::vector<double> res; //v
	std::vector<double> reshalf; //v with cap
	std::vector<double> steps; //h
	std::vector<double> ss; //S
	std::vector<int> hinc;  // total step increases
	std::vector<int> hdec;  // total step decreases
	double func(double x, double u)
	{
		return (((pow(x, 3) + 1) / (pow(x, 5) + 1))*pow(u, 2) + u - pow(u, 3)*sin(10 * x));
	}
public:
	task1(double _x0, double _v0, double _h, int _n, double _eps, double _xmax, double _prec)
	{
		x0 = _x0;
		v0 = _v0;
		h = _h;
		n = _n;
		eps = _eps;
		xmax = _xmax;
		prec = _prec;
	}
	void reset(double _x0, double _v0, double _h, int _n, double _eps, double _xmax, double _prec)
	{
		arg.clear();
		res.clear();
		reshalf.clear();
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
		v0 = _v0;
		h = _h;
		n = _n;
		eps = _eps;
		xmax = _xmax;
		prec = _prec;
	}
	double RK4(double xn, double vn, double h)
	{
		double k1 = func(xn, vn);
		double k2 = func(xn + h / 2.0, vn + h * k1 / 2.0);
		double k3 = func(xn + h / 2.0, vn + h * k2 / 2.0);
		double k4 = func(xn + h, vn + h * k3);
		vn = vn + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
		return vn;
	}
	std::vector<double> calculate()
	{
		arg.push_back(x0);
		res.push_back(v0);
		hinc.push_back(0);
		ss.push_back(0.0);
		hdec.push_back(0);
		steps.push_back(0.0);
		reshalf.push_back(0.0);
		double xn = x0;
		double vn = v0;
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
					vn = RK4(xn, vn, h);
					xn += h;
				}
				else
				{
					vn = RK4(xn, vn, h);
					xn += h;
				}
				arg.insert(arg.begin() + i + 1, xn);
				res.insert(res.begin() + i + 1, vn);
				steps.insert(steps.begin() + i + 1, h);
				hinc.insert(hinc.begin() + i + 1, inc);
				hdec.insert(hdec.begin() + i + 1, dec);
				reshalf.insert(reshalf.begin() + i + 1, 0.0);
				ss.insert(ss.begin() + i + 1, 0.0);
				i++;
			}
		}
		N = i;
		Xn = arg[i];
		Vn = res[i];
		if (i == n)
			hmin = 1;
		else hmin = i;
		hmax = 1;
		return res;
	}
	std::vector<double> calculate_w_error()
	{
		ss.push_back(0.0);
		arg.push_back(x0);
		res.push_back(v0);
		steps.push_back(0.0);
		reshalf.push_back(0.0);
		hinc.push_back(0);
		hdec.push_back(0);
		hinc.push_back(0);
		hdec.push_back(0);
		double xn = x0;
		double vn = v0;
		double xhalf = x0;
		double vhalf, reswcap, vnext;
		double S, Sabs;
		int i = 0;
		while (i < n)
		{
			if ((xn > (xmax - prec)) && (xn <= xmax))
			{
				break;
			}
			else
			{
				vhalf = RK4(xn, vn, h / 2.0);
				xhalf = xn + h / 2.0;
				reswcap = RK4(xhalf, vhalf, h / 2.0);
				vnext = RK4(xn, vn, h);
				S = (reswcap - vnext) / 15.0;
				Sabs = abs(S);
				if ((Sabs >= (eps / 32.0)) && (Sabs <= eps))
				{
					if ((xn + h) > xmax)
					{
						while ((!((xn > (xmax - prec)) && (xn <= xmax))) && ((xn + h) > xmax))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						vn = RK4(xn, vn, h);
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, ++dec);
					}
					else
					{
						vn = vnext;
						xn += h;
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, dec);
					}

				}
				else if (Sabs < (eps / 32.0))
				{
					if ((xn + h) > xmax)
					{
						while ((!((xn > (xmax - prec)) && (xn <= xmax))) && ((xn + h) > xmax))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						vn = RK4(xn, vn, h);
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, ++dec);
					}
					else {
						vn = vnext;
						steps.insert(steps.begin() + i + 1, h);
						xn += h;
						h *= 2.0;
						hinc.insert(hinc.begin() + i + 1, ++inc);
						hdec.insert(hdec.begin() + i + 1, dec);
					}
				}
				else if (Sabs > eps)
				{
					if ((xn + h) > xmax)
					{
						while ((!((xn > (xmax - prec)) && (xn <= xmax))) && ((xn + h) > xmax))
						{
							h /= 2.0;
						}
						xn += h;
						xhalf = xn + h / 2.0;
						vn = RK4(xn, vn, h);
						steps.insert(steps.begin() + i + 1, h);
						hinc.insert(hinc.begin() + i + 1, inc);
						hdec.insert(hdec.begin() + i + 1, ++dec);
					}
					else {
						h = h / 2.0;
						vn = RK4(xn, vn, h);
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
				else if (Sabs < ss[Smin])
					Smin = i;
				else if (Sabs > ss[Smax])
					Smax = i;
				if (h < steps[hmin])
					hmin = i;
				if (h > steps[hmax])
					hmax = i;
				reshalf.insert(reshalf.begin() + i, reswcap);
				ss.insert(ss.begin() + i, S);
				arg.insert(arg.begin() + i, xn);
				res.insert(res.begin() + i, vn);
			}
		}
		N = i;
		Xn = arg[i];
		Vn = res[i];
		inc = hinc[i];
		dec = hdec[i];
		return res;
	}
	friend std::ostream & operator<<(std::ostream &out, task1 &vc)
	{
		if (vc.res.empty())
			out << "There are no calculated results yet.";
		else {
			out << "Основная задача 1" << std::endl;
			out << std::setw(5) << "n" << std::setw(12) << "h n-1" << std::setw(10) << "x" << std::setw(10) << "vn" << std::setw(10) << "v^" << std::setw(15) << "S*" << std::setw(5) << "inc" << std::setw(4) << "dec" << std::endl;
			out << "-----------------------------------------------------------------------" << std::endl;

			for (int i = 0; i < vc.res.size(); i++)
			{
				out << std::setw(5) << i << std::setw(12) << vc.steps[i] << std::setw(10) << vc.arg[i] << std::setw(10) << vc.res[i] << std::setw(10) << vc.reshalf[i] << std::setw(15) << vc.ss[i] << std::setw(5) << vc.hinc[i] << std::setw(4) << vc.hdec[i] << std::endl;
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

	double retin(int count)
	{
		return res[count];
	}

	double retin_with_cap(int count)
	{
		return reshalf[count];
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