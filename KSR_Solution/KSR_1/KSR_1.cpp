#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

double functionG(double x) {
	return x + atan(x) + x / (1 + x * x) - 3;
}

int main()
{
	// x + ArcTan[x] + x/(1 + x^2) - 3
	// 16 digits solution -  1.5471210825582697446
	double k1 = 1;
	double k2 = 3;

	cout << "Input precision (digits after dot >= 1)" << endl;
	int precN;
	cin >> precN;
	cout << "Input start value (any number)" << endl;
	double x0;
	cin >> x0;
	
	cout << "computing x + ArcTan[x] + x/(1 + x^2) - 3 ..." << endl;

	double prec = pow(10., -precN);
	double alpha = 1. - k1 / k2;
	auto n = int ( floor(1. / log(alpha) * log(prec * (1. - alpha) / abs(1. / k2 * functionG(0)))) + 1. );
	double res = x0;
	for (int i = 0; i < n; i++)
	{
		res = res - 1. / k2 * functionG(res);
	}

	cout << std::fixed << std::setprecision(precN + 3) <<  "Result: " << res << endl;
	cout << "Iterations done: " << n << endl;

	return 0;
}