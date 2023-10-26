#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

// Solution: (-0.132, 0.705, -1.454)

int main()
{
	ifstream in("in2.txt");
	ofstream out("out2.txt");

	int n;
	in >> n;
	auto C = vector<vector<double>>(n, vector<double>(n, 0));
	auto D = vector<double>(n, 0);

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			in >> C[i][j];
		}
	}
	double normD = 0;
	for (size_t i = 0; i < n; i++)
	{
		in >> D[i];
		normD += abs(D[i]);
	}

	cout << "Input precision (digits after dot >= 1)" << endl;
	int precN;
	cin >> precN;

	double prec = pow(10., -precN);
	double alpha = 0.722082; // max eigen value of matrix C
	auto niter = int(floor(1. / log(alpha) * log(prec * (1. - alpha) / normD)) + 1.);

	cout << "computing system of LAE ..." << endl;

	auto res = vector<double>(n, 0);
	auto buf = vector<double>(n, 0);
	for (size_t iter = 0; iter < niter; iter++)
	{
		for (size_t i = 0; i < n; i++)
		{
			buf[i] = 0;
			for (size_t j = 0; j < n; j++)
			{
				buf[i] += C[i][j] * res[j];
			}
			buf[i] += D[i];
		}
		res = buf;
	}

	for (size_t i = 0; i < n; i++)
	{
		out << std::fixed << std::setprecision(precN + 2) << res[i] << endl;
	}

	out << "Iterations Done: " << niter << endl;

	cout << "Results in out2.txt" << endl;

	return 0;
}