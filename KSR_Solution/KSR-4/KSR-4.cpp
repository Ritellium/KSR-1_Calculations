#include <iostream>
#include <iomanip>
#include <vector>
#include <windows.h>

using namespace std;

// 
// 
/*
Task (var 5):
a = -1, b = 1,
x(t) - integrate [from 0 to t] & [(t - s) * x(s) * ds] = 1 + t

Fredgholm integral equasion
*/

template <class Kernel>
double integrate(Kernel func, vector<double> x, double a, double b) {
    double step = (b - a) / (double(x.size()) - 1.);
    double point = a;
    double res = func(a) * x[0];
    size_t iter = 1;
    for (iter = 1; iter < x.size() - 1; iter++)
    {
        if (iter / 2 == 1)
            res += 4. * func(point) * x[iter];
        else
            res += 2. * func(point) * x[iter];
        point += step;
    }
    res += func(b) * x[iter];
    res *= step / 3.;
    return res;
}

double kernel_s_1(double x) {
    return 1.;
}

double kernel_s_2(double x) {
    return x;
}

void output_list(vector<double> v) {
    cout << "{";
    for (size_t i = 0; i < v.size(); i++)
    {
        if (i < v.size() - 1)
        {
            cout << v[i] << ", ";
        }
        else {
            cout << v[i];
        }
    }
    cout << "}" << endl;
}

int main()
{
    double a = 0;
    double b = 1;

    cout << "Input amount of dots to approximate function (n >= 2). good if n is odd" << endl;
    int n;
    cin >> n;

    cout << "Input amount of iterations: (i >= 2)" << endl;
    int niter;
    cin >> niter;

    cout << "Calculating solution function..." << endl;

    auto grid = vector<double>(n, 0);
    double point = a;
    double step = (b - a) / double(n - 1);
    for (size_t i = 0; i < n; i++)
    {
        grid[i] = point;
        point += step;
    }

    auto f_x = vector<double>(n, 0);
    auto buf = vector<double>(n, 0);
    for (size_t iter = 0; iter < niter; iter++)
    {
        for (size_t i = 0; i < n; i++)
        {
            auto t = grid[i];
            buf[i] = 1 + t + t * integrate(kernel_s_1, f_x, a, t) - integrate(kernel_s_2, f_x, a, t);
        }
        f_x = buf;
    }

    cout << "iterations done: " << niter << endl;
    cout << "function: " << endl;
    for (size_t i = 0; i < n; i++)
    {
        cout << "t = " << grid[i] << ", x(t) = " << f_x[i] << endl;
    }

    cout << endl;
    cout << "datalists for wolfram: " << endl;
    output_list(grid);
    output_list(f_x);

    return 0;
}