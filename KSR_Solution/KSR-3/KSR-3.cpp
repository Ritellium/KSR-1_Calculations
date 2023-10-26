#include <iostream>
#include <iomanip>
#include <vector>
#include <windows.h>

using namespace std;

// checked on l = 0.1 and l = 0.4
// Accuracy at t=0 is not granted 
/*
Task (var 5):
a = -1, b = 1,
x(t) - l * integrate [from -1 to 1] & [(t^2 - 1) * (s^3) * x(s) * ds] = t

Fredgholm integral equasion
*/

template <class Kernel>
double integrate(Kernel func, vector<double> x, double a, double b) {
    double step = (b-a) / (double(x.size()) - 1.);
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

double kernel_s(double x) {
    return x * x * x;
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
    double a = -1;
    double b = 1;
    double lim = 1./2.;

    cout << "Input lambda (multiplicator for integral) < " << lim << endl;
    double l;
    cin >> l;

    cout << "Input amount of dots to approximate function (n >= 2). good if n is odd" << endl;
    int n;
    cin >> n;

    cout << "Input precision (digits after dot >= 1)" << endl;
    int precN;
    cin >> precN;

    cout << "Calculating solution function..." << endl;

    double prec = pow(10., -precN);
    auto alpha = l;
    // Norm (õ1 - x0) = 1;
    auto niter = int(log(prec*(1-alpha)) / log(alpha)) + 1;

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
            buf[i] = t + l * (t * t - 1) * integrate(kernel_s, f_x, a, b);
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