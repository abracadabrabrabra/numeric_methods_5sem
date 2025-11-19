#include <iostream>
#include <vector>
#include <iomanip>
#include <array>
#include <cmath>
using namespace std;
using callback = double (*)(double);

struct integral_vals {
    double rect;
    double trapeze;
    double simpson;
};

double calc_runge_romberg(double I_h, double I_h2, int order) {
    double denom = pow(2, order) - 1;
    return (I_h2 - I_h) / denom;
}

void print_itegrals(const integral_vals& vals) {
    cout << "Rectangle method: " << vals.rect << endl;
    cout << "Trapeze method: " << vals.trapeze << endl;
    cout << "Simpson's method: " << vals.simpson << endl;
}

double f(double x) { return 1 / (x * x * x * x + 16); }
double func_ex(double x) { return x / ((3 * x + 4) * (3 * x + 4)); }

void calculate_integral(callback func, integral_vals& dst, double start, double fin, double step) {
    if (fin < start || fabs(fin - start) < 1e-10 || step < 1e-10)
        throw invalid_argument("Incorrect range or step <= 0");

    int segments = static_cast<int>((fin - start) / step);
    if (segments % 2)
        throw invalid_argument("Odd number of segments for Simpson's method");

    dst.rect = 0.;
    dst.trapeze = 0.;
    dst.simpson = 0.;

    for (int i = 0; i < segments; ++i) {
        double x_mid = start + (i + 0.5) * step;
        dst.rect += func(x_mid);
    }
    dst.rect *= step;

    dst.trapeze = (func(start) + func(fin)) * 0.5;
    for (int i = 1; i < segments; ++i) {
        double x_i = start + i * step;
        dst.trapeze += func(x_i);
    }
    dst.trapeze *= step;

    dst.simpson = func(start) + func(fin);
    for (int i = 1; i < segments; ++i) {
        double x_i = start + i * step;
        if (i % 2 == 1) 
            dst.simpson += 4 * func(x_i);
        else 
            dst.simpson += 2 * func(x_i);
    }
    dst.simpson *= step / 3.0;
}

void calculate_err(array<double, 3>& dst, const integral_vals& vh1, const integral_vals& vh2) {
    dst[0] = calc_runge_romberg(vh1.rect, vh2.rect, 2);
    dst[1] = calc_runge_romberg(vh1.trapeze, vh2.trapeze, 2);
    dst[2] = calc_runge_romberg(vh1.simpson, vh2.simpson, 4);
}


int main() {

    integral_vals vh1, vh2;

    cout << fixed << setprecision(6);
    cout << "Calculating integral:\n\nStep 0.5:" << endl;
    calculate_integral(f, vh1, 0., 2., 0.5);
    print_itegrals(vh1);

    cout << "\nStep 0.25:" << endl;
    calculate_integral(f, vh2, 0., 2., 0.25);
    print_itegrals(vh2);

    array<double, 3> err;
    calculate_err(err, vh1, vh2);
    cout << "\nRunge_Romberg method:" << endl;
    cout << "Rectangle method: error = " << fabs(err[0]) << ", refined value = " << vh2.rect + err[0] << endl;
    cout << "Trapeze method: error = " << fabs(err[1]) << ", refined value = " << vh2.trapeze + err[1] << endl;
    cout << "Simpson's method: error = " << fabs(err[2]) << ", refined value = " << vh2.simpson + err[2] << endl;

    return 0;
}
