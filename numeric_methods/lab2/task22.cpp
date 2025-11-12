#include <iostream>
#include <cmath>
using namespace std;

double f1(double x, double y) { return x * x / 4. + y * y - 1.; }
double f2(double x, double y) { return 2. * y - exp(x) - x; }
double fi1(double x, double y) { return log(2 * y - x); }
double fi2(double x, double y) { return sqrt(1. - x * x / 4.); }

constexpr double q = 0.29;

pair<double, double> iter_method(double x0, double y0, double eps, size_t max_iters = 10000) {
    cout << "Simple iteration method:" << endl;
    size_t cnt_iters = 0;
    double x_new = x0, x_prev, y_new = y0, y_prev, koeff = q / (1 - q);

    do {
        x_prev = x_new;
        y_prev = y_new;
        ++cnt_iters;
        if (cnt_iters == max_iters) {
            cout << "Too many iterations" << endl;
            break;
        }

        x_new = fi1(x_prev, y_prev);
        y_new = fi2(x_prev, y_prev);

        if (!isfinite(x_new) || !isfinite(y_new)) {
            cout << "Nan or infinity detected at iteration " << cnt_iters << endl;
            cout << "x_prev = " << x_prev << ", y_prev = " << y_prev << endl;
            break;
        }
        if (fabs(x_new) > 1.5 ||  y_new > 1.2  || y_new < (x_new + 3.) / 4.) {
            cout << "Approximation on iteration " << cnt_iters << "isn't located in area G" << endl;
            cout << "x_prev = " << x_prev << ", y_prev = " << y_prev << endl;
            break;
        }

        cout << "Iter " << cnt_iters << ": x = " << x_new << ", y = " << y_new << endl;
    } while (koeff * (fabs(x_new - x_prev) + fabs(y_new - y_prev)) > eps);

    cout << "Method was stopped on iteration " << cnt_iters << endl;
    return { x_new, y_new };
}

int main() {
    //0.5 0.9
    double x0, y0, eps;
    cin >> x0 >> y0 >> eps;
    if (eps < 1e-16)
        throw runtime_error("Precision must be positive");

    auto res = iter_method(x0, y0, eps);
    cout << "Solution: x = " << res.first << ", y = " << res.second << endl;
    cout << "Back substitution:" << endl;
    cout << "f1 = " << f1(res.first, res.second) << endl;
    cout << "f2 = " << f2(res.first, res.second) << endl;

    return 0;
}