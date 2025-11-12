#include <iostream>
#include <cmath>
using namespace std;
using func = double (*)(double x, double y);

double f1(double x, double y) { return x * x / 4. + y * y - 1.; }
double f2(double x, double y) { return 2. * y - exp(x) - x; }
double df1_dx(double x, double y) { return x / 2.; }
double df1_dy(double x, double y) { return 2. * y; }
double df2_dx(double x, double y) { return -exp(x) - 1.; }
double df2_dy(double x, double y) { return 2.; }

double det(double x, double y, func a, func b, func c, func d) {
    return a(x, y) * d(x, y) - b(x, y) * c(x, y);
}

pair<double, double> newton_method(double x0, double y0, double eps, size_t max_iters = 10000) {
    cout << "Newton method:" << endl;
    size_t cnt_iters = 0;
    double x_new = x0, x_prev, y_new = y0, y_prev;

    do {
        x_prev = x_new;
        y_prev = y_new;
        ++cnt_iters;
        if (cnt_iters == max_iters) {
            cout << "Too many iterations" << endl;
            break;
        }
        double det_J = det(x_prev, y_prev, df1_dx, df1_dy, df2_dx, df2_dy);
        if (fabs(det_J) < 1e-12)
            throw runtime_error("Singular Jacobian, zero division error");
        x_new = x_prev - det(x_prev, y_prev, f1, df1_dy, f2, df2_dy) / det_J;
        y_new = y_prev - det(x_prev, y_prev, df1_dx, f1, df2_dx, f2) / det_J;

    } while (fabs(x_new - x_prev) + fabs(y_new - y_prev) > eps);

    cout << "Method was stopped on iteration " << cnt_iters << endl;
    return { x_new, y_new };
}


int main() {
    //0.5 ; 1
    double x0, y0, eps;
    cin >> x0 >> y0 >> eps;

    if (eps < 1e-12)
        throw runtime_error("Precision must be positive");

    auto res = newton_method(x0, y0, eps);
    cout << "x1 = " << res.first << ", x2 = " << res.second << endl;

    cout << "Back substitution" << endl;
    cout << "f1 = " << f1(res.first, res.second) << "; f2 = " << f2(res.first, res.second) << endl;

    return 0;
}