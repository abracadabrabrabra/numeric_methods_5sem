#include <iostream>
#include <cmath>
using namespace std;
using func = double (*)(double x);

double f(double x) { return log(x + 1.) - 2 * x + 0.5; }
double f_derivative(double x) { return (-2 * x - 1) / (x + 1); }
double F(double x) { return (log(x + 1.) + 0.5) / 2; }
double F_derivative(double x) { return 1. / (2 * x + 2); }
const double q = 0.5;

bool check_iter_convergence(func f, double x0) {
    return fabs(f(x0)) < 1;
}

pair<double, size_t> newton_method(func f, func df, const double& eps, const double& x0) {
    size_t cnt_iters = 0, max_iters = 10000;
    double x_new = x0, x_prev;
    do {
        x_prev = x_new;
        ++cnt_iters;
        if (cnt_iters == max_iters) {
            cout << "Too many iterations" << endl;
            break;
        }
        double deriv = df(x_prev);
        if (fabs(deriv) < 1e-10)
            throw runtime_error("Zero division error");
        x_new = x_prev - f(x_prev) / deriv;

    } while (fabs(x_new - x_prev) > eps);
    return { x_new, cnt_iters };
}

pair<double, size_t> iter_method(func f, func df, const double& eps, const double& x0) {
    size_t cnt_iters = 0, max_iters = 10000;
    double x_new = x0, x_prev, koeff = q / (1 - q);
    if (!check_iter_convergence(df, x0))
        cout << "Iteration method may not converge!" << endl;
    do {
        x_prev = x_new;
        ++cnt_iters;
        if (cnt_iters == max_iters) {
            cout << "Too many iterations" << endl;
            break;
        }
        x_new = f(x_prev);

    } while (fabs(x_new - x_prev) * koeff > eps);
    return { x_new, cnt_iters };
}

int main() {
    //x0 = 1
    double x0, eps;
    cin >> x0 >> eps;

    if (eps < 1e-12)
        throw runtime_error("Precision must be positive");

    auto res1 = newton_method(f, f_derivative, eps, x0),
        res2 = iter_method(F, F_derivative, eps, x0);
    cout << "For Newton method:" << endl;
    cout << "x = " << res1.first << ", was received on iteration " << res1.second << endl;

    cout << "For iteration method:" << endl;
    cout << "x = " << res2.first << ", was received on iteration " << res2.second << endl;

    return 0;
}