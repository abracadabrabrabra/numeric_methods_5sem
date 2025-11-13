#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

double f(double x) { return tan(x) + x; }
constexpr double point = 3 * M_PI / 16;
const double f_point = f(point);

double gorner_scheme(const vector<double>& pol, const double val) {
    double res = 0;
    auto it = pol.crbegin();
    while (it != pol.crend()) {
        res = res * val + *it;
        ++it;
    }
    return res;
}

double delta(const vector<double>& pol) {
    return fabs(f_point - gorner_scheme(pol, point));
}

void print_polynom(const vector<double>& pol) {
    if (pol.empty())  throw invalid_argument("Empty polynom");
    if (pol[0])   cout << pol[0] << " + ";
    for (size_t i = 0; i < pol.size(); ++i) {
        if (!pol[i])  continue;
        cout << " (" << pol[i] << ")*";
        if (i == 1)   cout << "x";
        else    cout << "x^" << i;
        if (i != pol.size() - 1)  cout << " +";
    }
    cout << endl;
}

vector<double> basis_lagrange(const size_t& i, const vector<double>& x) {
    size_t n = x.size();
    vector<double> coeffs(1, 1.);

    for (size_t j = 0; j < n; ++j) {
        if (j == i) continue;
        double denom = x[i] - x[j];
        vector<double> temp_coeffs(coeffs.size() + 1, 0.);
        for (size_t k = 0; k < coeffs.size(); ++k) {
            temp_coeffs[k] += coeffs[k] * (-x[j] / denom);
            temp_coeffs[k + 1] += coeffs[k] / denom;
        }
        coeffs = move(temp_coeffs);
    }

    return coeffs;
}

vector<double> lagrange_coeffs(const vector<double>& x, const vector<double>& y) {
    if (x.size() != y.size())   throw invalid_argument("Different size of vectors");
    size_t n = x.size();
    vector<double> res(n);

    for (size_t i = 0; i < n; ++i) {
        vector<double> basis_coeffs = basis_lagrange(i, x);
        for (size_t j = 0; j < min(n, basis_coeffs.size()); ++j)
            res[j] += y[i] * basis_coeffs[j];
    }

    return res;
}

void polynom_lagrange(const vector<double>& x, const vector<double>& y) {
    vector<double> res = lagrange_coeffs(x, y);
    cout << " Lagrange polynom:" << endl;
    print_polynom(res);
    cout << " Eps as delta funcs: " << delta(res) << endl;
}

vector<double> divided_diffs(const vector<double>& x, const vector<double>& y) {
    if (x.size() != y.size())   throw invalid_argument("Different size of vectors");
    size_t n = x.size();
    vector<vector<double>> res(n);
    res[0].resize(n);
    vector<double> dd(n);

    for (size_t i = 0; i < n; ++i)
        res[0][i] = y[i];
    dd[0] = res[0][0];

    for (size_t i = 1; i < n; ++i) {
        res[i].resize(n - i);
        for (size_t j = 0; j < n - i; ++j)
            res[i][j] = (res[i - 1][j] - res[i - 1][j + 1]) / (x[j] - x[j + i]);
        dd[i] = res[i][0];
    }

    return dd;
}

double aposteriory_error(const vector<double>& x, double& next_dd) {
    double omega = 1.;
    for (size_t i = 0; i < x.size() - 1; ++i)
        omega *= (point - x[i]);

    return fabs(next_dd * omega);
}

vector<double> newton_coeffs(const vector<double>& x, const vector<double>& y) {
    vector<double> dd = divided_diffs(x, y);
    size_t n = dd.size();
    cout << " Eps as aposteriory approximation for Newton polynom: " << aposteriory_error(x, dd[n - 1]) << endl;
    vector<double> standard_coeffs(n, 0.), cur_poly = { 1. };
    standard_coeffs[0] = dd[0];

    for (size_t i = 1; i < n; ++i) {
        vector<double> new_poly(cur_poly.size() + 1);
        for (size_t j = 0; j < cur_poly.size(); ++j) {
            new_poly[j] += cur_poly[j] * (-x[i - 1]);
            new_poly[j + 1] += cur_poly[j];
        }
        cur_poly = move(new_poly);

        for (size_t j = 0; j < min(cur_poly.size(), standard_coeffs.size()); ++j) {
            standard_coeffs[j] += dd[i] * cur_poly[j];
        }
    }

    return standard_coeffs;
}

void polynom_newton(const vector<double>& x, const vector<double>& y) {
    vector<double> res = newton_coeffs(x, y);
    cout << " Newton polynom:" << endl;
    print_polynom(res);
    cout << " Eps as delta funcs: " << delta(res) << endl;
}



int main() {

    vector<double> x_i{ 0, M_PI / 8, M_PI / 4, 3 * M_PI / 8 };
    vector<double> y_i{ f(x_i[0]), f(x_i[1]), f(x_i[2]), f(x_i[3]) };

    cout << "Variant A:" << endl;
    polynom_lagrange(x_i, y_i);
    x_i[2] = M_PI / 3;
    y_i[2] = f(x_i[2]);
    cout << "Variant B:" << endl;
    polynom_newton(x_i, y_i);

    return 0;
}
