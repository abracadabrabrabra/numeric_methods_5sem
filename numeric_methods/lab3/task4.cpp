#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

bool eq(const double& a, const double& b) {
    return fabs(a - b) < 1e-10;
}

void calculate_derivatives(const vector<double>& x, const vector<double>& y, const double& x_target) {

    for (int i = 0; i < x.size(); ++i) {
        if (eq(x_target, x[i])) {
            cout << "Target point is equal with node x[" << i << "] = " << x[i] << endl;
            if (i < x.size() - 1)
                cout << "Right-side derivative: " << (y[i + 1] - y[i]) / (x[i + 1] - x[i]) << endl;
            if (i > 0)
                cout << "Left-side  derivative: " << (y[i] - y[i - 1]) / (x[i] - x[i - 1]) << endl;
            if (i > 0 && i < x.size() - 1) {
                double central_deriv = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1]);
                cout << "Central derivative: " << central_deriv << endl;
            }

            if (i > 0 && i < x.size() - 1) {
                cout << "Second derivative: ";
                double h1 = x[i] - x[i - 1], h2 = x[i+1] - x[i];
                double second_deriv = 2 * (y[i - 1] / (h1 * (h1 + h2)) - y[i] / (h1 * h2) + y[i + 1] / (h2 * (h1 + h2)));
                cout << second_deriv << endl;
            }

            return;
        }
    }

    for (int i = 0; i < x.size() - 1; ++i) {
        if (x_target > x[i] && x_target < x[i + 1]) {
            cout << "Target point is in segment [" << x[i] << ", " << x[i + 1] << "]" << endl;
            cout << "First order derivative: " << (y[i + 1] - y[i]) / (x[i + 1] - x[i]) << endl;

            if (i > 0 && i < x.size() - 1) {
                cout << "Second derivative: ";
                    double h1 = x[i] - x[i - 1], h2 = x[i + 1] - x[i];
                double second_deriv = 2 * (y[i - 1] / (h1 * (h1 + h2)) - y[i] / (h1 * h2) + y[i + 1] / (h2 * (h1 + h2)));
                cout << second_deriv << endl;
            }

            return;
        }
    }

    throw invalid_argument("Target point out of range");
}

int main() {
    int n, i;
    cin >> n;
    vector<double> x(n);
    vector<double> y(n);
    double x_target;
    cin >> x_target;

    for (i = 0; i < n; ++i)
        cin >> x[i];
    for (i = 0; i < n; ++i)
        cin >> y[i];

    calculate_derivatives(x, y, x_target);

    return 0;
}