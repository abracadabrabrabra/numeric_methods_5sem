#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

void process_coeffs_p1(const vector<double>& x, const vector<double>& y, vector<double>& coeffs) {
    coeffs[0] = x.size();
    for (int i = 0; i < x.size(); ++i) {
        coeffs[1] += x[i];
        coeffs[2] += y[i];
        coeffs[3] += x[i] * x[i];
        coeffs[4] += y[i] * x[i];
    }
}

pair<double, double> pow_1(const vector<double>& coeffs) {
    pair<double, double> res;
    double numerator = (coeffs[4] - coeffs[1] * coeffs[2] / coeffs[0]);
    double denominator = (coeffs[3] - coeffs[1] * coeffs[1] / coeffs[0]);
    res.second = numerator / denominator;  
    res.first = (coeffs[2] - coeffs[1] * res.second) / coeffs[0];  
    return res;
}

struct quadro_coeffs {
    double a0;
    double a1;
    double a2;
};

void process_coeffs_p2(const vector<double>& x, const vector<double>& y,
    double& n, double& sum_x, double& sum_y,
    double& sum_x2, double& sum_x3, double& sum_x4,
    double& sum_xy, double& sum_x2y) {
    n = x.size();
    sum_x = sum_y = sum_x2 = sum_x3 = sum_x4 = sum_xy = sum_x2y = 0.0;

    for (int i = 0; i < n; ++i) {
        double xi = x[i];
        double xi2 = xi * xi;
        double xi3 = xi2 * xi;
        double xi4 = xi3 * xi;

        sum_x += xi;
        sum_y += y[i];
        sum_x2 += xi2;
        sum_x3 += xi3;
        sum_x4 += xi4;
        sum_xy += xi * y[i];
        sum_x2y += xi2 * y[i];
    }
}

quadro_coeffs pow2(const vector<double>& x, const vector<double>& y) {
    double n, sum_x, sum_y, sum_x2, sum_x3, sum_x4, sum_xy, sum_x2y;
    process_coeffs_p2(x, y, n, sum_x, sum_y, sum_x2, sum_x3, sum_x4, sum_xy, sum_x2y);

    double A[3][3] = {
        {n,    sum_x,  sum_x2},
        {sum_x, sum_x2, sum_x3},
        {sum_x2, sum_x3, sum_x4}
    };
    double B[3] = { sum_y, sum_xy, sum_x2y };

    double detA = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    if (fabs(detA) < 1e-12) {
        cerr << "Singular matrix" << endl;
        return { 0, 0, 0 };
    }

    quadro_coeffs coeffs;
    //Kramer
    double A0[3][3] = {
        {B[0],    sum_x,  sum_x2},
        {B[1], sum_x2, sum_x3},
        {B[2], sum_x3, sum_x4}
    };
    double detA0 = A0[0][0] * (A0[1][1] * A0[2][2] - A0[1][2] * A0[2][1])
        - A0[0][1] * (A0[1][0] * A0[2][2] - A0[1][2] * A0[2][0])
        + A0[0][2] * (A0[1][0] * A0[2][1] - A0[1][1] * A0[2][0]);
    coeffs.a0 = detA0 / detA;

    double A1[3][3] = {
        {n,    B[0],  sum_x2},
        {sum_x, B[1], sum_x3},
        {sum_x2, B[2], sum_x4}
    };
    double detA1 = A1[0][0] * (A1[1][1] * A1[2][2] - A1[1][2] * A1[2][1])
        - A1[0][1] * (A1[1][0] * A1[2][2] - A1[1][2] * A1[2][0])
        + A1[0][2] * (A1[1][0] * A1[2][1] - A1[1][1] * A1[2][0]);
    coeffs.a1 = detA1 / detA;

    double A2[3][3] = {
        {n,    sum_x,  B[0]},
        {sum_x, sum_x2, B[1]},
        {sum_x2, sum_x3, B[2]}
    };
    double detA2 = A2[0][0] * (A2[1][1] * A2[2][2] - A2[1][2] * A2[2][1])
        - A2[0][1] * (A2[1][0] * A2[2][2] - A2[1][2] * A2[2][0])
        + A2[0][2] * (A2[1][0] * A2[2][1] - A2[1][1] * A2[2][0]);
    coeffs.a2 = detA2 / detA;

    return coeffs;
}

int main() {
    int n, i;
    cin >> n;
    vector<double> x(n);
    vector<double> y(n);
    vector<double> coeffs_p1(5, 0.);

    for (i = 0; i < n; ++i)
        cin >> x[i];
    for (i = 0; i < n; ++i)
        cin >> y[i];

    process_coeffs_p1(x, y, coeffs_p1);
    pair<double, double> c1 = pow_1(coeffs_p1);

    cout << "Polynom pow 1:" << endl;
    cout << fixed << setprecision(6);
    cout << "Result: a0 = " << c1.first << ", a1 = " << c1.second << endl;
    double F1 = 0.;
    for (i = 0; i < n; ++i) {
        double delta = fabs(y[i] - (c1.first + c1.second * x[i]));
        F1 += delta * delta;
    }
    cout << "F = " << F1 << endl;

    cout << "Polynom pow 2:" << endl;
    quadro_coeffs c2 = pow2(x, y);
    cout << fixed << setprecision(6);
    cout << "Result: a0 = " << c2.a0 << ", a1 = " << c2.a1 << ", a2 = " << c2.a2 << endl;
    double F2 = 0.;
    for (i = 0; i < n; ++i) {
        double delta = fabs(y[i] - (c2.a0 + c2.a1 * x[i] + c2.a2 * x[i] * x[i]));
        F2 += delta * delta;
    }
    cout << "F = " << F2 << endl;

    return 0;
}