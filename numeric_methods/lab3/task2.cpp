#include <iostream>
#include <vector>
#include <cstdint>
using namespace std;

int main() {

    int n, i;
    cin >> n;
    vector<double> x(n);
    double x_target, x_prev, x_cur;
    cin >> x_target;

    vector<double> h(n-1);
    vector<double> f(n);

    cin >> x_prev;
    x[0] = x_prev;
    for (i = 1; i < n; ++i) {
        cin >> x_cur;
        x[i] = x_cur;
        h[i - 1] = x_cur - x_prev;
        x_prev = x_cur;
    }
    for (i = 0; i < n; ++i)
        cin >> f[i];

    vector<double> c(n);
    c[0] = 0.;
    c[n - 1] = 0.;
    vector<double> p(n), q(n);

    //tridiagonal matrix algo
    //direct move
    double denominator = 2 * (h[0] + h[1]);
    p[1] = -h[1] / denominator;
    q[1] = 3 * ((f[2] - f[1]) / h[1] - (f[1] - f[0]) / h[0]) / denominator;

    for (i = 1; i < n - 3; ++i) {
        denominator = 2 * (h[i] + h[i+1]) + h[i] * p[i];
        p[i+1] = -h[i+1] / denominator;
        q[i+1] = (3 * ((f[i+2] - f[i+1]) / h[i+1] - (f[i+1] - f[i]) / h[i]) - h[i] * q[i]) / denominator;
    }
    

    denominator = 2 * (h[n-3] + h[n-2]) + h[n-3] * p[n-3];
    p[n-2] = -h[n-2] / denominator;
    q[n-2] = (3 * ((f[n-1] - f[n-2]) / h[n-2] - (f[n-2] - f[n-3]) / h[n-3]) - h[n-3] * q[n-3]) / denominator;

    //reverse move
    c[n-2] = q[n-2];
    for (i = n - 3; i > 0; --i) 
        c[i] = p[i] * c[i+1] + q[i];

    vector<double> a(n - 1), b(n - 1), d(n - 1);
    for (i = 0; i < n - 1; ++i) {
        a[i] = f[i];
        b[i] = (f[i + 1] - f[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }

    cout << "Cube spline coeffs:" << endl;
    for (i = 0; i < n - 1; ++i) {
        cout << "Segment [" << x[i] << ", " << x[i + 1] << "]:" << endl;
        cout << "  a" << (int)i << " = " << a[i] << endl;
        cout << "  b" << (int)i << " = " << b[i] << endl;
        cout << "  c" << (int)i << " = " << c[i] << endl;
        cout << "  d" << (int)i << " = " << d[i] << endl;
        cout << "  S" << (int)i << "(x) = " << a[i] << " + " << b[i] << "(x - " << x[i] << ") + "
            << c[i] << "(x - " << x[i] << ")^2 + " << d[i] << "(x - " << x[i] << ")^3" << endl;
    }

    cout << "X_target = " << x_target << ", it lies in segment [0.9, 1.8]" << endl;
    cout << "Value in target point" << ": ";
    double dx = x_target - 0.9;
    cout << a[1] + b[1] * dx + c[1] * dx * dx + d[1] * dx * dx * dx << endl;

    return 0;
}