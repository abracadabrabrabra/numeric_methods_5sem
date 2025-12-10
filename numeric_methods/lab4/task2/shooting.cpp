#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>

using namespace std;

double exact_solution(double x) {
    return exp(x) - 1;
}

double exact_derivative(double x) {
    return exp(x);
}

void system_ode(double x, const double* y, double* dydx) {
    double exp_x = exp(x);
    dydx[0] = y[1];  // y' = z
    dydx[1] = (2.0 * y[1] + exp_x * y[0]) / (exp_x + 1.0);  
}

vector<vector<double>> runge_kutta_4(
    function<void(double, const double*, double*)> system,
    double x0, double y0, double z0, double x_end, int n) {

    double h = (x_end - x0) / n;
    vector<double> x_vals(n + 1);
    vector<double> y_vals(n + 1);
    vector<double> z_vals(n + 1);

    x_vals[0] = x0;
    y_vals[0] = y0;
    z_vals[0] = z0;

    double y[2] = { y0, z0 };
    double k1[2], k2[2], k3[2], k4[2];
    double y_temp[2];

    for (int i = 0; i < n; i++) {
        double x = x_vals[i];

        system(x, y, k1);

        y_temp[0] = y[0] + 0.5 * h * k1[0];
        y_temp[1] = y[1] + 0.5 * h * k1[1];
        system(x + 0.5 * h, y_temp, k2);

        y_temp[0] = y[0] + 0.5 * h * k2[0];
        y_temp[1] = y[1] + 0.5 * h * k2[1];
        system(x + 0.5 * h, y_temp, k3);

        y_temp[0] = y[0] + h * k3[0];
        y_temp[1] = y[1] + h * k3[1];
        system(x + h, y_temp, k4);

        y[0] += (h / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
        y[1] += (h / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);

        x_vals[i + 1] = x + h;
        y_vals[i + 1] = y[0];
        z_vals[i + 1] = y[1];
    }

    return { x_vals, y_vals, z_vals };
}

double boundary_residual(double s, double x_end, int n) {
    auto result = runge_kutta_4(system_ode, 0.0, s, 1.0, x_end, n);

    double y_end = result[1].back();    // y(1)
    double z_end = result[2].back();    // y'(1)

    // y'(1) - y(1) should equal 1
    return (z_end - y_end) - 1.0;
}

// Bisection method for finding s = y(0)
double bisection_method(double a, double b, double tol, int max_iter, int n) {
    double fa = boundary_residual(a, 1.0, n);
    double fb = boundary_residual(b, 1.0, n);

    if (fa * fb > 0) {
        cerr << "Bisection method not applicable: f(a) and f(b) have same sign" << endl;
        return 0.5 * (a + b);
    }

    for (int i = 0; i < max_iter; i++) {
        double c = 0.5 * (a + b);
        double fc = boundary_residual(c, 1.0, n);

        if (abs(fc) < tol || abs(b - a) < tol) {
            cout << "Bisection method converged in " << i + 1 << " iterations" << endl;
            return c;
        }

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        }
        else {
            a = c;
            fa = fc;
        }
    }

    cout << "Bisection method did not converge in " << max_iter << " iterations" << endl;
    return 0.5 * (a + b);
}

double runge_romberg_error(const vector<double>& y_h, const vector<double>& y_h2, int p) {
    int n = y_h.size() - 1;
    double max_error = 0.0;

    for (int i = 0; i <= n; i++) {
        double error = abs(y_h[i] - y_h2[2 * i]) / (pow(2, p) - 1);
        if (error > max_error) {
            max_error = error;
        }
    }

    return max_error;
}

int main() {
    cout << fixed << setprecision(10);

    // Parameters
    double x0 = 0.0;
    double x_end = 1.0;
    int n1 = 20;      
    int n2 = 40;    

    cout << "\n=== Searching for proper interval ===" << endl;
    double step = 0.5;
    double a_test = -5.0;  
    double b_test = 5.0;

    double fa = boundary_residual(a_test, 1.0, n1);
    cout << "f(" << a_test << ") = " << fa << endl;

    for (double test = a_test + step; test <= b_test; test += step) {
        double fb = boundary_residual(test, 1.0, n1);
        cout << "f(" << test << ") = " << fb << endl;

        if (fa * fb < 0) {
            a_test = test - step;
            b_test = test;
            cout << "\nFound interval with sign change: ["
                << a_test << ", " << b_test << "]" << endl;
            break;
        }
        fa = fb;
    }

    double y0 = bisection_method(a_test, b_test, 1e-12, 100, n1);
        ;
    cout << "Found y(0) = " << y0 << endl;
    cout << "Exact y(0) = " << exact_solution(0.0) << endl;
    cout << "Error in y(0): " << abs(y0 - exact_solution(0.0)) << endl << endl;

    auto result1 = runge_kutta_4(system_ode, x0, y0, 1.0, x_end, n1);
    vector<double> x1 = result1[0];
    vector<double> y1 = result1[1];
    vector<double> z1 = result1[2];

    auto result2 = runge_kutta_4(system_ode, x0, y0, 1.0, x_end, n2);
    vector<double> x2 = result2[0];
    vector<double> y2 = result2[1];
    vector<double> z2 = result2[2];

    // Error estimation using Runge-Romberg method (p = 4 for RK4)
    double rr_error = runge_romberg_error(y1, y2, 4);

    // Compute actual error compared to exact solution
    double max_actual_error = 0.0;
    double max_relative_error = 0.0;
    int max_error_index = 0;

    cout << "=== Results ===" << endl;
    cout << setw(10) << "x" << setw(20) << "y_num" << setw(20) << "y_exact"
        << setw(20) << "Error" << setw(20) << "Rel. Error" << endl;
    cout << string(90, '-') << endl;

    for (int i = 0; i <= n1; i++) {
        double x = x1[i];
        double y_num = y1[i];
        double y_exact = exact_solution(x);
        double error = abs(y_num - y_exact);
        double rel_error = error / abs(y_exact);

        if (error > max_actual_error) {
            max_actual_error = error;
            max_relative_error = rel_error;
            max_error_index = i;
        }

        if (i % 2 == 0) {  // Print every second point
            cout << setw(10) << x
                << setw(20) << y_num
                << setw(20) << y_exact
                << setw(20) << error
                << setw(20) << rel_error << endl;
        }
    }

    cout << endl;
    cout << "=== Error analysis ===" << endl;
    cout << "Maximum absolute error: " << max_actual_error << endl;
    cout << "Maximum relative error: " << max_relative_error << endl;
    cout << "Runge-Romberg error estimate: " << rr_error << endl;
    cout << "x where max error occurs: " << x1[max_error_index] << endl;

    // Check boundary conditions
    cout << endl << "=== Boundary conditions check ===" << endl;
    cout << "Left boundary - y'(0):" << endl;
    cout << "  Numerical: " << z1[0] << endl;
    cout << "  Exact: " << exact_derivative(0.0) << endl;
    cout << "  Error: " << abs(z1[0] - exact_derivative(0.0)) << endl;

    cout << endl << "Right boundary - y'(1) - y(1):" << endl;
    double boundary_value = z1.back() - y1.back();
    cout << "  Numerical: " << boundary_value << endl;
    cout << "  Should be: 1" << endl;
    cout << "  Error: " << abs(boundary_value - 1.0) << endl;

    return 0;
}