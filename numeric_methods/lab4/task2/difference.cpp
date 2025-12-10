#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

double exact_solution(double x) {
    return exp(x) - 1;
}

double exact_derivative(double x) {
    return exp(x);
}

// Finite difference method implementation
vector<double> solve_bvp_fd(int n) {
    double a = 0.0, b = 1.0;
    double h = (b - a) / n;

    vector<double> x(n + 1);
    vector<double> y(n + 1);

    // Create grid
    for (int i = 0; i <= n; i++) {
        x[i] = a + i * h;
    }


    vector<double> A(n + 1), B(n + 1), C(n + 1), D(n + 1);

    // Interior points (i = 1, ..., n-1)
    for (int i = 1; i < n; i++) {
        double xi = x[i];
        double exp_xi = exp(xi);
        double denom = exp_xi + 1.0;

        double p = -2.0 / denom;   
        double q = -exp_xi / denom; 

        // Central difference approximations:
        A[i] = 1.0 / (h * h) - p / (2.0 * h);
        B[i] = -2.0 / (h * h) + q;
        C[i] = 1.0 / (h * h) + p / (2.0 * h);
        D[i] = 0.0;
    }

    A[0] = -3.0;
    B[0] = 4.0;
    C[0] = -1.0;
    D[0] = 2.0 * h;

    A[n] = 1.0;           
    B[n] = -4.0;          
    C[n] = 3.0 - 2.0 * h; 
    D[n] = 2.0 * h;       

    // Thomas algorithm for tridiagonal system
    vector<double> alpha(n + 1), beta(n + 1);

    // Forward sweep
    alpha[0] = C[0] / B[0];
    beta[0] = D[0] / B[0];

    for (int i = 1; i < n; i++) {
        double denom = B[i] - A[i] * alpha[i - 1];
        alpha[i] = C[i] / denom;
        beta[i] = (D[i] - A[i] * beta[i - 1]) / denom;
    }

    int m = n + 1;
    vector<double> a_vec(m), b_vec(m), c_vec(m), d_vec(m);

    for (int i = 1; i <= n - 1; i++) {
        a_vec[i] = A[i];
        b_vec[i] = B[i];
        c_vec[i] = C[i];
        d_vec[i] = D[i];
    }

    a_vec[0] = 0;
    b_vec[0] = B[0];
    c_vec[0] = C[0];
    d_vec[0] = D[0];


    b_vec[n] = -4.0 / (2.0 * h) - 1.0;  
    a_vec[n] = 1.0 / (2.0 * h);        
    c_vec[n] = 3.0 / (2.0 * h);         
    d_vec[n] = 1.0;                 

    
    vector<double> diag(n + 1, 0.0), sub(n + 1, 0.0), sup(n + 1, 0.0), rhs(n + 1, 0.0);

    // Interior points
    for (int i = 1; i < n; i++) {
        double xi = x[i];
        double exp_xi = exp(xi);
        double denom = exp_xi + 1.0;
        double p = -2.0 / denom;
        double q = -exp_xi / denom;

        sub[i] = 1.0 / (h * h) - p / (2.0 * h);    
        diag[i] = -2.0 / (h * h) + q;          
        sup[i] = 1.0 / (h * h) + p / (2.0 * h);    
        rhs[i] = 0.0;
    }


    double x0 = x[0];
    double exp_x0 = exp(x0);
    double denom0 = exp_x0 + 1.0;
    double p0 = -2.0 / denom0;
    double q0 = -exp_x0 / denom0;


    diag[0] = -2.0 / (h * h) + q0;
    sup[0] = 2.0 / (h * h);
    rhs[0] = 2.0 / h - p0;


    double xn = x[n];
    double exp_xn = exp(xn);
    double denomn = exp_xn + 1.0;
    double pn = -2.0 / denomn;
    double qn = -exp_xn / denomn;

    sub[n] = 2.0 / (h * h);
    diag[n] = -2.0 / (h * h) + 2.0 / h + pn + qn;
    rhs[n] = -2.0 / h - pn;

    // Thomas algorithm for tridiagonal system
    vector<double> c_prime(n + 1, 0.0), d_prime(n + 1, 0.0);

    c_prime[0] = sup[0] / diag[0];
    d_prime[0] = rhs[0] / diag[0];

    for (int i = 1; i < n; i++) {
        double denom = diag[i] - sub[i] * c_prime[i - 1];
        c_prime[i] = sup[i] / denom;
        d_prime[i] = (rhs[i] - sub[i] * d_prime[i - 1]) / denom;
    }

    // Last row
    d_prime[n] = (rhs[n] - sub[n] * d_prime[n - 1]) / (diag[n] - sub[n] * c_prime[n - 1]);

    // Back substitution
    y[n] = d_prime[n];
    for (int i = n - 1; i >= 0; i--) {
        y[i] = d_prime[i] - c_prime[i] * y[i + 1];
    }

    return y;
}


double runge_romberg_error(const vector<double>& y1, const vector<double>& y2, int p) {
    int n = y1.size() - 1;
    double max_error = 0.0;

    for (int i = 0; i <= n; i++) {
        double error = abs(y1[i] - y2[2 * i]) / (pow(2, p) - 1);
        if (error > max_error) {
            max_error = error;
        }
    }

    return max_error;
}

int main() {
    cout << fixed << setprecision(10);

    // Parameters
    double a = 0.0, b = 1.0;
    int n1 = 20;    
    int n2 = 40;    

    cout << "=== Finite Difference Method for BVP ===" << endl;
    cout << "Equation: (e^x + 1)y'' - 2y' - e^x*y = 0" << endl;
    cout << "Boundary conditions: y'(0) = 1, y'(1) - y(1) = 1" << endl;
    cout << "Exact solution: y(x) = e^x - 1" << endl << endl;

    // Grid points
    vector<double> x1(n1 + 1);
    double h1 = (b - a) / n1;
    for (int i = 0; i <= n1; i++) {
        x1[i] = a + i * h1;
    }

    // Solve using finite differences
    vector<double> y_fd1 = solve_bvp_fd(n1);
    vector<double> y_fd2 = solve_bvp_fd(n2);

    // Get shooting method solution for comparison
    vector<double> y_exact1(n1 + 1);
    for (int i = 0; i <= n1; i++) {
        y_exact1[i] = exact_solution(x1[i]);
    }

    // Calculate errors
    double max_abs_error = 0.0;
    double max_rel_error = 0.0;
    int max_error_idx = 0;

    cout << "=== Comparison with Exact Solution ===" << endl;
    cout << setw(10) << "x" << setw(20) << "y_FD" << setw(20) << "y_exact"
        << setw(20) << "Abs Error" << setw(20) << "Rel Error" << endl;
    cout << string(90, '-') << endl;

    for (int i = 0; i <= n1; i++) {
        double abs_error = abs(y_fd1[i] - y_exact1[i]);
        double rel_error = (abs(y_exact1[i]) > 1e-12) ? abs_error / abs(y_exact1[i]) : abs_error;

        if (abs_error > max_abs_error) {
            max_abs_error = abs_error;
            max_rel_error = rel_error;
            max_error_idx = i;
        }

        if (i % 2 == 0) {
            cout << setw(10) << x1[i]
                << setw(20) << y_fd1[i]
                << setw(20) << y_exact1[i]
                << setw(20) << abs_error
                << setw(20) << rel_error << endl;
        }
    }

    cout << endl << "=== Error Analysis ===" << endl;
    cout << "Grid size: " << n1 << " (h = " << h1 << ")" << endl;
    cout << "Maximum absolute error: " << max_abs_error << endl;
    cout << "Maximum relative error: " << max_rel_error << endl;
    cout << "Runge-Romberg error estimate (p=2): "
        << runge_romberg_error(y_fd1, y_fd2, 2) << endl;
    cout << "x where max error occurs: " << x1[max_error_idx] << endl;

    return 0;
}