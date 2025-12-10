#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <iomanip>
using namespace std;

using func2 = double(*)(double, double, double);

//y'' = -2y'ctg(x) - 3y
double f2(double x, double y, double z) {
    return -2 * z / tan(x) - 3 * y;
}

double exact_solve(double x) {
    return (0.4776 * sin(2 * x) - 0.9783 * cos(2 * x)) / sin(x);
}

struct system_input_data {
    const double x_0, y_0, z_0, x_n;
    func2 f;
};

vector<tuple<double, double, double>> euler_system(system_input_data& input, vector<double>& y, const double h) {
    const double eps = 1e-12;
    if (h <= eps)
        throw invalid_argument("Step too small");
    if (input.x_0 > input.x_n || fabs(input.x_0 - input.x_n) < eps)
        throw invalid_argument("Incorrect range: x_0 should be less than x_n");

    int n = static_cast<int>((input.x_n - input.x_0) / h) + 1;
    if (n < 2)
        throw invalid_argument("Not enough elements in net function");

    vector<tuple<double, double, double>> net(n);
    y.resize(n);

    net[0] = { input.x_0, input.y_0, input.z_0 };
    y[0] = input.y_0;

    for (int i = 1; i < n; ++i) {
        double x_prev = get<0>(net[i - 1]);
        double y_prev = get<1>(net[i - 1]);
        double z_prev = get<2>(net[i - 1]);

        // first half step
        double x_mid = x_prev + h / 2.0;
        double y_mid = y_prev + (h / 2.0) * z_prev;
        double z_mid = z_prev + (h / 2.0) * input.f(x_prev, y_prev, z_prev);

        // second half step
        double x_next = x_prev + h;
        double new_y = y_prev + h * z_mid;
        double new_z = z_prev + h * input.f(x_mid, y_mid, z_mid);

        y[i] = new_y;
        net[i] = { x_next, new_y, new_z };
    }

    return net;
}

using method = vector<tuple<double, double, double>>(*)(system_input_data&, vector<double>&, const double);

void estimate_error(method m, system_input_data& input, const double h, int order) {
    vector<double> y_h, y_h2;

    auto net_h = m(input, y_h, h);
    auto net_h2 = m(input, y_h2, h / 2);

    cout << "x\t\ty_h\t\ty_h/2\t\ty_exact\t\tabs_error\tRR_error" << endl;

    for (size_t i = 0; i < net_h.size(); ++i) {
        double x = get<0>(net_h[i]);
        double y_h_val = get<1>(net_h[i]);
        double y_exact = exact_solve(x);
        double abs_error = fabs(y_h_val - y_exact);

        size_t j = i * 2;
        if (j < net_h2.size()) {
            double y_h2_val = get<1>(net_h2[j]);

            double rr_error = fabs(y_h_val - y_h2_val) / (pow(2, order) - 1);

            cout << fixed << setprecision(7) << x << "\t" << y_h_val << "\t" << y_h2_val
                << "\t" << y_exact << "\t" << abs_error << "\t" << rr_error << endl;
        }
    }
}

int main() {
    system_input_data input{ 1., 1., 1., 2., f2 };

    cout << "=== First improved Euler method ===" << endl;
    estimate_error(euler_system, input, 0.1, 2);

    return 0;
}