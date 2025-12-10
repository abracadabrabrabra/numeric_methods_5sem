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

vector<tuple<double, double, double>> rk4_system(system_input_data& input, vector<double>& y, const double h) {
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

        double k1, k2, k3, k4;
        double l1, l2, l3, l4;

        k1 = h * z_prev;
        l1 = h * input.f(x_prev, y_prev, z_prev);

        k2 = h * (z_prev + l1 / 2);
        l2 = h * input.f(x_prev + h / 2, y_prev + k1 / 2, z_prev + l1 / 2);

        k3 = h * (z_prev + l2 / 2);
        l3 = h * input.f(x_prev + h / 2, y_prev + k2 / 2, z_prev + l2 / 2);

        k4 = h * (z_prev + l3);
        l4 = h * input.f(x_prev + h, y_prev + k3, z_prev + l3);

        double delta_y = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        double delta_z = (l1 + 2 * l2 + 2 * l3 + l4) / 6;

        double new_x = x_prev + h;
        double new_y = y_prev + delta_y;
        double new_z = z_prev + delta_z;

        y[i] = new_y;
        net[i] = { new_x, new_y, new_z };
    }

    return net;
}

vector<tuple<double, double, double>> adams_system(system_input_data& input, vector<double>& y, const double h) {
    const double eps = 1e-12;
    if (h <= eps)
        throw invalid_argument("Step too small");
    if (input.x_0 > input.x_n || fabs(input.x_0 - input.x_n) < eps)
        throw invalid_argument("Incorrect range: x_0 should be less than x_n");

    int n = static_cast<int>((input.x_n - input.x_0) / h) + 1;
    if (n < 5)
        throw invalid_argument("Not enough elements for 4-step Adams method");

    vector<tuple<double, double, double>> net(n);
    y.resize(n);

    net[0] = { input.x_0, input.y_0, input.z_0 };
    y[0] = input.y_0;

    //start filling: rk4
    for (int i = 1; i < 4; ++i) {
        double x_prev = get<0>(net[i - 1]);
        double y_prev = get<1>(net[i - 1]);
        double z_prev = get<2>(net[i - 1]);

        double k1 = h * z_prev;

        double l1 = h * input.f(x_prev, y_prev, z_prev);

        double k2 = h * (z_prev + l1 / 2);
        double l2 = h * input.f(x_prev + h / 2, y_prev + k1 / 2, z_prev + l1 / 2);

        double k3 = h * (z_prev + l2 / 2);
        double l3 = h * input.f(x_prev + h / 2, y_prev + k2 / 2, z_prev + l2 / 2);

        double k4 = h * (z_prev + l3);
        double l4 = h * input.f(x_prev + h, y_prev + k3, z_prev + l3);

        double delta_y = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        double delta_z = (l1 + 2 * l2 + 2 * l3 + l4) / 6;

        double new_x = x_prev + h;
        double new_y = y_prev + delta_y;
        double new_z = z_prev + delta_z;

        y[i] = new_y;
        net[i] = { new_x, new_y, new_z };
    }

    //main loop: adams
    for (int i = 4; i < n; ++i) {
        auto [x_n, y_n, z_n] = net[i - 1];
        auto [x_n1, y_n1, z_n1] = net[i - 2];
        auto [x_n2, y_n2, z_n2] = net[i - 3];
        auto [x_n3, y_n3, z_n3] = net[i - 4];

        double f_n = input.f(x_n, y_n, z_n);
        double f_n1 = input.f(x_n1, y_n1, z_n1);
        double f_n2 = input.f(x_n2, y_n2, z_n2);
        double f_n3 = input.f(x_n3, y_n3, z_n3);

        double y_next = y_n + h / 24. * (55 * z_n - 59 * z_n1 + 37 * z_n2 - 9 * z_n3);
        double z_next = z_n + h / 24. * (55 * f_n - 59 * f_n1 + 37 * f_n2 - 9 * f_n3);
        double x_next = x_n + h;

        y[i] = y_next;
        net[i] = { x_next, y_next, z_next };
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

    cout << "=== Adams 4th order ===" << endl;
    estimate_error(adams_system, input, 0.1, 4);

    return 0;
}