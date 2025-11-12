#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

/*
-6 -5 -3 -8 101
5 -1 -5 -4 51
-6 0 5 5 -53
-7 -2 8 5 -63
*/

void print_matrix(const vector<vector<double>>& m, const int& sz) {
    for (int i = 0; i < sz; ++i) {
        for (int j = 0; j < sz; ++j)
            cout << m[i][j] << " ";
        cout << endl;
    }
}

void gauss(const vector<vector<double>>& a, const vector<vector<double>>& l,
    const vector<double>& right, vector<double>& res, int sz) {
    vector<double> z(sz, 0.);

    for (int i = 0; i < sz; ++i) {
        z[i] = right[i];
        for (int k = 0; k < i; ++k)
            z[i] -= l[i][k] * z[k];
    }

    for (int i = sz - 1; i >= 0; --i) {
        res[i] = z[i];
        for (int k = sz - 1; k > i; --k)
            res[i] -= a[i][k] * res[k];
        res[i] /= a[i][i];
    }
}

int main() {


    double det = 1.;
    int sz, cnt_perms = 0;
    cout << "Input size of system: ";
    cin >> sz;
    if (!sz)
        throw invalid_argument("System size cannot be zero");
    vector<vector<double>> a(sz, vector<double>(sz, 0.)), u(sz, vector<double>(sz, 0.)), l(sz, vector<double>(sz, 0.));
    vector<double> right_parts(sz, 0.), res(sz, 0.);
    vector<int> perms;

    cout << endl << "Input vector of right parts:" << endl;
    for (int k = 0; k < sz; ++k)
        cin >> right_parts[k];

    cout << endl << "Input matrix of system:" << endl;
    for (int k = 0; k < sz; ++k) {
        perms.push_back(k);
        for (int i = 0; i < sz; ++i) {
            cin >> u[k][i];
            a[k][i] = u[k][i];
        }
    }

    for (int k = 0; k < sz; ++k) {
        int lead_coeff = abs(u[k][k]), lead_string = k;
        l[k][k] = 1.;
        for (int i = k + 1; i < sz; ++i) {
            if (abs(u[i][k]) > lead_coeff) {
                lead_coeff = abs(u[i][k]);
                lead_string = i;
            }
        }
        if (lead_coeff < 1e-10)
            throw runtime_error("Zero div error");
        if (k != lead_string) {
            ++cnt_perms;
            u[k].swap(u[lead_string]);
            swap(perms[k], perms[lead_string]);
            swap(right_parts[k], right_parts[lead_string]);
            for (int j = 0; j < k; ++j)
                swap(l[k][j], l[lead_string][j]);
        }
        for (int i = k + 1; i < sz; ++i) {
            double m = u[i][k] / u[k][k];
            l[i][k] = m;
            for (int j = k + 1; j < sz; ++j)
                u[i][j] -= m * u[k][j];
            u[i][k] = 0.;
        }
        det *= u[k][k];
    }
    if (fabs(det) < 1e-10)
        throw invalid_argument("Singular matrix");
    if (cnt_perms % 2)
        det = -det;

    cout << endl << "Upper triangle matrix:" << endl;
    print_matrix(u, sz);
    cout << endl << "Lower triangle matrix:" << endl;
    print_matrix(l, sz);
    cout << endl << "Determinant = " << det << endl;

    cout << endl << "System solution:" << endl;
    gauss(u, l, right_parts, res, sz);
    for (int i = 0; i < sz; ++i)
        cout << "x" << i + 1 << " = " << res[i] << endl;

    vector<vector<double>> inverse(sz, vector<double>(sz, 0.));
    for (int i = 0; i < sz; ++i) {
        vector<double> unit_col(sz, 0.);
        for (int j = 0; j < sz; ++j) {
            if (perms[j] == i)
                unit_col[j] = 1.;
        }
        gauss(u, l, unit_col, inverse[i], sz);
    }

    cout << endl << "Inverse matrix:" << endl;
    for (int i = 0; i < sz; ++i) {
        for (int k = 0; k < sz; ++k)
            cout << inverse[k][i] << " ";
        cout << endl;
    }

    cout << endl << "System matrix" << endl;
    print_matrix(a, sz);

    cout << endl << "Inverse multiply check" << endl;
    for (int i = 0; i < sz; ++i) {
        for (int j = 0; j < sz; ++j) {
            double x = 0.;
            for (int k = 0; k < sz; ++k)
                x += a[i][k] * inverse[j][k];
            cout << x << " ";
        }
        cout << endl;
    }

    return 0;
}