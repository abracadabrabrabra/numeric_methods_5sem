#define MAX_CNT 10000
#define DEBUG
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
using namespace std;

double matrix_norm(const vector<vector<double>>& m) {
    double max = -1.;
    for (size_t i = 0; i < m.size(); ++i) {
        double cur = 0.;
        for (size_t j = 0; j < m.size(); ++j)
            cur += fabs(m[i][j]);
        if (cur > max)
            max = cur;
    }
    return max;
}

double norm(const vector<double>& x_new, const vector<double>& x_old) {
    double res = 0.;
    for (size_t i = 0; i < x_new.size(); ++i)
        res += fabs(x_new[i] - x_old[i]);
    return res;
}

void print_matrix(const vector<vector<double>>& m) {
    for (int i = 0; i < m.size(); ++i) {
        for (int j = 0; j < m.size(); ++j)
            cout << m[i][j] << " ";
        cout << endl;
    }
}

bool diagonal_dominance(vector<vector<double>>& A) {
    size_t sz = A.size();
    for (size_t i = 0; i < sz; ++i) {
        double row_sum = 0.0;
        for (size_t j = 0; j < sz; ++j) {
            if (i != j)
                row_sum += fabs(A[i][j]);
        }
        if (fabs(A[i][i]) <= row_sum)
            return false;
    }
    return true;
}

void iter_method(const vector<vector<double>>& A, const vector<double>& right_parts,
    vector<double>& x_new, const double epsilon) {

    size_t sz = right_parts.size();
    vector<vector<double>> alpha(sz, vector<double>(sz, 0.));
    vector<double> beta(sz, 0.);

    for (size_t i = 0; i < sz; ++i) {
        alpha[i][i] = 0.;
        beta[i] = right_parts[i] / A[i][i];
        for (size_t j = 0; j < sz; ++j) {
            if (i != j)
                alpha[i][j] = -A[i][j] / A[i][i];
        }
    }

    size_t iter_cnt = 0;
    double precision;
    vector<double> x_old(beta);
    cout << endl;
    for (size_t i = 0; i < sz; ++i)
        cout << "b" << i + 1 << " = " << beta[i] << endl;
    cout << "\nalpha:\n";
    print_matrix(alpha);
    double norm_alpha = matrix_norm(alpha), koeff = norm_alpha / (1 - norm_alpha);

    while (iter_cnt < MAX_CNT) {
        x_new.assign(sz, 0.);
        ++iter_cnt;

        #pragma omp parallel for
        for (size_t i = 0; i < sz; ++i) {
            x_new[i] = beta[i];
            for (size_t j = 0; j < i; ++j)
                x_new[i] += alpha[i][j] * x_old[j];
            for (size_t j = i + 1; j < sz; ++j)
                x_new[i] += alpha[i][j] * x_old[j];

        }
        precision = norm(x_new, x_old) * koeff;

        #ifdef DEBUG
        cout << "\n iter " << iter_cnt << endl;
        for (size_t i = 0; i < sz; ++i)
            cout << "x" << i + 1 << " = " << x_new[i] << endl;
        cout << "prec " << precision << endl;
        #endif

        if (precision < epsilon)
            break;
        x_old = x_new;
    }
    cout << endl << "Solution after " << iter_cnt << " iterations" << endl;
}

int main() {

    double eps;
    cout << "Input precision: ";
    cin >> eps;

    size_t sz;
    cout << "Input size of system: ";
    cin >> sz;
    if (!sz)
        throw invalid_argument("System size cannot be zero");

    vector<vector<double>> A;
    A.reserve(sz);
    vector<double> right_parts, res;
    double x;

    cout << "Input matrix of system and right parts" << endl;
    for (size_t i = 0; i < sz; ++i) {
        for (size_t j = 0; j < sz; ++j) {
            cin >> x;
            if (i == j && fabs(x) < 1e-10)
                throw runtime_error("Diagonal element is zero");
            A[i].push_back(x);
        }
        cin >> x;
        right_parts.push_back(x);
    }

    if (!diagonal_dominance(A))
        throw invalid_argument("No diagonal dominance");

    iter_method(A, right_parts, res, eps);

    for (size_t i = 0; i < sz; ++i)
        cout << "x" << i + 1 << " = " << res[i] << endl;

    return 0;
}