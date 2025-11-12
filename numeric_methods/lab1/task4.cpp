#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
using namespace std;
using matrix = vector<vector<double>>;

int sign(double x) { return copysign(1., x); }

void print_matrix(const matrix& m){
    for (int i = 0; i < m.size(); ++i){
        for (int j = 0; j < m.size(); ++j)
            cout << m[i][j] << " ";
        cout << endl;
    }
}

void unit_diag(matrix& m){
    for (size_t i = 0; i < m.size(); ++i)
        m[i][i] = 1.;
}

void matrix_multiply(const matrix& a, const matrix& b, matrix& res){
    for (size_t i = 0; i < a.size(); ++i){
        for (size_t j = 0; j < a.size(); ++j){
            res[i][j] = 0.;
            for (size_t k = 0; k < a.size(); ++k)
                res[i][j] += a[i][k] * b[k][j];
        }
    }
}

void qr_householder(const matrix& A, matrix& Q, matrix& R){
    vector<double> v(A.size());
    matrix A_prev = A, H_prev(A.size(), vector<double>(A.size())), H = H_prev;
    unit_diag(H_prev);
    
    for (size_t i = 0; i < A.size() - 1; ++i){
        for (size_t j = 0; j < i; ++j)
            v[j] = 0.;
        double sum = A_prev[i][i] * A_prev[i][i];
        for (size_t j = i + 1; j < A.size(); ++j){
            v[j] = A_prev[j][i];
            sum += v[j] * v[j];
        }
        v[i] = A_prev[i][i] + sign(A_prev[i][i])*pow(sum, 0.5);
        sum -= A_prev[i][i] * A_prev[i][i] - v[i] * v[i];
        
        for (size_t j = 0; j < A.size(); ++j){
            for (size_t k = 0; k < A.size(); ++k){
                H[j][k] = -2 / sum * v[j] * v[k];
                if (j == k)
                    H[j][k] += 1.;
            }
        }
        matrix_multiply(H, A_prev, R);
        A_prev = R;
        matrix_multiply(H_prev, H, Q);
        H_prev = H;
    }
}

pair<double, double> block_2x2(double a, double b, double c, double d) {
    double trace = a + d, det = a * d - b * c, discriminant = trace * trace - 4 * det;
    double real_part = trace / 2;
    double imag_part = sqrt(-discriminant) / 2;
    return { real_part, imag_part };
}

bool is_convergence(const matrix& A, const matrix& A_prev, double eps, vector<complex<double>>& eigens) {
    size_t sz = A.size(), i;
    for (i = 0; i < sz - 1;) {
        double sum = 0.;
        for (size_t j = i + 1; j < sz; ++j)
            sum += A[j][i] * A[j][i];
        sum = pow(sum, 0.5);
        if (sum > eps){
            auto p1 = block_2x2(A[i][i], A[i][i+1], A[i+1][i], A[i+1][i+1]);
            auto p2 = block_2x2(A_prev[i][i], A_prev[i][i+1], A_prev[i+1][i], A_prev[i+1][i+1]);
            if (fabs(p1.first - p2.first) > eps || fabs(p1.second - p2.second) > eps)
                return false;
            eigens[i] = complex<double>(p1.first, p1.second);
            eigens[i+1] = complex<double>(p1.first, -(p1.second));
            i += 2;
        }
        else{
            eigens[i] = complex<double>(A[i][i], 0);
            ++i;
        }
    }
    if (i == sz - 1) {
        if (fabs(A[i][i] - A_prev[i][i]) > eps)
            return false;
        eigens[i] = complex<double>(A[i][i], 0);
    }
    
    return true;
}

vector<complex<double>> qr_algo(matrix& A, double eps) {
    size_t sz = A.size(), MAX_CNT = 1000;
    matrix Q(sz, vector<double>(sz)), R(sz, vector<double>(sz)), A_prev(sz, vector<double>(sz));
    vector<complex<double>> eigens(sz);
    
    for (size_t iter = 0; iter < MAX_CNT; ++iter) {
        A_prev = A;
        qr_householder(A, Q, R);
        matrix_multiply(R, Q, A);
        
        if (is_convergence(A, A_prev, eps, eigens)) {
            cout << "Algo was stopped on iteration " << iter + 1 << endl;
            cout << "Final matrix: " << endl;
            print_matrix(A);
            return eigens;
        }
        
    }
    
    throw runtime_error("Iterations limit exceed. Convergence not reached"); 
}


int main(){
    
    //{{-1., 2., 9.}, {9., 3., 4.}, {8., -4., -6.}}
    
    size_t sz;
    cin >> sz;
    double eps;
    cin >> eps;
    
    matrix A(sz, vector<double>(sz));
    for (size_t i = 0; i < sz; ++i){
        for (size_t j = 0; j < sz; ++j)
            cin >> A[i][j];
    }
    
    matrix Q(sz, vector<double>(sz)), R(sz, vector<double>(sz));
    qr_householder(A, Q, R);
    
    cout << endl << "Q:" << endl;
    print_matrix(Q);
    cout << endl << "R:" << endl;
    print_matrix(R);
    
    matrix B(sz, vector<double>(sz));
    matrix_multiply(Q, R, B);
    cout << endl << "Multiplication check of qr_householder:" << endl;
    print_matrix(B);
    
    cout << endl << "Launch of qr_algo:" << endl;
    matrix A0 = A;
    auto eigens = qr_algo(A0, eps);
    
    cout << endl << "Founded eigen values" << endl;
    for (size_t i = 0; i < eigens.size(); ++i) {
        cout << "lambda" << i << " = " << eigens[i].real();
        if (fabs(eigens[i].imag()) > 1e-10) 
            cout << " + (" << eigens[i].imag() << ") * i";
        cout << endl;
    }

    return 0;
}
