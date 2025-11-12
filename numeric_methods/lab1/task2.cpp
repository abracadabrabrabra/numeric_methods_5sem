#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

/*
14 9 125
-8 14 6 -56
-5 -17 8 144
1 5 -2 36
-4 -10 70
*/

int main() {

    int sz;
    vector<double> P, Q, X;
    double a, b, c, d;

    cout << "Input size of system: " << endl;
    cin >> sz;
    if (sz == 0)
        throw invalid_argument("System size cannot be zero");
    cout << "Input non-zero diagonal elements and right parts: " << endl;
    X.resize(sz);

    for (int i = 0; i < sz; ++i) {

        if (i == 0) {
            a = 0.;
            cin >> b >> c >> d;
            if (abs(b) < abs(c))
                throw invalid_argument("No diagonal dominance");
            P.push_back(-c / b);
            Q.push_back(d / b);

        }
        else if (i == sz - 1) {
            c = 0.;
            cin >> a >> b >> d;
            if (abs(b) < abs(a))
                throw invalid_argument("No diagonal dominance");
            P.push_back(0.);
            Q.push_back((d - a * Q[i - 1]) / (b + a * P[i - 1]));
        }
        else {
            cin >> a >> b >> c >> d;
            if (abs(b) < abs(a) + abs(c))
                throw invalid_argument("No diagonal dominance");
            P.push_back(-c / (b + a * P[i - 1]));
            Q.push_back((d - a * Q[i - 1]) / (b + a * P[i - 1]));
        }

    }

    X[sz - 1] = Q[sz - 1];
    cout << "\nSolution:\nx" << sz << "=" << X[sz - 1] << endl;
    for (int i = sz - 2; i >= 0; --i) {
        X[i] = P[i] * X[i + 1] + Q[i];
        cout << "x" << i + 1 << "=" << X[i] << endl;
    }

    return 0;
}