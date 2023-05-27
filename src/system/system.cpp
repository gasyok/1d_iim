#include <system/system.h>
#include <fstream>

Vector2d System::GetValue(int i) {
    return u[i];
}
Vector2d System::equation(int i, bool flag) {
    Matrix2d _A = A_minus;

    if (flag) {
        _A = A_plus;
    }
    return u[i] - 0.5 * (tau / h) * _A * (u[i + 1] - u[i - 1]) + 0.5 * (tau / h) * (tau / h) * _A * _A * (u[i + 1] - 2 * u[i] + u[i - 1]);
}
Vector2d System::irrEquation(int i, bool flag) {
    Vector2d res (0, 0);
    const vector<int> offsets = {i - 1, i, i + 1};
    for (int l = 0; l < 3; ++l) {
        int k = l;
        if (flag) {
            k = 2 - l;
        }
        int new_i = offsets[k];
        res += gamma_matrices[i][k] * u[new_i];
    }
    return u[i] + (tau / h) * res; 
}
void System::solve(double t) {
    vector<Vector2d> new_u = u;
    for (int i = 1; i < Mx - 1; ++i) {
        if (i == J) {
            new_u[i] = irrEquation(i, false);
            new_u[i + 1] = irrEquation(i, true);
            i = i + 1;
        }
        else {
            if (i > J + 1) {
                new_u[i] = equation(i, true);
            }
            else {
                new_u[i] = equation(i, false);
            }
        }
    }
    u = new_u;
}
System::System(double _tau, double _h, int _M, double _x0, double _A, double _omega)
    : PreProcess(_tau, _h, _M, _x0, _A, _omega) {}
