#include <preprocess/preprocess.h>
#include <fstream>
#include <cstdlib>
#include <sstream>

// flag - true -> plus 
double PreProcess::get_alpha(int i, int j, int point, bool flag) {
    if (!flag) {
        return pow(((point + j - 1) * h - alpha) / h, i);
    }
    return pow(((point - (j - 1)) * h - alpha) / h, i);
}
Matrix2d PreProcess::GetDefaultQ (int l, int point, bool flag) {
    double rho_temp, k_temp, c_temp;
    const double eps = 1e-12;
    rho_temp = rho_minus / rho_plus;
    k_temp = k_minus / k_plus;
    c_temp = c_minus / c_plus;
    if (flag) {
        rho_temp = 1 / rho_temp;
        k_temp = 1 / k_temp;
        c_temp = 1 / c_temp;
    }
    Vector2d q2_diag (k_temp, 1 / rho_temp);

    Matrix2d q1 = Matrix2d::Identity();
    Matrix2d q2 = q2_diag.asDiagonal();
    Matrix2d q3 = c_temp * c_temp * Matrix2d::Identity();

    vector<Matrix2d> result {q1, q2, q3};
    return result[l];
}
Matrix2d PreProcess::GetQmatrix(int i, int l, int point, bool flag) {
    Matrix2d q1 = get_alpha(i, l, point, flag) * Matrix2d::Identity();
    Matrix2d q2 = get_alpha(i, l, point, flag) * Matrix2d::Identity();
    Matrix2d q3 = get_alpha(i, l, point, flag) * GetDefaultQ(i, point, flag);

    vector<Matrix2d> q_s = {q1, q2, q3};
    return q_s[l];
}
Matrix2d PreProcess::GetFmatrix(int i, int point, bool flag) {
    Matrix2d _A = A_minus;

    if (flag) {
        _A = A_plus;
    }
    Matrix2d f1 = Matrix2d::Zero();
    Matrix2d f2 = -_A;
    Matrix2d f3 = -2 * get_alpha(1, 1, point, flag) * _A + (tau / h) * _A * _A;
    vector<Matrix2d> matrices = {f1, f2, f3};
    return matrices[i];
}
vector<Matrix2d> PreProcess::CalcGammaMatrices(int point, bool flag) {
    Matrix<double, 6, 6> Q;
    Eigen::Matrix<double, 6, 2> F;
    Q.setZero();
    F.setZero();

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Q.block<2, 2>(2 * i, 2 * j) = GetQmatrix(i, j, point, flag).transpose();
        }
    }
    for (int i = 0; i < 3; ++i) {
        F.block<2, 2>(2 * i, 0) = GetFmatrix(i, point, flag).transpose();
    }
    vector<Matrix2d> block_matrices;
    Matrix<double, 6, 2> G = Q.colPivHouseholderQr().solve(F);

    for (int i = 0; i < 6; i += 2) {
        Matrix2d block = G.block<2, 2>(i, 0);
        block_matrices.push_back(block.transpose());
    }
    return block_matrices;
}
void PreProcess::Solve() {
    for (int i = 0; i < 2; ++i) {
        vector<Matrix2d> gammas = CalcGammaMatrices(i + J, i);
        gamma_matrices[i + J] = gammas;
    }
}
PreProcess::PreProcess(double _cir_left, int _M, double _x0, double _A, double _omega)
    : InitValues(_cir_left, _M, _x0, _A, _omega) {
    Solve();
}
