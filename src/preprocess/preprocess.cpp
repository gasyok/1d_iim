#include <preprocess/preprocess.h>
#include <fstream>
#include <cstdlib>
#include <sstream>

// flag - true -> plus 
double PreProcess::get_alpha(int i, int j, int point) {
    return pow(((point + j - 2) * h - alpha) / h, i);
}
Matrix2d PreProcess::GetDefaultQ (int l, int point) {
    double rho_temp, k_temp, c_temp;
    const double eps = 1e-12;
    rho_temp = rho_minus / rho_plus;
    k_temp = k_minus / k_plus;
    c_temp = c_minus / c_plus;

    if (point > J) {
        rho_temp = 1 / rho_temp;
        k_temp = 1 / k_temp;
        c_temp = 1 / c_temp;
    }
    Vector2d q2_diag (k_temp, 1 / rho_temp);
    Vector2d q4_diag (k_temp * k_temp / rho_temp, k_temp / (rho_temp * rho_temp));

    Matrix2d q1 = Matrix2d::Identity();
    Matrix2d q2 = q2_diag.asDiagonal();
    Matrix2d q3 = c_temp * c_temp * Matrix2d::Identity();
    Matrix2d q4 = q4_diag.asDiagonal();
    Matrix2d q5 = c_temp * c_temp * c_temp * c_temp * Matrix2d::Identity();

    vector<Matrix2d> result {q1, q2, q3, q4, q5};
    return result[l];
}
Matrix2d PreProcess::GetQmatrix(int i, int l, int point) {
    // const vector<int> offsets {point - 2, point - 1, point, point + 1};
    if ((point * h - alpha) * ((point - 2 + l) * h - alpha) < 0) {
        return get_alpha(i, l, point) * GetDefaultQ(i, point);
    }
    return get_alpha(i, l, point) * Matrix2d::Identity();
}
Matrix2d PreProcess::GetFmatrix(int i, int point) {
    Matrix2d _A = A_minus;

    if (point > J) {
        _A = A_plus;
    }
    Matrix2d f1 = Matrix2d::Zero();
    Matrix2d f2 = -_A;
    Matrix2d f3 = -2 * get_alpha(1, 2, point) * _A + (tau / h) * _A * _A;
    Matrix2d f4 = -3 * get_alpha(2, 2, point) * _A + 3 * (tau / h) * get_alpha(1, 2, point) * _A * _A - (tau / h) * (tau / h) * _A * _A * _A;
    Matrix2d f5 = -4 * get_alpha(3, 2, point) * _A + 6 * (tau / h) * get_alpha(2, 2, point) * _A * _A - 4 * (tau / h) * (tau / h) * _A * _A * _A * get_alpha(1, 2, point) + (tau / h) * (tau / h) * (tau / h) * _A * _A * _A * _A;

    vector<Matrix2d> matrices = {f1, f2, f3, f4, f5};
    return matrices[i];
}
vector<Matrix2d> PreProcess::CalcGammaMatrices(int point) {
    Matrix<double, 10, 10> Q;
    Eigen::Matrix<double, 10, 2> F;
    Q.setZero();
    F.setZero();

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            Q.block<2, 2>(2 * i, 2 * j) = GetQmatrix(i, j, point).transpose();
        }
    }
    for (int i = 0; i < 5; ++i) {
        F.block<2, 2>(2 * i, 0) = GetFmatrix(i, point).transpose();
    }
    vector<Matrix2d> block_matrices;
    // Eigen::JacobiSVD<Eigen::MatrixXd> svd(Q, Eigen::ComputeThinU | Eigen::ComputeThinV);
    // Eigen::MatrixXd G = svd.solve(F);
    Matrix<double, 10, 2> G = Q.colPivHouseholderQr().solve(F);

    for (int i = 0; i < 10; i += 2) {
        Matrix2d block = G.block<2, 2>(i, 0);
        block_matrices.push_back(block.transpose());
    }

    return block_matrices;
}
void PreProcess::Solve() {
    for (int i = J-1; i <= J + 2; ++i) {
        vector<Matrix2d> gammas = CalcGammaMatrices(i);
        gamma_matrices[i] = gammas;
    }
}
PreProcess::PreProcess(double _cir_left, int _M, double _x0, double _A, double _omega)
    : InitValues(_cir_left, _M, _x0, _A, _omega) {
    Solve();

    // for (int i = J; i <= J + 2; ++i) {
    //     std::cout << "\n\nGamma Matrices for i = " << i << "\n\n";
    //     for (auto c : gamma_matrices[i]) {
    //         std::cout << c << "\n\n";
    //     }
    // }
    // std::cout << "\n\nChecker:\n";
}
