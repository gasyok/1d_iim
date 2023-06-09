#pragma once
#include <init/init.h>

class PreProcess : public InitValues {
private:
    Matrix2d GetDefaultQ (int l, int point);
    Matrix2d GetQmatrix(int i, int l, int point);
    Matrix2d GetFmatrix(int i, int point);
    vector<Matrix2d> CalcGammaMatrices(int point);
    double get_alpha(int i, int j, int point);
    void Solve();
public:
    std::unordered_map<int, std::vector<Eigen::Matrix2d>> gamma_matrices;
    vector<Matrix2d> gamma_minus, gamma_plus;
    PreProcess(double _tau, int _M, double _x0, double _A, double _omega);
};
