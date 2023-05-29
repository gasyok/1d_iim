#pragma once
#include <init/init.h>

class PreProcess : public InitValues {
private:
    Matrix2d GetDefaultQ (int l, int point, bool flag);
    Matrix2d GetQmatrix(int i, int l, int point, bool flag);
    Matrix2d GetFmatrix(int i, int point, bool flag);
    vector<Matrix2d> CalcGammaMatrices(int point, bool flag);
    double get_alpha(int i, int j, int point, bool flag);
    void Solve();
public:
    std::unordered_map<int, std::vector<Eigen::Matrix2d>> gamma_matrices;
    PreProcess(double _tau, int _M, double _x0, double _A, double _omega);
};
