#pragma once
#include <preprocess/preprocess.h>
#include <unordered_map>

using namespace Eigen;

class System : public PreProcess {
public:
    System(double _tau, double _h, int _M, double _x0, double _A, double _omega);
    Vector2d GetValue(int i);
    Vector2d equation(int i, bool flag);
    Vector2d irrEquation(int i, bool flag);
    void solve(double t);
};
