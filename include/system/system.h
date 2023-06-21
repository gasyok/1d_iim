#pragma once
#include <preprocess/preprocess.h>
#include <unordered_map>

using namespace Eigen;

class System : public PreProcess {
private:
    double du1(int m);
    double du2(int m);
public:
    double l1, lmax;
    double total_time;
    int total_steps;
    System(double _tau, int _M, double _x0, double _A, double _omega);
    Vector2d GetValue(int i);
    Vector2d step_one(int i, bool is_forward);
    Vector2d step_two(int i);
    Vector2d equation(int i);
    double equation_1(int i);
    double equation_2(int i);
    Vector2d irrEquation(int i);
    void solve(double t);
    void sample();
    double l_1(double t);
    double l_inf(double t);
    double tmp_1(int i);
    double tmp_2(int i);
};
