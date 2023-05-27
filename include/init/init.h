#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <Eigen/Dense>
using namespace Eigen;
using std::vector;

class InitValues {
private:
    // vector<double> coord_x, coord_y, velocity_x, velocity_y;
    double A, omega, x0, y0;
public:
    vector<Vector2d> u;
    Matrix2d A_minus, A_plus;
    double rho_minus, rho_plus, c_minus, c_plus, k_minus, k_plus;
    double tau, h, alpha;

    int J;
    int Mx;

    InitValues();
    InitValues(double _tau, double _h, int _M, double _x0, double _A, double _omega);
    
    void SetInitRadU(double x0, double omega);
    void PrintInit();
};
