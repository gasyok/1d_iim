#include <init/init.h>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <sstream>

void InitValues::PrintInit() {
    std::ostringstream fileinit;
    fileinit << "../bin/animation/init_out.bin";
    std::ofstream file(fileinit.str(), std::ios::binary);
    for (int i = 0; i < Mx; ++i) {
        float x_coord = i * h;
        float pressure_value = u[i](1);
        // Запись x, y и p(x, y) (каждый параметр - 4 байта, little-endian)
        file.write(reinterpret_cast<char*>(&x_coord), sizeof(x_coord));
        file.write(reinterpret_cast<char*>(&pressure_value), sizeof(pressure_value));
    }
    file.close();
}
double InitValues::foo(double x, bool flag = false) {
    double _omega = omega;
    if (flag) {
        _omega = omega * c_plus / c_minus;
    }
    if (x < _omega && x > 0.0) return pow(sin(2 * M_PI * x / (2 * _omega)), 4);
    return 0.0;
}
Vector2d InitValues::GetExactSol(int i, double x0, double omega, double t) {
    double x = i * h;
    double _pressure = 0.0;
    double omega_t = omega * c_plus / c_minus;

    if (i <= J) {
        _pressure = A * foo(x - x0 - c_minus * t);
        _pressure += A_r * foo(x + x0 + omega - 2 * alpha + c_minus * t);
        return Vector2d(1 / (rho_minus * c_minus) * _pressure, _pressure);
    }
    else {
        _pressure = A_t * foo(x - alpha + omega_t + c_plus / c_minus * (alpha - x0 - omega) - c_plus * t, true);
        return Vector2d(1 / (rho_plus * c_plus) * _pressure, _pressure);
    }
}
void InitValues::SetInitRadU(double x0, double omega) {
    double _pressure, vel;
    for (int i = 0; i < Mx; ++i) {
        if (i <= J) {
            _pressure = foo(i * h - x0);
            vel = 1 / (rho_minus * c_minus) * _pressure;
        }
        else {
            _pressure = 0;
            vel = 0.0;
        }
        u.push_back(Vector2d(vel, _pressure));
    }
}
InitValues::InitValues(double cir_left, int _M, double _x0, double _A, double _omega)
: h(2.0 / _M), Mx(_M), x0(_x0), omega(_omega), A(_A) {

    rho_minus = 1;
    c_minus = 2;

    rho_plus = 0.8;
    c_plus = 1.5;
    k_minus = c_minus * c_minus * rho_minus;
    k_plus = c_plus * c_plus * rho_plus;

    z_l = (rho_minus * c_minus);
    z_r = (rho_plus * c_plus);
    A_t = A * 2 * z_r / (z_l + z_r);
    A_r = A * (z_r - z_l) / (z_l + z_r);
    tau = h * cir_left / c_minus;

    std::cout << "A_t = " << A_t << std::endl;
    std::cout << "A_r = " << A_r << std::endl;
    std::cout << "TAU: " << tau << std::endl;
    std::cout << "CIR LEFT: " << c_minus * tau / h << std::endl;
    std::cout << "CIR RIGHT: " << c_plus * tau / h << std::endl;

    A_minus << 0, 1 / rho_minus,
        k_minus, 0;
    A_plus << 0, 1 / rho_plus,
        k_plus, 0;

    alpha = 1.0 + 1e-10;
    J = alpha / h;
    if (h * (J + 1) <= alpha) {
        std::cout << "Error in J\n";
        exit(-1);
    }

    SetInitRadU(x0, omega);
}
