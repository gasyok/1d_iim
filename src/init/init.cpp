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
        _omega = omega * (A * A - A_r * A_r) / (A_t * A_t) * (rho_plus * c_plus * c_plus) / (rho_minus * c_minus * c_minus);
    }
    if (x <= _omega && x >= 0.0) return 0.5 * (1 - cos(2 * M_PI * x / _omega));
    return 0.0;
}
Vector2d InitValues::GetExactSol(int i, double x0, double omega, double t) {
    double x = i * h;
    double _pressure = 0.0;
    double omega_t = omega * (A * A - A_r * A_r) / (A_t * A_t) * (rho_plus * c_plus * c_plus) / (rho_minus * c_minus * c_minus);

    // if (x < alpha) {
    //     _pressure = foo(x - c_minus * t);
    //     _pressure += A * A_r * foo(x - 2 * alpha + omega + c_minus * t);
    //     // _pressure += A_r * foo(x - 2 * alpha + c_minus * t);
    //     return Vector2d(1 / (rho_minus * c_minus) * _pressure, _pressure);
    // } else {
    //     // _pressure = A_t * foo(x - (alpha - omega) - c_plus * (t - (alpha - omega) / c_minus), true);
    //     _pressure = A_t * foo(x - (alpha - omega_t) - c_plus * (t - (alpha - omega) / c_minus), true);
    //     return Vector2d(1 / (rho_plus * c_plus) * _pressure, _pressure);
    // }
    if (t < (alpha - omega) / c_minus) {
        _pressure = foo(x - c_minus * t);
        return Vector2d(1 / (rho_minus * c_minus) * _pressure, _pressure);
    }
    else {
        if (i <= J) {
            _pressure = foo(x - c_minus * t);
            _pressure += A * A_r * foo(x - 2 * alpha + omega + c_minus * t);
            return Vector2d(1 / (rho_minus * c_minus) * _pressure, _pressure);
        }
        else {
            _pressure = A_t * foo(x - (alpha - omega_t) - c_plus * (t - (alpha - omega) / c_minus), true);
            return Vector2d(1 / (rho_plus * c_plus) * _pressure, _pressure);
        }
    }
}
// Vector2d InitValues::GetExactSol(int i, double x0, double omega, double t) {
//     double x = i * h;
//     double _pressure = 0.0;
//     if (x < alpha - c_minus * t) {
//         _pressure = A * foo(x - c_minus * t);
//         return Vector2d(1 / (rho_minus * c_minus) * _pressure, _pressure);
//     } else if (x < alpha + c_plus * t) {
//         _pressure = A * A_r * foo(x + c_minus * t) + A * A_t * foo(x - c_plus * t);
//         return Vector2d(1 / (rho_minus * c_minus) * _pressure, _pressure);
//     } else {
//         _pressure = A * A_t * foo(x - c_plus * t);
//         return Vector2d(1 / (rho_plus * c_plus) * _pressure, _pressure);
//     }
// }

// Vector2d InitValues::GetExactSol(int i, double x0, double omega, double t) {
//     double _pressure, _c, _rho;
//     double _n, _A, xi;
//     // double _A, A_r, A_t;
//     _pressure = 0.0;
//     if (t < alpha / c_minus) {
//         _pressure =  foo(i * h - c_minus * t);
//         return Vector2d(1 / (_rho * _c) * _pressure, _pressure);
//     }
//     else {
//         if (i <= J) {
//             _pressure = A_r * foo(i * h - 2 * alpha + omega + c_minus * t);
//             return Vector2d(1 / (_rho * _c) * _pressure, _pressure);
//         }
//         else {
//             _pressure =  A_t * foo(i * h - alpha - c_plus * (t - alpha / c_minus));
//             return Vector2d(1 / (_rho * _c) * _pressure, _pressure);
//         }
//     }
    // if (i <= J) {
    //     _c = c_minus;
    //     _rho = rho_minus;
    // }
    // else {
    //     _c = c_plus;
    //     _rho = rho_plus;
    // }
    // if (i * h <= alpha - c_minus * t) {
    //     _pressure =  foo(i * h - c_minus * t);
    //     return Vector2d(1 / (_rho * _c) * _pressure, _pressure);
    // }
    // if ((i * h <= alpha + c_plus * t) && (i * h >= alpha - c_minus * t)) {
    //     _pressure = A_r * foo(i * h + c_minus * t) + A_t * foo(i * h - c_plus * t);
    //     return Vector2d(1 / (_rho * _c) * _pressure, _pressure);
    // }
    // if (i * h > alpha + c_plus * t) {
    //     _pressure = A_t * foo(i * h - c_plus * t);
    //     return Vector2d(1 / (_rho * _c) * _pressure, _pressure);
    // }
    // std::cout << "ERROR\n";
    // return Vector2d(1 / (_rho * _c) * _pressure, _pressure);
    // if (i <= J) {
    //     _c = c_minus;
    //     _rho = rho_minus;
    // }
    // else {
    //     _c = c_plus;
    //     _rho = rho_plus;
    // }
    // xi = (i * h - x0) - _c * t;
    // if (xi <= 1 / omega && xi >= 0) {
    //     _pressure = 0.5 * _A * (1 - cos(2 * M_PI * omega * xi));
    // }
    // else _pressure = 0;
    //
    // return Vector2d(_n / (_c * _rho) * _pressure, _pressure);
// }
void InitValues::SetInitRadU(double x0, double omega) {
    double _pressure, vel;
    for (int i = 0; i < Mx; ++i) {
        if (i <= J) {
            _pressure = foo(i * h);
            vel = 1 / (rho_minus * c_minus) * _pressure;
        }
        else {
            _pressure = 0;
            vel = 0.0;
        }
        u.push_back(Vector2d(vel, _pressure));
        // u.push_back(GetExactSol(i, x0, omega, 0));
        // double xi = i * h - x0;
        // if (xi <= 1 / omega && xi >= 0) {
        //     _pressure = 0.5 * (1 - cos(2 * M_PI * omega * xi));
        // }
        // else {
        //     _pressure = 0;
        // }
        // u.push_back(Vector2d(1 / (c_minus * rho_minus) * _pressure, _pressure));
    }
}
InitValues::InitValues(double cir_left, int _M, double _x0, double _A, double _omega)
: h(2.0 / _M), Mx(_M), x0(_x0), omega(_omega), A(_A) {

    rho_minus = 1;
    rho_plus = 2;
    c_minus = 1;
    c_plus = 2;
    k_minus = c_minus * c_minus * rho_minus;
    k_plus = c_plus * c_plus * rho_plus;

    z_l = (rho_minus * c_minus);
    z_r = (rho_plus * c_plus);
    A_t = 2 * z_r / (z_l + z_r);
    A_r = (z_r - z_l) / (z_l + z_r);
    // A_r = (z_r - z_l) / (z_l + z_r);

    std::cout << "A_t = " << A_t << std::endl;
    std::cout << "A_r = " << A_r << std::endl;

    // double cir_left = 0.4;
    tau = h * cir_left / c_minus;
    std::cout << "TAU: " << tau << std::endl;

    // double cir_right = 0.6;

    std::cout << "CIR LEFT: " << c_minus * tau / h << std::endl;
    std::cout << "CIR RIGHT: " << c_plus * tau / h << std::endl;

    A_minus << 0, 1 / rho_minus,
        k_minus, 0;
    A_plus << 0, 1 / rho_plus,
        k_plus, 0;

    J = int(Mx / 2);
    alpha = J * h + h / 2;
    // alpha = (sqrt(2) - 1);
    // for (int i = 0; i < _M; ++i) {
    //     if (i * h >= alpha) {
    //         J = i;
    //         break;
    //     }
    // }

    SetInitRadU(x0, omega);
    PrintInit();
}
