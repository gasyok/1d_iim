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
void InitValues::SetInitRadU(double x0, double omega) {
    double _pressure;
    for (int i = 0; i < Mx; ++i) {
        double xi = i * h - x0;
        if (xi <= 1 / omega && xi >= 0) {
            _pressure = 0.5 * (1 - cos(2 * M_PI * omega * xi));
        }
        else {
            _pressure = 0;
        }
        u.push_back(Vector2d(1 / (c_minus * rho_minus) * _pressure, _pressure));
    }
}
InitValues::InitValues() : InitValues(0.0004, 0.01, 100, 0.5, 1, 10) {}
InitValues::InitValues(double _tau, double _h, int _M, double _x0, double _A, double _omega)
: tau(_tau), h(_h), Mx(_M), x0(_x0), omega(_omega), A(_A) {

    rho_minus = 1.;
    rho_plus = 1.5;
    c_minus = 1.;
    c_plus = 1.5;
    k_minus = c_minus * c_minus * rho_minus;
    k_plus = c_plus * c_plus * rho_plus;

    A_minus << 0, 1 / rho_minus,
        k_minus, 0;
    A_plus << 0, 1 / rho_plus,
        k_plus, 0;

    J = int(Mx / 3);
    alpha = J * h + h / 2;

    SetInitRadU(x0, omega);
    PrintInit();
}
