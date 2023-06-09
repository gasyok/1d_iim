#include <system/system.h>
#include <fstream>

Vector2d System::GetValue(int i) {
    return u[i];
}
Vector2d System::equation(int i, bool flag) {
    Matrix2d _A = A_minus;

    if (flag) {
        _A = A_plus;
    }
    return u[i] - 0.5 * (tau / h) * _A * (u[i + 1] - u[i - 1]) + 0.5 * (tau / h) * (tau / h) * _A * _A * (u[i + 1] - 2 * u[i] + u[i - 1]);
}
Vector2d System::irrEquation(int i, bool flag) {
    Vector2d res (0, 0);
    const vector<int> offsets = {i - 1, i, i + 1};
    for (int l = 0; l < 3; ++l) {
        if (flag) {
            res += gamma_plus[l] * u[offsets[2 - l]];
        }
        else {
            res += gamma_minus[l] * u[offsets[l]];
        }
    }
    return u[i] + (tau / h) * res; 
}
void System::solve(double t) {
    vector<Vector2d> new_u = u;
    new_u[0] = GetExactSol(0, x0, omega, t);
    new_u[Mx - 1] = GetExactSol(Mx - 1, x0, omega, t);
    for (int i = 1; i < Mx - 1; ++i) {
        if (i == J) {
            new_u[i] = irrEquation(i, false);
            new_u[i + 1] = irrEquation(i + 1, true);
            i = i + 1;
        }
        else {
            if (i > J + 1) {
                new_u[i] = equation(i, true);
            }
            else {
                new_u[i] = equation(i, false);
            }
        }
    }
    u = new_u;
}
double System::l_1(double t) {
    double res = 0;
    for (int i = J + 1; i < Mx - 1; ++i) {
        res += abs(GetValue(i)(1) - GetExactSol(i, x0, omega,  t)(1));
    }
    return res / (Mx - J - 2);
}
double System::l_inf(double t) {
    double res = 0.0;
    for (int i = J + 1; i < Mx - 1; ++i) {
         double temp = abs(GetValue(i)(1) - GetExactSol(i, x0, omega, t)(1));
         if (temp >= res) res = temp;
    }
    return res;
}
void System::sample() {
    int n = 0;
    while (n < total_steps) {
        std::ostringstream filename;
        std::ostringstream exact;
        filename << "../bin/animation/velocity_out_" << std::setfill('0') << std::setw(5) << n << ".bin";
        exact << "../bin/animation/exact_out_" << std::setfill('0') << std::setw(5) << n << ".bin";

        std::ofstream file_velocity(filename.str(), std::ios::binary);
        std::ofstream file_exact(exact.str(), std::ios::binary);
        // Запись количества точек данных для текущего временного шага (4 байта, little-endian)
        for (int i = 0; i < Mx; ++i) {
            float x_coord = h * i;
            float pressure_value = GetValue(i)(1);
            float exact_pressure = GetExactSol(i, x0, omega, n * tau)(1);
            // Запись x, y и p(x, y) (каждый параметр - 4 байта, little-endian)
            file_velocity.write(reinterpret_cast<char*>(&x_coord), sizeof(x_coord));
            file_velocity.write(reinterpret_cast<char*>(&pressure_value), sizeof(pressure_value));

            file_exact.write(reinterpret_cast<char*>(&x_coord), sizeof(x_coord));
            file_exact.write(reinterpret_cast<char*>(&exact_pressure), sizeof(exact_pressure));
        }
        file_velocity.close();
        file_exact.close();
        n++;
        solve(n * tau);
    }
    l1 = l_1(n * tau);
    lmax = l_inf(n * tau);
    std::cout << "--------Compute ready!----------\n";
}
System::System(double _tau, int _M, double _x0, double _A, double _omega)
: PreProcess(_tau, _M, _x0, _A, _omega) {
    // total_time = (alpha - _x0) / c_minus + 0.5 * (2 - alpha) / c_plus;
    total_time = (alpha - _x0) / c_minus + 0.5 / c_plus;
    total_steps = total_time / tau;
    std::cout << "Total Steps: " << total_steps << "\nTotal Time: " << total_time << std::endl;
}
