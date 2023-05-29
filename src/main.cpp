#include <system/system.h>
#include <fstream>
#include <cstdlib>
#include <sstream>

using std::ofstream;

int main() {
    // CIR_left, Mx, x0, A, omega
    // double tau = 0.0004;
    double cir_left = 0.6;
    int N = 1000;
    int Mx = 100;
    double x0 = 0;
    double A = 1;
    double omega = 5;
    System mesh(cir_left, Mx, x0, A, omega);
    mesh.sample();

    // for (int n = 0; n < N; ++n) {
    //     mesh.solve(n * mesh.tau);
    //     // Формирование имени файла с индексом временного шага
    //     std::ostringstream filename;
    //     filename << "../bin/animation/velocity_out_" << std::setfill('0') << std::setw(5) << n << ".bin";
    //     std::ofstream file_velocity(filename.str(), std::ios::binary);
    //     // Запись количества точек данных для текущего временного шага (4 байта, little-endian)
    //     for (int i = 0; i < Mx; ++i) {
    //         float x_coord = mesh.h * i;
    //         float pressure_value = mesh.GetValue(i)(1);
    //         // Запись x, y и p(x, y) (каждый параметр - 4 байта, little-endian)
    //         file_velocity.write(reinterpret_cast<char*>(&x_coord), sizeof(x_coord));
    //         file_velocity.write(reinterpret_cast<char*>(&pressure_value), sizeof(pressure_value));
    //     }
    //     file_velocity.close();
    // }
    return 0;
}
