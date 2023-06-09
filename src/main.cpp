#include <system/system.h>
#include <fstream>
#include <cstdlib>
#include <sstream>

using std::ofstream;

int main() {
    // CIR_left, Mx, x0, A, omega
    // double tau = 0.0004;
    double cir_left = 0.3;
    int Mx = 100;
    double x0 = 0.2;
    double A = 1;
    double omega = 0.3;
    // System mesh(cir_left, Mx, x0, A, omega);
    // mesh.sample();
    ofstream file;
    file.open("../bin/result.data");
    double prev_l1 = 0, prev_linf = 0, l1 = 0, linf = 0;
    for (int k = 0; k < 5; ++k) {
        System mesh(cir_left, Mx, x0, A, omega);
        mesh.sample();
        l1 = mesh.l1;
        linf = mesh.lmax;
        double p1 = 0, pinf = 0;
        if (k > 0) {
            p1 = (log(l1) - log(prev_l1)) / log(0.5);
            pinf = (log(linf) - log(prev_linf)) / log(0.5);
        }
        // double h = mesh.h;
        file << Mx << "    " << l1 << "    " << p1 << "    " << linf << "    " << pinf << std::endl;
        prev_l1 = l1;
        prev_linf = linf;
        Mx *= 2;
    }

    file.close();

    return 0;
}
