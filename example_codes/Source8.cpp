#include"Header.h"
#include <iomanip>
#include <chrono>
#include <random>
#include<time.h>




int main()
{
    using namespace panagiotis;
    clock_t begin = std::clock();
    double a = -1.0, b = 1.0, h = 0.1, step_size = (b - a) / h + 1, step = a;
    auto f = [](double x)-> double { return std::exp(std::exp(2 * x)); };
    size_t n = 12;
    std::vector<std::vector<   double>> A(n, std::vector<   double>(n, 0));
    std::vector<   double> F(n, 0);
    //std::vector< double> C(n, 0);

    // Build left hand side
    long double temp = 0;
    for (unsigned int i = 0; i < n; i++) {
        A[i][i] = gauss_legendre(i, i, a, b);
        for (size_t j = 0; j < n; j++) {

        }
        F[i] = gauss_legendre_f(i, a, b, f);
    }
    auto op = backwards_substitution(A, F);
    if (!op.has_value()) {
        return 1;
    }
    std::vector<  double>C = std::move(op.value());

    /*  print2DVector(A);
     print_Vector(F);  */

    step = a;
    step_size = (b - a) / h + 1;
    std::cout << "x  Aproximation    Real " << '\n';
    for (size_t o = 0; o < step_size; o++) {
        std::cout << step << " " << Pnf(step, C) << " " << f(step) << " " << '\n';
        step += h;
    }


    clock_t end = std::clock();
    long double time_spent = (long double)(end - begin) / CLOCKS_PER_SEC;
    std::cout << std::fixed;

    std::cout << std::setprecision(10) << "L2 norm: " << error_gauss(f, Pnf<  double>, a, b + h, C) << " time: " << time_spent << '\n';
    //print2DVector(A);
    print_Vector(F);
    print_Vector(C);
    return 0;
}
