#include"Header4.h"


int main()
{
     double a = 0, b = 1, dx = 0.1;
 std::cout << trapezoid_integral(a, b, dx, [](const double& x) {return std::exp(2 * x); });

 return 0;
}
