#include"Header4.h"


int main()
{
     double a = 0, b = 1, dx = 0.1;
     if(auto u=trapezoid_integral(a, b, dx, [](const double& x) {return std::exp(2 * x); });u.has_value()){
        std::cout<<u.value();
     }
}
