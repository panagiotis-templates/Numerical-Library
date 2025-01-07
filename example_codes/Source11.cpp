#include"Header11.h"

int main()
{
    finite_diff_central(0.0, 1.0, 0.1, 
        [](double t,double q) {return ( (-0.25 * std::exp(-2 * t)) * (-(std::pow(std::numbers::pi, 2) - 16) * std::sin((std::numbers::pi / 2) * t) - 8 * std::numbers::pi * std::cos((std::numbers::pi / 2) * t)) + q * std::exp(-2 * t) * std::sin((std::numbers::pi / 2) * t)); },
        [](double t) { return ( std::exp(-2 * t) * std::sin((std::numbers::pi / 2) * t)); });
    
    
    return 0;
}
