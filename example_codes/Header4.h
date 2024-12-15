#pragma once
#include <iostream>
#include <cmath>
#include<concepts>
#include<cassert>
template<std::floating_point _Floating>
[[nodiscard]] inline _Floating f(const _Floating& x) noexcept
{
    //{for pi  #include<numbers> => std::numbers::pi}
    //return pow(2,x) + 5 * x +3;
    return exp(2 * x);
    //return sin(2 * M_PI *x);
}
template<std::floating_point _Floating>
[[nodiscard]] inline _Floating trapezoid_integral(const _Floating& a,const  _Floating& b,const  _Floating& dx) {
    assert(b > a);
    assert(dx > 0);
    const size_t &iterations = static_cast<size_t>((b - a) / dx);

    _Floating result = 0;
    _Floating xi = a;
    for (size_t i = 0; i < iterations; i++) {
        result += static_cast<_Floating>((f(xi + dx) + f(xi)) * (dx * 0.5));
        std::cout << xi << "-" << xi + dx << " " << result << '\n';
        xi += dx;
    }

    return result;
}






