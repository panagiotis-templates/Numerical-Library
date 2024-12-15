#pragma once
#include <iostream>
#include <cmath>
#include<cassert>
#include<type_traits>
template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
[[nodiscard]] inline _Ty f(const _Ty& x) noexcept
{
    //{for pi  #include<numbers> => std::numbers::pi}
    //return pow(2,x) + 5 * x +3;
    return std::exp(2 * x);
    //return sin(2 * M_PI *x);
}
template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
[[nodiscard]] inline _Ty trapezoid_integral(const _Ty& a,const  _Ty& b,const  _Ty& dx) {
    assert(b > a);
    assert(dx > 0);
    const size_t &iterations = static_cast<size_t>((b - a) / dx);

    _Ty result = 0;
    _Ty xi = a;
    for (size_t i = 0; i < iterations; i++) {
        result += static_cast<_Ty>((f(xi + dx) + f(xi)) * (dx * 0.5));
        std::cout << xi << "-" << xi + dx << " " << result << '\n';
        xi += dx;
    }

    return result;
}







