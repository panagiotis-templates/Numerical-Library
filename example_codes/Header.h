#pragma once
#include<cmath>
#include<vector>
#include<iostream>
#include<cassert>
#include<concepts>
template<std::floating_point _Floating>
inline[[nodiscard]] _Floating  f(_Floating x)
{
    return exp(2 * x);
}
template<std::floating_point _Floating>
inline [[nodiscard]] _Floating ddf(_Floating x)
{
    return 4 * exp(2 * x);
}
template<std::floating_point _Floating>
inline void finite_difference(_Floating a, _Floating b, _Floating h) {
    _Floating xi{}, res{}, hsq{};
    size_t n{ static_cast<size_t>((b - a) / h + 1) };
    xi = a;
    hsq= std::pow(h, 2);
    for (size_t i = 0; i < n; i++) {
        res = (f(xi - h) - 2 * f(xi) + f(xi + h)) / hsq;
        std::cout << xi << " " << res << " " << ddf(xi) << " " << fabs(ddf(xi) - res) << '\n';
        xi +=  h;
    }
    
}
 
