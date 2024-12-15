#pragma once
#include<concepts>
#include<iostream>
#include<cmath>
//#include<concepts>
template<std::floating_point _Floating>
inline [[nodiscard]] _Floating  f(const _Floating& x) noexcept
{
    return std::exp(2 * x);
}
template<std::floating_point _Floating>
inline [[nodiscard]] _Floating ddf(const _Floating& x) noexcept
{
    return 4 * std::exp(2 * x);
}
template<std::floating_point _Floating>
inline void finite_difference(const _Floating& a, const  _Floating& b, const  _Floating& h) {
   
    _Floating xi{}, res{}, hsq{ std::pow(h, 2) };
    const size_t& n{ static_cast<size_t>((b - a) / h + 1) };
    xi = a;
    for (size_t i = 0; i < n; i++) {
        res = (f(xi - h) - 2 * f(xi) + f(xi + h)) / hsq;
        std::cout << xi << " " << res << " " << ddf(xi) << " " << std::abs(ddf(xi) - res) << '\n';
        xi += h;
    }

}
